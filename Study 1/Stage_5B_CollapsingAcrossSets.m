%% ===== Stage 5B: Condition-level Aggregation (across SETS; files in different folders) =====
% INPUT:
%   Select N setAgg files one-by-one (each can be in a different folder):
%       *_stage5A_setAgg.mat
%   They must belong to ONE participant + ONE condition.
%
% OUTPUT:
%   *_stage5B_condAgg.mat  (saved into a user-chosen output folder)
%
% Stage 5B meaning:
%   For each bin (early / late / all), average the Stage5A outputs across sets.

clear; close all;

%%  (1) User parameters
expected_nSets = 2;              % <-- SET THIS (e.g., 2 or 3)
collapse_method_sets = "mean";   % "mean" or "median"
doQCplot = true;

bins = {'early','late','all'};

reqCurveFields = {'torque_101','angle_101','emgMeanNorm_101','centroidXc_101','centroidYc_101'};
reqBipField    = 'emgRMS_bipNorm_101';
reqMapField    = 'emgRMS_bipNorm_map28';

%%  (2) Load Stage5A setAgg files one-by-one (from different folders)
sets = cell(expected_nSets,1); % Creates a cell array to store the loaded setAgg structures from each set
sourceFiles = strings(expected_nSets,1); %Creates a string array to store the file paths of the sets you load

for s = 1:expected_nSets
    [f, p] = uigetfile('*_stage5A_setAgg.mat', sprintf('Select Stage5A setAgg file (%d/%d)', s, expected_nSets));
    if isequal(f,0)
        error('User cancelled at set %d/%d. Need exactly %d sets.', s, expected_nSets, expected_nSets);
    end

    fullp = fullfile(p,f);
    L = load(fullp);
    assert(isfield(L,'setAgg'), 'File missing variable "setAgg": %s', fullp);

    sets{s} = L.setAgg; %Store the loaded set in the cell array
    sourceFiles(s) = string(fullp); % Record which file was used
end

fprintf('✅ Loaded %d sets.\n', expected_nSets); % ✅ Loaded 3 sets

% Template
tmpl = sets{1}; %Select the first set as a template
%The template is used to determine: number of time points /field structure
%/coordinate definitions, as all sets share the same struct

% Determine nT
assert(isfield(tmpl,'params') && isfield(tmpl.params,'nTime'), 'Template setAgg missing params.nTime');
nT = tmpl.params.nTime;
xAxis = linspace(0,100,nT);


%%  (3) Consistency checks
% Check that all sets have the same nTime
for s = 1:expected_nSets
    if ~isfield(sets{s},'params') || ~isfield(sets{s}.params,'nTime')
        error('Set %d missing params.nTime.', s);
    end
    if sets{s}.params.nTime ~= nT
        error('nTime mismatch: set %d has %d but template has %d.', s, sets{s}.params.nTime, nT);
    end
end

% Check same participantName / conditionName if present
if isfield(tmpl,'meta')
    if isfield(tmpl.meta,'participantName')
        p0 = string(tmpl.meta.participantName);
    else
        p0 = "";
    end

    if isfield(tmpl.meta,'conditionName')
        c0 = string(tmpl.meta.conditionName);
    else
        c0 = "";
    end

    for s = 2:expected_nSets
        if isfield(sets{s},'meta')
            if strlength(p0) > 0 && isfield(sets{s}.meta,'participantName')
                if string(sets{s}.meta.participantName) ~= p0
                    error('participantName mismatch between set1 and set%d.', s);
                end
            end

            if strlength(c0) > 0 && isfield(sets{s}.meta,'conditionName')
                if string(sets{s}.meta.conditionName) ~= c0
                    error('conditionName mismatch between set1 and set%d.', s);
                end
            end
        end
    end
end

%%  (4) Collapse helpers

collapseAcrossSets_1x101 = @(curvesCell) local_collapse(cat(1, curvesCell{:}), collapse_method_sets, 1); % -> 1x101, torque(angle ,EMG mean ,centroid X ,centroid Y)
collapseAcrossSets_28x101 = @(matsCell)  local_collapse(cat(3, matsCell{:}),  collapse_method_sets, 3); % -> 28x101 (Based on Each set contributes one spatial-temporal EMG matrix)
collapseAcrossSets_28x1  = @(vecsCell)   local_collapse(cat(2, vecsCell{:}),  collapse_method_sets, 2); % -> 28x1 (This is used for the time-averaged channel vector)

%%  (5) Build condAgg structure
condAgg = struct();

condAgg.params = struct( ...
    'expected_nSets', expected_nSets, ...
    'collapse_method_sets', string(collapse_method_sets), ...
    'nTime', nT, ...
    'time_axis_label', '%CV (0–100)', ...
    'time_axis_values', xAxis, ...
    'definition', 'Stage5B collapses Stage5A outputs across sets into condition-level early/late/all curves' );

condAgg.sourceFiles = sourceFiles;

% Keep mapping + coords (assume consistent)
if isfield(tmpl,'meta'), condAgg.meta = tmpl.meta; else, condAgg.meta = struct(); end
if isfield(tmpl,'mapInfo'), condAgg.mapInfo = tmpl.mapInfo; else, condAgg.mapInfo = struct(); end
if isfield(tmpl,'coords'), condAgg.coords = tmpl.coords; else, condAgg.coords = struct(); end
if isfield(tmpl,'centroid_params'), condAgg.centroid_params = tmpl.centroid_params; else, condAgg.centroid_params = struct(); end

%%  (6) Aggregate per bin (early/late/all)

for b = 1:numel(bins)
    bin = bins{b};

    % ---- CURVES (1x101) ----
    for f = 1:numel(reqCurveFields)
        fn = reqCurveFields{f};
        C = local_gather_1x101(sets, bin, fn, nT);  % cell of 1x101

        if isempty(C)
            condAgg.curve.(bin).(fn) = [];
        else
            condAgg.curve.(bin).(fn) = collapseAcrossSets_1x101(C);
        end

        condAgg.nSetsUsed.(bin).(fn) = numel(C);
    end

    % ---- BIPOLAR (28x101) ----
    B = local_gather_28x101(sets, bin, reqBipField, nT);
    if isempty(B)
        condAgg.bipolar.(bin).(reqBipField) = [];
    else
        condAgg.bipolar.(bin).(reqBipField) = collapseAcrossSets_28x101(B);
    end
    condAgg.nSetsUsed.(bin).bipolar = numel(B);

    % ---- MAP (28x1) ----
    V = local_gather_28x1(sets, bin, reqMapField);
    if isempty(V)
        condAgg.map.(bin).(reqMapField) = [];
    else
        condAgg.map.(bin).(reqMapField) = collapseAcrossSets_28x1(V);
    end
    condAgg.nSetsUsed.(bin).map28 = numel(V);
end

%%  (7) QC plot: Early vs Late (condition-level)

if doQCplot
    figure('Color','w','Name','Stage5B QC: Condition-level Early vs Late');
    tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    % Torque
    nexttile;
    lg = {};
    if ~isempty(condAgg.curve.early.torque_101)
        plot(xAxis, condAgg.curve.early.torque_101,'LineWidth',2); hold on;
        lg{end+1} = 'Early';
    end
    if ~isempty(condAgg.curve.late.torque_101)
        plot(xAxis, condAgg.curve.late.torque_101,'LineWidth',2);
        lg{end+1} = 'Late';
    end
    title('Torque'); xlabel('%CV'); ylabel('Torque'); box off;
    if ~isempty(lg), legend(lg,'Box','off'); end

    % Angle
    nexttile;
    if ~isempty(condAgg.curve.early.angle_101), plot(xAxis, condAgg.curve.early.angle_101,'LineWidth',2); hold on; end
    if ~isempty(condAgg.curve.late.angle_101),  plot(xAxis, condAgg.curve.late.angle_101,'LineWidth',2); end
    title('Angle'); xlabel('%CV'); ylabel('Angle'); box off;

    % EMG mean (%MVC)
    nexttile;
    if ~isempty(condAgg.curve.early.emgMeanNorm_101), plot(xAxis, condAgg.curve.early.emgMeanNorm_101,'LineWidth',2); hold on; end
    if ~isempty(condAgg.curve.late.emgMeanNorm_101),  plot(xAxis, condAgg.curve.late.emgMeanNorm_101,'LineWidth',2); end
    title('EMG mean (%MVC)'); xlabel('%CV'); ylabel('%MVC'); box off;

    % Centroid X centered
    nexttile;
    if ~isempty(condAgg.curve.early.centroidXc_101), plot(xAxis, condAgg.curve.early.centroidXc_101,'LineWidth',2); hold on; end
    yline(0,'--');
    if ~isempty(condAgg.curve.late.centroidXc_101),  plot(xAxis, condAgg.curve.late.centroidXc_101,'LineWidth',2); end
    title('Centroid X (centered, mm)'); xlabel('%CV'); ylabel('mm'); box off;

    % Centroid Y centered
    nexttile;
    if ~isempty(condAgg.curve.early.centroidYc_101), plot(xAxis, condAgg.curve.early.centroidYc_101,'LineWidth',2); hold on; end
    yline(0,'--');
    if ~isempty(condAgg.curve.late.centroidYc_101),  plot(xAxis, condAgg.curve.late.centroidYc_101,'LineWidth',2); end
    title('Centroid Y (centered, mm)'); xlabel('%CV'); ylabel('mm'); box off;

    nexttile; axis off;

    tt = sprintf('Stage5B QC | %d sets | collapse=%s', expected_nSets, string(collapse_method_sets));
    if isfield(condAgg,'meta') && isfield(condAgg.meta,'participantName') && isfield(condAgg.meta,'conditionName')
        tt = sprintf('%s | %s | Stage5B QC (%d sets, %s)', ...
            string(condAgg.meta.participantName), string(condAgg.meta.conditionName), expected_nSets, string(collapse_method_sets));
    end
    sgtitle(tt, 'FontWeight','bold','Interpreter','none');
end

%%  (8) Save Stage 5B output (choose output folder)

outFolder = uigetdir(pwd, 'Select output folder for Stage5B condAgg');
if isequal(outFolder,0), outFolder = pwd; end

baseTag = "condAgg";
if isfield(condAgg,'meta') && isfield(condAgg.meta,'participantName') && isfield(condAgg.meta,'conditionName')
    baseTag = sprintf('%s_%s', string(condAgg.meta.participantName), string(condAgg.meta.conditionName));
end

baseTag = regexprep(string(baseTag), '[^\w-]', '_');

outName = fullfile(outFolder, sprintf('%s_stage5B_condAgg.mat', baseTag));
save(outName, 'condAgg', '-v7.3');
fprintf('✅ Saved Stage 5B condition aggregation:\n%s\n', outName);




%% =========================
%  Local functions
%% =========================
function Y = local_collapse(X, method, dim)
    if isempty(X)
        Y = [];
        return
    end
    switch lower(string(method))
        case "mean"
            Y = mean(X, dim, 'omitnan');
        case "median"
            Y = median(X, dim, 'omitnan');
        otherwise
            error('Unknown collapse method: %s', method);
    end
end

function C = local_gather_1x101(sets, binName, fieldName, nT)
    C = {};
    for s = 1:numel(sets)
        A = sets{s};
        if ~isfield(A,'curve'), continue; end
        if ~isfield(A.curve, binName), continue; end
        if ~isfield(A.curve.(binName), fieldName), continue; end
        v = A.curve.(binName).(fieldName);
        if isempty(v), continue; end
        v = reshape(v, 1, []);
        if numel(v) ~= nT, continue; end
        C{end+1} = v; %#ok<AGROW>
    end
end

function C = local_gather_28x101(sets, binName, fieldName, nT)
    C = {};
    for s = 1:numel(sets)
        A = sets{s};
        if ~isfield(A,'bipolar'), continue; end
        if ~isfield(A.bipolar, binName), continue; end
        if ~isfield(A.bipolar.(binName), fieldName), continue; end
        X = A.bipolar.(binName).(fieldName);
        if isempty(X), continue; end
        if ~ismatrix(X) || size(X,1)~=28 || size(X,2)~=nT, continue; end
        C{end+1} = X; %#ok<AGROW>
    end
end

function C = local_gather_28x1(sets, binName, fieldName)
    C = {};
    for s = 1:numel(sets)
        A = sets{s};
        if ~isfield(A,'map'), continue; end
        if ~isfield(A.map, binName), continue; end
        if ~isfield(A.map.(binName), fieldName), continue; end
        v = A.map.(binName).(fieldName);
        if isempty(v), continue; end
        v = v(:);
        if numel(v) ~= 28, continue; end
        C{end+1} = v; %#ok<AGROW>
    end
end












% %% ===== Stage 5B STRUCT CHECK =====
% 
% fprintf('\n===== Stage 5B STRUCT CHECK =====\n')
% 
% % Basic info
% fprintf('Expected sets: %d\n', expected_nSets);
% 
% if isfield(condAgg,'meta') && isfield(condAgg.meta,'participantName')
%     fprintf('Participant: %s\n', string(condAgg.meta.participantName));
% end
% if isfield(condAgg,'meta') && isfield(condAgg.meta,'conditionName')
%     fprintf('Condition: %s\n', string(condAgg.meta.conditionName));
% end
% 
% if isfield(condAgg,'params') && isfield(condAgg.params,'collapse_method_sets')
%     fprintf('Collapse method across sets: %s\n', string(condAgg.params.collapse_method_sets));
% end
% 
% % Determine nTime safely
% if isfield(condAgg,'curve') && isfield(condAgg.curve,'all') && isfield(condAgg.curve.all,'torque_101') ...
%         && ~isempty(condAgg.curve.all.torque_101)
%     nT_check = numel(condAgg.curve.all.torque_101);
% else
%     nT_check = nT;
% end
% fprintf('nTime = %d\n', nT_check);
% 
% binsToCheck = {'all','early','late'};
% 
% for b = 1:numel(binsToCheck)
%     bin = binsToCheck{b};
% 
%     fprintf('\n--- BIN: %s ---\n', upper(bin));
% 
%     % nSetsUsed
%     if isfield(condAgg,'nSetsUsed') && isfield(condAgg.nSetsUsed, bin)
%         disp(condAgg.nSetsUsed.(bin))
%     else
%         fprintf('nSetsUsed.%s missing.\n', bin);
%     end
% 
%     % --- Curves (1x101 expected) ---
%     if isfield(condAgg,'curve') && isfield(condAgg.curve, bin)
%         curveFields = {'torque_101','angle_101','emgMeanNorm_101','centroidXc_101','centroidYc_101'};
% 
%         for f = 1:numel(curveFields)
%             fn = curveFields{f};
% 
%             if isfield(condAgg.curve.(bin), fn) && ~isempty(condAgg.curve.(bin).(fn))
%                 v = condAgg.curve.(bin).(fn);
% 
%                 fprintf('%s size: ', fn);
%                 disp(size(v))
% 
%                 fprintf('%s numel: %d\n', fn, numel(v));
% 
%                 if numel(v) ~= nT_check
%                     warning('curve.%s.%s does not have %d elements.', bin, fn, nT_check);
%                 end
% 
%                 fprintf('%s range: %.2f → %.2f\n', fn, ...
%                     min(v,[],'omitnan'), max(v,[],'omitnan'));
%             else
%                 fprintf('%s is missing/empty.\n', fn);
%             end
%         end
%     else
%         fprintf('curve.%s missing.\n', bin);
%     end
% 
%     % --- Bipolar map over time (28x101 expected) ---
%     if isfield(condAgg,'bipolar') && isfield(condAgg.bipolar, bin)
%         fn = 'emgRMS_bipNorm_101';
% 
%         if isfield(condAgg.bipolar.(bin), fn) && ~isempty(condAgg.bipolar.(bin).(fn))
%             X = condAgg.bipolar.(bin).(fn);
% 
%             fprintf('%s size: ', fn);
%             disp(size(X))
% 
%             if ~all(size(X) == [28 nT_check])
%                 warning('bipolar.%s.%s is not 28x%d.', bin, fn, nT_check);
%             end
% 
%             fprintf('%s range: %.2f → %.2f\n', fn, ...
%                 min(X(:),[],'omitnan'), max(X(:),[],'omitnan'));
%         else
%             fprintf('%s is missing/empty.\n', fn);
%         end
%     else
%         fprintf('bipolar.%s missing.\n', bin);
%     end
% 
%     % --- Time-averaged channel vector (28x1 expected) ---
%     if isfield(condAgg,'map') && isfield(condAgg.map, bin)
%         fn = 'emgRMS_bipNorm_map28';
% 
%         if isfield(condAgg.map.(bin), fn) && ~isempty(condAgg.map.(bin).(fn))
%             v = condAgg.map.(bin).(fn);
% 
%             fprintf('%s size: ', fn);
%             disp(size(v))
% 
%             fprintf('%s numel: %d\n', fn, numel(v));
% 
%             if numel(v) ~= 28
%                 warning('map.%s.%s does not have 28 elements.', bin, fn);
%             end
% 
%             fprintf('%s range: %.2f → %.2f\n', fn, ...
%                 min(v,[],'omitnan'), max(v,[],'omitnan'));
%         else
%             fprintf('%s is missing/empty.\n', fn);
%         end
%     else
%         fprintf('map.%s missing.\n', bin);
%     end
% end
% 
% 
% 
% fprintf('\n===== Stage 5B SET-USAGE CHECK =====\n')
% 
% binsToCheck = {'all','early','late'};
% 
% for b = 1:numel(binsToCheck)
%     bin = binsToCheck{b};
% 
%     fprintf('\n--- %s ---\n', upper(bin));
% 
%     if isfield(condAgg,'nSetsUsed') && isfield(condAgg.nSetsUsed, bin)
%         S = condAgg.nSetsUsed.(bin);
%         disp(S)
% 
%         vals = struct2cell(S);
%         vals = cellfun(@double, vals);
% 
%         fprintf('Min sets used: %d\n', min(vals));
%         fprintf('Max sets used: %d\n', max(vals));
% 
%         if any(vals < expected_nSets)
%             warning('%s bin used fewer than expected sets in at least one field.', bin);
%         end
%     else
%         warning('nSetsUsed.%s missing.', bin);
%     end
% end