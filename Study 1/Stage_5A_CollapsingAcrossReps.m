%% ===== Stage 5A: Set-level Aggregation (within ONE set file) =====
% INPUT :
%   (1) *_stage3_EMGnorm.mat      (rep-level: torque_101, angle_101, EMG, bipolar 28x101)
%   (2) *_stage4_centroid_mm.mat  (rep-level: centroidX/Y (absolute + centered), per-rep 101-pt)
%
% OUTPUT:
%   *_stage5A_setAgg.mat  (set-level: one trace per bin: early / late / all)
%
% WHAT "collapse reps" MEANS HERE:
%   - EARLY  = average across first nEarly VALID reps  -> 1x101 (and 28x101 for bipolar)
%   - LATE   = average across last  nLate  VALID reps  -> 1x101 (and 28x101 for bipolar)
%   - ALL    = average across ALL VALID reps           -> 1x101 (and 28x101 for bipolar)
%
% Notes:
%   - Stage 5A NEVER creates new rep categories. It only aggregates rep-level outputs into set-level curves.
%   - Robust to whether rep curves are stored as 101x1 or 1x101 (we force row vectors before stacking).

clear;


%%  (1) User parameters

nEarly = 2;                 % number of early reps to average
nLate  = 2;                 % number of late reps to average
collapse_method = "mean";   % "mean" or "median" (rep->set)
modeIfTooFew = "error";     % "error" | "shrink" | "skipLate"
doQCplot = true;

%%  (2) Load Stage 3 + Stage 4 (same set)
[file3, path3] = uigetfile('*_stage3_EMGnorm.mat', 'Select Stage-3 EMGnorm (THIS set)');
if isequal(file3,0), return; end
S3 = load(fullfile(path3, file3));
assert(isfield(S3,'stage3'), 'Stage-3 file does not contain variable "stage3".');
stage3 = S3.stage3;

[file4, path4] = uigetfile('*_stage4_centroid_mm.mat', 'Select Stage-4 centroid (SAME set)');
if isequal(file4,0), return; end
S4 = load(fullfile(path4, file4));
assert(isfield(S4,'stage4'), 'Stage-4 file does not contain variable "stage4".');
stage4 = S4.stage4;

% Consistency
nReps3 = numel(stage3.rep);
nReps4 = numel(stage4.rep);
assert(nReps3 == nReps4, 'Stage3 and Stage4 rep counts differ (%d vs %d).', nReps3, nReps4); % Only reps valid in both stages are used

nT = stage3.params.timeNorm_points;
xAxis = linspace(0,100,nT);

%%  (3) Determine valid reps and pick early/late
valid3 = arrayfun(@(s) logical(s.is_valid), stage3.rep);
valid4 = arrayfun(@(s) logical(s.is_valid), stage4.rep);
validMask = valid3 & valid4; % both steps must be consistent at this line

repValid = find(validMask);
if isempty(repValid)
    error('No valid reps in this set (validMask is empty).');
end

[repEarly, repLate] = pickEarlyLate(repValid, nEarly, nLate, modeIfTooFew);

fprintf('Valid reps: %s | Early: %s | Late: %s\n', mat2str(repValid(:)), mat2str(repEarly), mat2str(repLate));


%%  (4) Robust stacking helpers (force 1x101 rows)
asRow = @(v) reshape(v, 1, []);  % works for 101x1 and 1x101

getMat   = @(repStruct, reps, field) local_cat_dim( ...
    arrayfun(@(r) asRow(repStruct(r).(field)), reps, 'UniformOutput', false), 1); %For each rep index, extracts a field, forces it to 1×N, then concatenates along dim 1, so R × 101

getMat28 = @(repStruct, reps, field) local_cat_dim( ...
    arrayfun(@(r) repStruct(r).(field), reps, 'UniformOutput', false), 3); %Concatenates per-rep matrices along dim 3

% Collapse across reps (dim 1 for [R x 101], dim 3 for [28 x 101 x R])
collapseR = @(X) local_collapse(X, collapse_method, 1);  % -> 1 × 101 = one time-normalized curve across the contraction cycle 
collapseB = @(X) local_collapse(X, collapse_method, 3);  % -> 28x101, so reps averaged live in D3



%%  (5) Build setAgg structure
setAgg = struct();

% Keep meta + mapping + coordinate conventions for traceability
if isfield(stage3,'meta')
    setAgg.meta = stage3.meta;
else
    setAgg.meta = struct();
end

setAgg.params = struct( ...
    'nTime', nT, ...
    'nEarly', nEarly, ...
    'nLate', nLate, ...
    'collapse_method', string(collapse_method), ...
    'modeIfTooFew', string(modeIfTooFew), ...
    'time_axis_label', '%CV (0–100)', ...
    'time_axis_values', xAxis, ...
    'definition', 'Stage5A collapses rep-level outputs into set-level early/late/all curves' );

setAgg.reps = struct( ...
    'valid', repValid(:), ...
    'early', repEarly(:), ...
    'late', repLate(:));

% Keep grid mapping & centroid coordinate definitions (from Stage 3/4)
if isfield(stage3,'map')
    setAgg.mapInfo = stage3.map;
else
    setAgg.mapInfo = struct();
end

if isfield(stage4,'coords')
    setAgg.coords = stage4.coords;
else
    setAgg.coords = struct();
end

if isfield(stage4,'params')
    setAgg.centroid_params = stage4.params;
else
    setAgg.centroid_params = struct();
end
%%  (6) Build curves: ALL / EARLY / LATE

% ---- ALL (all valid reps)
TOR_all = getMat(stage3.rep, repValid, 'torque_101');                 % [R x 101]
ANG_all = getMat(stage3.rep, repValid, 'angle_101');
EMG_all = getMat(stage3.rep, repValid, 'emgRMSmeanNorm_101');         % EMG mean (%MVC)
CX_all  = getMat(stage4.rep, repValid, 'centroidX_centered_mm_101');  % centered (0,0=grid center)
CY_all  = getMat(stage4.rep, repValid, 'centroidY_centered_mm_101');

setAgg.curve.all.torque_101        = collapseR(TOR_all);
setAgg.curve.all.angle_101         = collapseR(ANG_all);
setAgg.curve.all.emgMeanNorm_101   = collapseR(EMG_all);
setAgg.curve.all.centroidXc_101    = collapseR(CX_all);
setAgg.curve.all.centroidYc_101    = collapseR(CY_all);

% bipolar 28x101
B_all = getMat28(stage3.rep, repValid, 'emgRMS_bipNorm_101');         % [28 x 101 x R]
setAgg.bipolar.all.emgRMS_bipNorm_101 = collapseB(B_all);             % [28 x 101]
setAgg.map.all.emgRMS_bipNorm_map28   = mean(setAgg.bipolar.all.emgRMS_bipNorm_101, 2, 'omitnan'); % [28 x 1]

% ---- EARLY
TOR_e = getMat(stage3.rep, repEarly, 'torque_101');
ANG_e = getMat(stage3.rep, repEarly, 'angle_101');
EMG_e = getMat(stage3.rep, repEarly, 'emgRMSmeanNorm_101');
CX_e  = getMat(stage4.rep, repEarly, 'centroidX_centered_mm_101');
CY_e  = getMat(stage4.rep, repEarly, 'centroidY_centered_mm_101');

setAgg.curve.early.torque_101      = collapseR(TOR_e);
setAgg.curve.early.angle_101       = collapseR(ANG_e);
setAgg.curve.early.emgMeanNorm_101 = collapseR(EMG_e);
setAgg.curve.early.centroidXc_101  = collapseR(CX_e);
setAgg.curve.early.centroidYc_101  = collapseR(CY_e);

B_e = getMat28(stage3.rep, repEarly, 'emgRMS_bipNorm_101');
setAgg.bipolar.early.emgRMS_bipNorm_101 = collapseB(B_e);
setAgg.map.early.emgRMS_bipNorm_map28   = mean(setAgg.bipolar.early.emgRMS_bipNorm_101, 2, 'omitnan');

% ---- LATE (may be empty if modeIfTooFew="skipLate")
if ~isempty(repLate)
    TOR_l = getMat(stage3.rep, repLate, 'torque_101');
    ANG_l = getMat(stage3.rep, repLate, 'angle_101');
    EMG_l = getMat(stage3.rep, repLate, 'emgRMSmeanNorm_101');
    CX_l  = getMat(stage4.rep, repLate, 'centroidX_centered_mm_101');
    CY_l  = getMat(stage4.rep, repLate, 'centroidY_centered_mm_101');

    setAgg.curve.late.torque_101      = collapseR(TOR_l);
    setAgg.curve.late.angle_101       = collapseR(ANG_l);
    setAgg.curve.late.emgMeanNorm_101 = collapseR(EMG_l);
    setAgg.curve.late.centroidXc_101  = collapseR(CX_l);
    setAgg.curve.late.centroidYc_101  = collapseR(CY_l);

    B_l = getMat28(stage3.rep, repLate, 'emgRMS_bipNorm_101');
    setAgg.bipolar.late.emgRMS_bipNorm_101 = collapseB(B_l);
    setAgg.map.late.emgRMS_bipNorm_map28   = mean(setAgg.bipolar.late.emgRMS_bipNorm_101, 2, 'omitnan');
else
    setAgg.curve.late   = struct();
    setAgg.bipolar.late = struct();
    setAgg.map.late     = struct();
end



%%  (7) Sanity check sizes (must be 1x101 and 28x101)

assert(numel(setAgg.curve.early.torque_101) == nT, 'Early torque curve length is not %d.', nT);
assert(numel(setAgg.curve.early.angle_101)  == nT, 'Early angle curve length is not %d.', nT);
assert(all(size(setAgg.bipolar.early.emgRMS_bipNorm_101) == [28 nT]), 'Early bipolar size is not 28x%d.', nT);

assert(numel(setAgg.curve.all.torque_101) == nT, 'All torque curve length is not %d.', nT);
assert(numel(setAgg.curve.all.angle_101)  == nT, 'All angle curve length is not %d.', nT);
assert(all(size(setAgg.bipolar.all.emgRMS_bipNorm_101) == [28 nT]), 'All bipolar size is not 28x%d.', nT);

if isfield(setAgg.curve,'late') && isfield(setAgg.curve.late,'torque_101') && ~isempty(setAgg.curve.late.torque_101)
    assert(numel(setAgg.curve.late.torque_101) == nT, 'Late torque curve length is not %d.', nT);
    assert(numel(setAgg.curve.late.angle_101)  == nT, 'Late angle curve length is not %d.', nT);
    assert(all(size(setAgg.bipolar.late.emgRMS_bipNorm_101) == [28 nT]), 'Late bipolar size is not 28x%d.', nT);
end

%% (8) QC plot (Torque + Angle + EMG + Centroids)

if doQCplot
    figure('Color','w','Name','Stage5A QC: Early vs Late (THIS set)');
    tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    % Torque (ALL + EARLY + LATE)
    nexttile;
    lg = {};
    if isfield(setAgg.curve,'all') && isfield(setAgg.curve.all,'torque_101') && ~isempty(setAgg.curve.all.torque_101)
        plot(xAxis, setAgg.curve.all.torque_101,'LineWidth',2); hold on;
        lg{end+1} = 'All';
    end
    if isfield(setAgg.curve,'early') && isfield(setAgg.curve.early,'torque_101') && ~isempty(setAgg.curve.early.torque_101)
        plot(xAxis, setAgg.curve.early.torque_101,'LineWidth',2);
        lg{end+1} = 'Early';
    end
    if isfield(setAgg.curve,'late') && isfield(setAgg.curve.late,'torque_101') && ~isempty(setAgg.curve.late.torque_101)
        plot(xAxis, setAgg.curve.late.torque_101,'LineWidth',2);
        lg{end+1} = 'Late';
    end
    title('Torque'); xlabel('%CV'); box off;
    if ~isempty(lg), legend(lg,'Box','off'); end

    % Angle
    nexttile;
    plot(xAxis, setAgg.curve.early.angle_101,'LineWidth',2); hold on;
    if isfield(setAgg.curve,'late') && isfield(setAgg.curve.late,'angle_101') && ~isempty(setAgg.curve.late.angle_101)
        plot(xAxis, setAgg.curve.late.angle_101,'LineWidth',2);
    end
    title('Angle'); xlabel('%CV'); box off;

    % EMG mean (%MVC)
    nexttile;
    plot(xAxis, setAgg.curve.early.emgMeanNorm_101,'LineWidth',2); hold on;
    if isfield(setAgg.curve,'late') && isfield(setAgg.curve.late,'emgMeanNorm_101') && ~isempty(setAgg.curve.late.emgMeanNorm_101)
        plot(xAxis, setAgg.curve.late.emgMeanNorm_101,'LineWidth',2);
    end
    title('EMG mean (%MVC)'); xlabel('%CV'); box off;

    % Centroid X centered (mm)
    nexttile;
    plot(xAxis, setAgg.curve.early.centroidXc_101,'LineWidth',2); hold on;
    yline(0,'--');
    if isfield(setAgg.curve,'late') && isfield(setAgg.curve.late,'centroidXc_101') && ~isempty(setAgg.curve.late.centroidXc_101)
        plot(xAxis, setAgg.curve.late.centroidXc_101,'LineWidth',2);
    end
    title('Centroid X (centered, mm)'); xlabel('%CV'); box off;

    % Centroid Y centered (mm)
    nexttile;
    plot(xAxis, setAgg.curve.early.centroidYc_101,'LineWidth',2); hold on;
    yline(0,'--');
    if isfield(setAgg.curve,'late') && isfield(setAgg.curve.late,'centroidYc_101') && ~isempty(setAgg.curve.late.centroidYc_101)
        plot(xAxis, setAgg.curve.late.centroidYc_101,'LineWidth',2);
    end
    title('Centroid Y (centered, mm)'); xlabel('%CV'); box off;

    nexttile; axis off;

    tt = 'Stage5A (set-level) QC';
    if isfield(setAgg,'meta') && isfield(setAgg.meta,'participantName') && isfield(setAgg.meta,'conditionName')
        tt = sprintf('%s | %s | Stage5A QC', string(setAgg.meta.participantName), string(setAgg.meta.conditionName));
    end
    sgtitle(tt, 'FontWeight','bold', 'Interpreter','none');
end
%%  (9) Save Stage 5A output
[~, baseName, ~] = fileparts(file3);
baseName = strrep(baseName, '_stage3_EMGnorm', '');
outName = fullfile(path3, sprintf('%s_stage5A_setAgg.mat', baseName));
save(outName, 'setAgg', '-v7.3');
fprintf('✅ Saved Stage 5A setAgg:\n%s\n', outName);

%% =========================
%  Local functions
%% =========================
function M = local_cat_dim(C, dim)
    % C is a cell array of numeric matrices to concatenate
    if isempty(C)
        M = [];
        return
    end
    M = cat(dim, C{:});
end




function [earlyReps, lateReps] = pickEarlyLate(validReps, nEarly, nLate, modeIfTooFew)
    validReps = validReps(:).';   % row
    nV = numel(validReps);

    % If not enough reps for early+late
    if nV < (nEarly + nLate)
        switch lower(string(modeIfTooFew))
            case "error"
                error('Not enough valid reps (%d) for nEarly=%d and nLate=%d.', nV, nEarly, nLate);

            case "shrink"
                % shrink proportionally: take as many as possible
                nEarlyUse = min(nEarly, nV);
                remaining = nV - nEarlyUse;
                nLateUse  = min(nLate, remaining);
                earlyReps = validReps(1:nEarlyUse);
                lateReps  = validReps(max(1, nV-nLateUse+1):nV);

            case "skiplate"
                earlyReps = validReps(1:min(nEarly, nV));
                lateReps  = [];

            otherwise
                error('Unknown modeIfTooFew: %s (use "error"|"shrink"|"skipLate")', modeIfTooFew);
        end
    else
        earlyReps = validReps(1:nEarly);
        lateReps  = validReps(end-nLate+1:end);
    end
end

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








% %% ===== Stage 5A STRUCT CHECK =====
% 
% fprintf('\n===== Stage 5A STRUCT CHECK =====\n')
% 
% % Basic info
% fprintf('nTime = %d\n', nT);
% fprintf('Valid reps: %s\n', mat2str(setAgg.reps.valid(:).'));
% fprintf('Early reps: %s\n', mat2str(setAgg.reps.early(:).'));
% fprintf('Late reps:  %s\n', mat2str(setAgg.reps.late(:).'));
% 
% binsToCheck = {'all','early','late'};
% 
% for b = 1:numel(binsToCheck)
%     bin = binsToCheck{b};
% 
%     fprintf('\n--- BIN: %s ---\n', upper(bin));
% 
%     % Check whether bin exists
%     hasCurve   = isfield(setAgg,'curve')   && isfield(setAgg.curve, bin);
%     hasBipolar = isfield(setAgg,'bipolar') && isfield(setAgg.bipolar, bin);
%     hasMap     = isfield(setAgg,'map')     && isfield(setAgg.map, bin);
% 
%     if ~hasCurve && ~hasBipolar && ~hasMap
%         fprintf('Bin "%s" is missing completely.\n', bin);
%         continue
%     end
% 
%     %% --- Curves (1x101 expected) ---
%     if hasCurve
%         curveFields = {'torque_101','angle_101','emgMeanNorm_101','centroidXc_101','centroidYc_101'};
% 
%         for f = 1:numel(curveFields)
%             fn = curveFields{f};
% 
%             if isfield(setAgg.curve.(bin), fn) && ~isempty(setAgg.curve.(bin).(fn))
%                 v = setAgg.curve.(bin).(fn);
%                 fprintf('%s size: ', fn);
%                 disp(size(v))
% 
%                 fprintf('%s numel: %d\n', fn, numel(v));
% 
%                 if numel(v) ~= nT
%                     warning('%s.%s does not have %d elements.', bin, fn, nT);
%                 end
%             else
%                 fprintf('%s is missing/empty.\n', fn);
%             end
%         end
%     else
%         fprintf('curve.%s missing.\n', bin);
%     end
% 
%     %% --- Bipolar map over time (28x101 expected) ---
%     if hasBipolar
%         fn = 'emgRMS_bipNorm_101';
% 
%         if isfield(setAgg.bipolar.(bin), fn) && ~isempty(setAgg.bipolar.(bin).(fn))
%             X = setAgg.bipolar.(bin).(fn);
%             fprintf('%s size: ', fn);
%             disp(size(X))
% 
%             if ~all(size(X) == [28 nT])
%                 warning('bipolar.%s.%s is not 28x%d.', bin, fn, nT);
%             end
% 
%             fprintf('%s range: %.2f → %.2f\n', fn, min(X(:),[],'omitnan'), max(X(:),[],'omitnan'));
%         else
%             fprintf('%s is missing/empty.\n', fn);
%         end
%     else
%         fprintf('bipolar.%s missing.\n', bin);
%     end
% 
%     %% --- Time-averaged channel vector (28x1 expected) ---
%     if hasMap
%         fn = 'emgRMS_bipNorm_map28';
% 
%         if isfield(setAgg.map.(bin), fn) && ~isempty(setAgg.map.(bin).(fn))
%             v = setAgg.map.(bin).(fn);
%             fprintf('%s size: ', fn);
%             disp(size(v))
% 
%             fprintf('%s numel: %d\n', fn, numel(v));
% 
%             if numel(v) ~= 28
%                 warning('map.%s.%s does not have 28 elements.', bin, fn);
%             end
% 
%             fprintf('%s range: %.2f → %.2f\n', fn, min(v,[],'omitnan'), max(v,[],'omitnan'));
%         else
%             fprintf('%s is missing/empty.\n', fn);
%         end
%     else
%         fprintf('map.%s missing.\n', bin);
%     end
% end





