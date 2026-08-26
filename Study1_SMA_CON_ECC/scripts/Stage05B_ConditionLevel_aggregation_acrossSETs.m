%% ===== Stage 05B: Condition-level aggregation across sets =====
%
% PURPOSE
%   Average the two sets recorded for each participant and condition into a
%   single condition-level curve, which is the unit entering the group
%   analysis.
%
% SUPPORTS
%   The condition profiles described in the Data Analysis section of the
%   manuscript, where condition profiles came from averaging both valid sets.
%
% INPUT
%   <set folder>/<tag>_stage5A_setAgg.mat        (Stage 05A)
%
% OUTPUT
%   <ROOT>/Stage5B/<P>_<COND>_stage5B_condAgg.mat
%   <ROOT>/QC_Stage5B/<P>_<COND>.png
%   <ROOT>/stage5B_summary.csv
%
% SET SELECTION
%   Every Set_* folder present is treated as retained data.
% The sets contributing to each condition are recorded in
%   condAgg.sourceFiles and in the summary table, so the grouping can be
%   checked rather than assumed.
%
% ON THE PAIRED DIFFERENCES IN THE SUMMARY
%   The CON minus ECC values printed at the end are descriptive only. They
%   reproduce what Stage 6B tests formally, and are shown here so that a
%   problem in the aggregation is visible before the group stages run.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped. A group
%   that errors is logged and the batch continues.
%
% DEPENDENCIES
%   None beyond base MATLAB. Requires R2023a or later for xregion in the QC
%   figure.

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'collapse_method_sets', "mean");  % "mean" or "median"
CFG = set_default(CFG, 'expected_nSets',       2);       % warn on any other count
CFG = set_default(CFG, 'cvp_lo',               20);      % analysis window, % of CV
CFG = set_default(CFG, 'cvp_hi',               80);
CFG = set_default(CFG, 'saveFigures',          true);

bins           = {'early','late','all'};
reqCurveFields = {'torque_101','angle_101','emgMeanNorm_101', ...
                  'centroidXc_101','centroidYc_101'};
reqBipField    = 'emgRMS_bipNorm_101';
reqMapField    = 'emgRMS_bipNorm_map';

%% ---------------------------
%  (2) Data root and file grouping
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 5B: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

outDir = fullfile(ROOT,'Stage5B');
if ~exist(outDir,'dir'), mkdir(outDir); end
qcDir = fullfile(ROOT,'QC_Stage5B');
if CFG.saveFigures && ~exist(qcDir,'dir'), mkdir(qcDir); end

F = dir(fullfile(ROOT, '**', '*_stage5A_setAgg.mat'));
assert(~isempty(F), 'Stage 5B: no Stage 5A files found under %s', ROOT);

fprintf('\n===== Stage 5B: condition-level aggregation =====\n');
fprintf('Root: %s\n', ROOT);

% Labels are read from inside each file rather than from the path, so a folder
% renamed after processing cannot regroup the data incorrectly.
keys = strings(numel(F),1);
for k = 1:numel(F)
    L = load(fullfile(F(k).folder, F(k).name),'setAgg');
    if isfield(L,'setAgg') && isfield(L.setAgg,'meta') && ...
       isfield(L.setAgg.meta,'participantName')
        keys(k) = string(L.setAgg.meta.participantName) + "|" + ...
                  string(L.setAgg.meta.conditionName);
    else
        [condPath,~]     = fileparts(F(k).folder);
        [partPath,cTmp]  = fileparts(condPath);
        [~,pTmp]         = fileparts(partPath);
        keys(k) = string(pTmp) + "|" + string(cTmp);
    end
end

groups = unique(keys,'stable');
fprintf('Found %d Stage 5A files in %d participant and condition groups.\n\n', ...
    numel(F), numel(groups));

%% ---------------------------
%  (3) Batch loop
%% ---------------------------
nG  = numel(groups);
res = table('Size',[nG 10], ...
    'VariableTypes',{'string','string','double','string','double','double', ...
                     'double','double','double','string'}, ...
    'VariableNames',{'participant','condition','nSets','sets','emgAll_pctMVA', ...
                     'centroidX_mm','centroidY_mm','torqueAll_Nm','dEMG_pct','status'});

for g = 1:nG
    parts = split(groups(g),"|");
    pName = parts(1); cName = parts(2);
    tag = pName + "_" + cName;

    try
        idx   = find(keys == groups(g));
        nSets = numel(idx);

        if nSets ~= CFG.expected_nSets
            warning('%s has %d sets, expected %d.', tag, nSets, CFG.expected_nSets);
        end

        sets        = cell(nSets,1);
        sourceFiles = strings(nSets,1);
        setNames    = strings(nSets,1);

        for s = 1:nSets
            fp = fullfile(F(idx(s)).folder, F(idx(s)).name);
            L  = load(fp,'setAgg');
            assert(isfield(L,'setAgg'), 'File missing setAgg: %s', fp);
            sets{s}        = L.setAgg;
            sourceFiles(s) = string(fp);
            if isfield(L.setAgg,'meta') && isfield(L.setAgg.meta,'setName')
                setNames(s) = string(L.setAgg.meta.setName);
            else
                [~, sn] = fileparts(F(idx(s)).folder);
                setNames(s) = string(sn);
            end
        end

        %% --- template and consistency ---
        tmpl = sets{1};
        assert(isfield(tmpl,'params') && isfield(tmpl.params,'nTime'), ...
            'Template missing params.nTime.');
        nT = tmpl.params.nTime;

        % Channel count taken from the data rather than assumed
        if isfield(tmpl.params,'nBipolar')
            nB = tmpl.params.nBipolar;
        elseif isfield(tmpl,'mapInfo') && isfield(tmpl.mapInfo,'bipolar_rowcol')
            nB = size(tmpl.mapInfo.bipolar_rowcol,1);
        else
            nB = 28;
        end

        if isfield(tmpl.params,'time_axis_values')
            xAxis = tmpl.params.time_axis_values;
        else
            xAxis = linspace(0,100,nT);
        end
        if isfield(tmpl.params,'meanAngle_101')
            meanAngle_101 = tmpl.params.meanAngle_101;
        else
            meanAngle_101 = [];
        end

        for s = 1:nSets
            assert(isfield(sets{s},'params') && isfield(sets{s}.params,'nTime') && ...
                   sets{s}.params.nTime == nT, 'nTime mismatch in set %d.', s);
        end

        %% --- structure ---
        condAgg = struct();
        condAgg.params = struct( ...
            'nSets', nSets, ...
            'collapse_method_sets', string(CFG.collapse_method_sets), ...
            'nTime', nT, 'nBipolar', nB, ...
            'time_axis_label', "%CV (0 to 100)", ...
            'time_axis_values', xAxis, ...
            'meanAngle_101', meanAngle_101, ...
            'createdOn', string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'definition', ["collapses Stage 5A outputs across sets into " ...
                           "condition-level curves"]);

        condAgg.sourceFiles = sourceFiles;
        condAgg.setNames    = setNames;

        if isfield(tmpl,'meta'),            condAgg.meta            = tmpl.meta;            else, condAgg.meta            = struct(); end
        if isfield(tmpl,'mapInfo'),         condAgg.mapInfo         = tmpl.mapInfo;         else, condAgg.mapInfo         = struct(); end
        if isfield(tmpl,'coords'),          condAgg.coords          = tmpl.coords;          else, condAgg.coords          = struct(); end
        if isfield(tmpl,'centroid_params'), condAgg.centroid_params = tmpl.centroid_params; else, condAgg.centroid_params = struct(); end

        % setName is set-specific, so it does not belong in a condition-level file
        if isfield(condAgg.meta,'setName'), condAgg.meta = rmfield(condAgg.meta,'setName'); end

        %% --- collapse across sets, one bin at a time ---
        for b = 1:numel(bins)
            bin = bins{b};

            for f = 1:numel(reqCurveFields)
                fn = reqCurveFields{f};
                C  = gather_curves(sets, bin, fn, nT);
                if isempty(C)
                    condAgg.curve.(bin).(fn) = [];
                else
                    condAgg.curve.(bin).(fn) = collapse_dim(cat(1,C{:}), ...
                        CFG.collapse_method_sets, 1);
                end
                condAgg.nSetsUsed.(bin).(fn) = numel(C);
            end

            B = gather_bipolar(sets, bin, reqBipField, nB, nT);
            if isempty(B)
                condAgg.bipolar.(bin).(reqBipField) = [];
            else
                condAgg.bipolar.(bin).(reqBipField) = collapse_dim(cat(3,B{:}), ...
                    CFG.collapse_method_sets, 3);
            end
            condAgg.nSetsUsed.(bin).bipolar = numel(B);

            V = gather_maps(sets, bin, reqMapField, nB);
            if isempty(V)
                condAgg.map.(bin).(reqMapField) = [];
            else
                condAgg.map.(bin).(reqMapField) = collapse_dim(cat(2,V{:}), ...
                    CFG.collapse_method_sets, 2);
            end
            condAgg.nSetsUsed.(bin).map = numel(V);
        end

        %% --- summary values ---
        inCVP = xAxis >= CFG.cvp_lo & xAxis <= CFG.cvp_hi;

        emgAll = mean(condAgg.curve.all.emgMeanNorm_101(inCVP),'omitnan');
        cxAll  = mean(condAgg.curve.all.centroidXc_101(inCVP),'omitnan');
        cyAll  = mean(condAgg.curve.all.centroidYc_101(inCVP),'omitnan');
        tqAll  = mean(condAgg.curve.all.torque_101(inCVP),'omitnan');

        if ~isempty(condAgg.curve.late.emgMeanNorm_101) && ...
           ~isempty(condAgg.curve.early.emgMeanNorm_101)
            dEMG = mean(condAgg.curve.late.emgMeanNorm_101(inCVP),'omitnan') - ...
                   mean(condAgg.curve.early.emgMeanNorm_101(inCVP),'omitnan');
        else
            dEMG = NaN;
        end

        %% --- QC figure ---
        if CFG.saveFigures
            plot_stage5b(condAgg, xAxis, CFG.cvp_lo, CFG.cvp_hi, tag, nSets, ...
                reqMapField, qcDir);
        end

        %% --- save ---
        outName = fullfile(outDir, sprintf('%s_stage5B_condAgg.mat', ...
            regexprep(tag,'[^\w-]','_')));
        save(outName, 'condAgg', '-v7.3');

        res(g,:) = {pName, cName, nSets, strjoin(setNames,";"), emgAll, cxAll, ...
                    cyAll, tqAll, dEMG, "ok"};

        fprintf(['[%3d/%3d] %-14s %d sets (%-15s) | EMG %5.1f %%MVA | ' ...
                 'X %+5.2f Y %+5.2f mm | torque %5.1f Nm\n'], ...
            g, nG, tag, nSets, strjoin(setNames,";"), emgAll, cxAll, cyAll, tqAll);

    catch ME
        res(g,:) = {pName, cName, NaN, "", NaN, NaN, NaN, NaN, NaN, string(ME.message)};
        fprintf('[%3d/%3d] %-14s FAILED: %s\n', g, nG, tag, ME.message);
    end
end

%% ---------------------------
%  (4) Summary
%% ---------------------------
ok = res.status == "ok";
fprintf('\n===== STAGE 5B SUMMARY =====\n');
fprintf('Groups %d | ok %d | failed %d\n', nG, nnz(ok), nnz(~ok));

if any(ok)
    odd = res(ok & res.nSets ~= CFG.expected_nSets, :);
    if isempty(odd)
        fprintf('Every group contains %d sets.\n', CFG.expected_nSets);
    else
        fprintf('\nGroups with a different number of sets:\n');
        disp(odd(:, {'participant','condition','nSets','sets'}));
    end

    fprintf('\nCondition-level means over the %d-%d%% CVP:\n', CFG.cvp_lo, CFG.cvp_hi);
    fprintf('  %-8s %14s %14s %14s %12s\n', ...
        'cond','EMG (%MVA)','X (mm)','Y (mm)','torque (Nm)');
    for cc = ["CON_75","ECC_75","CON_90","ECC_90"]
        sel = ok & res.condition == cc;
        if any(sel)
            fprintf('  %-8s %7.1f +/-%5.1f %7.2f +/-%5.2f %7.2f +/-%5.2f %7.1f +/-%4.1f\n', cc, ...
                mean(res.emgAll_pctMVA(sel)),  std(res.emgAll_pctMVA(sel)), ...
                mean(res.centroidX_mm(sel)),   std(res.centroidX_mm(sel)), ...
                mean(res.centroidY_mm(sel)),   std(res.centroidY_mm(sel)), ...
                mean(res.torqueAll_Nm(sel)),   std(res.torqueAll_Nm(sel)));
        end
    end

    % Descriptive only, see the header note.
    P = unique(res.participant(ok));
    fprintf('\nPaired CON minus ECC differences (n = %d participants):\n', numel(P));
    for I = ["75","90"]
        dE = nan(numel(P),1); dX = dE; dY = dE;
        for i = 1:numel(P)
            ci = ok & res.participant == P(i) & res.condition == "CON_"+I;
            ei = ok & res.participant == P(i) & res.condition == "ECC_"+I;
            if any(ci) && any(ei)
                dE(i) = res.emgAll_pctMVA(ci) - res.emgAll_pctMVA(ei);
                dX(i) = res.centroidX_mm(ci)  - res.centroidX_mm(ei);
                dY(i) = res.centroidY_mm(ci)  - res.centroidY_mm(ei);
            end
        end
        fprintf('  %s%%: EMG %+6.2f +/- %5.2f | X %+6.3f +/- %5.3f | Y %+6.3f +/- %5.3f\n', ...
            I, mean(dE,'omitnan'), std(dE,'omitnan'), ...
               mean(dX,'omitnan'), std(dX,'omitnan'), ...
               mean(dY,'omitnan'), std(dY,'omitnan'));
    end
end

if any(~ok)
    fprintf('\nFailures:\n');
    disp(res(~ok, {'participant','condition','status'}));
end

writetable(res, fullfile(ROOT,'stage5B_summary.csv'));
fprintf('\nSaved summary  : %s\n', fullfile(ROOT,'stage5B_summary.csv'));
fprintf('condAgg files  : %s\n', outDir);
if CFG.saveFigures, fprintf('QC figures     : %s\n', qcDir); end

%% =========================================================
%  Local functions
%% =========================================================
function S = set_default(S, fieldName, value)
% Fills a config field only when the caller has not already supplied it.
    if ~isfield(S, fieldName) || isempty(S.(fieldName))
        S.(fieldName) = value;
    end
end

function validate_root(ROOT)
% Confirms the chosen folder looks like the expected layout, so that a wrong
% selection fails immediately with a useful message rather than partway through
% a batch.
    assert(isfolder(ROOT), 'Data root not found at %s', ROOT);

    D = dir(ROOT);
    D = D([D.isdir] & ~startsWith({D.name},'.'));

    hasParticipant = false;
    for i = 1:numel(D)
        pPath = fullfile(ROOT, D(i).name);
        C = dir(pPath);
        C = C([C.isdir]);
        if any(~cellfun(@isempty, regexp({C.name},'^(CON|ECC)_\d+$','once')))
            hasParticipant = true;
            break;
        end
    end

       % A repository containing only the processed group file is a valid root for
    % the Stage 6 analyses, which read nothing else.
    hasGroupFile = isfile(fullfile(ROOT,'groupData_stage6A.mat'));

    assert(hasParticipant || hasGroupFile, ...
        ['The selected folder does not look like the expected data root.\n' ...
         'Expected either the raw layout\n' ...
         '  <ROOT>/<participant>/<CON_75|CON_90|ECC_75|ECC_90>/<Set_N>/\n' ...
         'or a folder containing groupData_stage6A.mat\n' ...
         'Selected: %s'], ROOT);
end

function Y = collapse_dim(X, method, dim)
    if isempty(X), Y = []; return; end
    switch lower(string(method))
        case "mean",   Y = mean(X, dim, 'omitnan');
        case "median", Y = median(X, dim, 'omitnan');
        otherwise, error('Unknown collapse method: %s', method);
    end
end

function C = gather_curves(sets, bin, fn, nT)
% Collects one curve field from every set that carries it at the right length.
% A set missing the field is skipped rather than failing the group, and the
% count of contributing sets is recorded by the caller.
    C = {};
    for s = 1:numel(sets)
        A = sets{s};
        if ~isfield(A,'curve') || ~isfield(A.curve,bin) || ~isfield(A.curve.(bin),fn), continue; end
        v = A.curve.(bin).(fn);
        if isempty(v), continue; end
        v = reshape(v,1,[]);
        if numel(v) ~= nT, continue; end
        C{end+1} = v; %#ok<AGROW>
    end
end

function C = gather_bipolar(sets, bin, fn, nB, nT)
    C = {};
    for s = 1:numel(sets)
        A = sets{s};
        if ~isfield(A,'bipolar') || ~isfield(A.bipolar,bin) || ~isfield(A.bipolar.(bin),fn), continue; end
        X = A.bipolar.(bin).(fn);
        if isempty(X) || ~ismatrix(X) || size(X,1) ~= nB || size(X,2) ~= nT, continue; end
        C{end+1} = X; %#ok<AGROW>
    end
end

function C = gather_maps(sets, bin, fn, nB)
    C = {};
    for s = 1:numel(sets)
        A = sets{s};
        if ~isfield(A,'map') || ~isfield(A.map,bin) || ~isfield(A.map.(bin),fn), continue; end
        v = A.map.(bin).(fn);
        if isempty(v), continue; end
        v = v(:);
        if numel(v) ~= nB, continue; end
        C{end+1} = v; %#ok<AGROW>
    end
end

function plot_stage5b(condAgg, xAxis, lo, hi, tag, nSets, mapField, qcDir)
% The five condition-level curves, with the time-averaged spatial map in the
% sixth tile so that the amplitude distribution can be seen alongside them.
    hasLate = isfield(condAgg.curve,'late') && ...
              isfield(condAgg.curve.late,'torque_101') && ...
              ~isempty(condAgg.curve.late.torque_101);

    f = figure('Color','w','Visible','off','Position',[100 100 1100 800]);
    tiledlayout(f,3,2,'TileSpacing','compact','Padding','compact');

    panels = {'torque_101','Torque (Nm)'; 'angle_101','Angle (deg)'; ...
              'emgMeanNorm_101','EMG (%MVA)'; ...
              'centroidXc_101','Centroid X (mm)'; ...
              'centroidYc_101','Centroid Y (mm)'};

    for i = 1:size(panels,1)
        ax = nexttile; hold(ax,'on');
        plot(ax, xAxis, condAgg.curve.all.(panels{i,1}), 'k','LineWidth',2);
        plot(ax, xAxis, condAgg.curve.early.(panels{i,1}), 'Color',[0 .45 .74],'LineWidth',1.5);
        if hasLate
            plot(ax, xAxis, condAgg.curve.late.(panels{i,1}), 'Color',[.85 .33 .1],'LineWidth',1.5);
        end
        if contains(panels{i,1},'centroid'), yline(ax, 0, '--'); end
        xregion(ax, lo, hi, 'FaceAlpha',0.07);
        ylabel(ax, panels{i,2},'FontWeight','bold');
        if i >= 4, xlabel(ax,'%CV','FontWeight','bold'); end
        set(ax,'Box','off','TickDir','out');
        if i == 1
            if hasLate, legend(ax,{'All','Early','Late'},'Box','off','Location','best');
            else,       legend(ax,{'All','Early'},       'Box','off','Location','best');
            end
        end
    end

    ax6 = nexttile;
    if isfield(condAgg,'mapInfo') && isfield(condAgg.mapInfo,'bipolar_rowcol') && ...
       isfield(condAgg.map,'all') && isfield(condAgg.map.all, mapField) && ...
       ~isempty(condAgg.map.all.(mapField))
        rc = condAgg.mapInfo.bipolar_rowcol;
        m  = condAgg.map.all.(mapField);
        M  = nan(max(rc(:,1)), max(rc(:,2)));
        for ch = 1:numel(m), M(rc(ch,1), rc(ch,2)) = m(ch); end
        imagesc(ax6, M); colorbar(ax6); axis(ax6,'image');
        xlabel(ax6,'Column (1 = lateral)'); ylabel(ax6,'Row (1 = proximal)');
        title(ax6,'Mean spatial map (%MVA)','FontWeight','bold');
    else
        axis(ax6,'off');
    end

    sgtitle(f, tag + sprintf("  |  %d sets averaged", nSets), ...
        'FontWeight','bold','Interpreter','none');

    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + ".png"), 'Resolution',110);
    close(f);
end
