%% ===== Stage 05A: Set-level aggregation across repetitions =====
%
% PURPOSE
%   Collapse the repetition-level curves from Stages 3 and 4 into set-level
%   curves, for the whole set and for its first and last repetitions (i.e., early, late).
%
% SUPPORTS
%   The condition profiles and the early to late contrast described in the Data
%   Analysis section of the manuscript.
%
% INPUT
%   <set folder>/<tag>_stage3_EMGnorm.mat        (Stage 03)
%   <set folder>/<tag>_stage4_centroid_mm.mat    (Stage 04)
%
% OUTPUT
%   <set folder>/<tag>_stage5_setAgg.mat
%   <ROOT>/QC_Stage5/<tag>.png
%   <ROOT>/stage5_summary.csv
%
% BINS
%   all    mean across every valid repetition
%   early  mean across the first CFG.nEarly valid repetitions
%   late   mean across the last  CFG.nLate  valid repetitions
%
%   No new repetition categories are created here. This stage only averages
%   what Stages 3 and 4 produced, and a repetition must be valid in both to be
%   included.
%
% ON THE EARLY TO LATE CONTRAST
%   Sets at 75% MVC hold 8 repetitions, so early (1 and 2) and late (7 and 8)
%   are separated by four intervening contractions. Sets at 90% hold 4, so
%   early (1 and 2) and late (3 and 4) are adjacent.
%
% ON THE TWO MAP FIELDS
%   setAgg.mapInfo carries the electrode map from Stage 03, describing where
%   each channel sits. setAgg.map carries the time-averaged amplitude map for
%   each bin. They are different things despite the similar names.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped. A file
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

CFG = set_default(CFG, 'nEarly',          2);
CFG = set_default(CFG, 'nLate',           2);
CFG = set_default(CFG, 'collapse_method', "mean");    % "mean" or "median"
CFG = set_default(CFG, 'modeIfTooFew',    "error");   % "error" | "shrink" | "skipLate"
CFG = set_default(CFG, 'cvp_lo',          20);        % analysis window, % of the CV phase
CFG = set_default(CFG, 'cvp_hi',          80);
CFG = set_default(CFG, 'saveFigures',     true);

%% ---------------------------
%  (2) Data root and file list
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 5: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

qcDir = fullfile(ROOT,'QC_Stage5');
if CFG.saveFigures && ~exist(qcDir,'dir'), mkdir(qcDir); end

F3 = dir(fullfile(ROOT, '**', '*_stage3_EMGnorm.mat'));
assert(~isempty(F3), 'Stage 5: no Stage 3 files found under %s', ROOT);

fprintf('\n===== Stage 5: set-level aggregation =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Found %d Stage 3 files.\n\n', numel(F3));

%% ---------------------------
%  (3) Batch loop
%% ---------------------------
nF  = numel(F3);
res = table('Size',[nF 11], ...
    'VariableTypes',{'string','string','string','double','string','string', ...
                     'double','double','double','double','string'}, ...
    'VariableNames',{'participant','condition','set','nValidReps','earlyReps', ...
                     'lateReps','emgEarly_pctMVA','emgLate_pctMVA','dEMG_pct', ...
                     'dCentroidY_mm','status'});

for k = 1:nF
    fp3 = fullfile(F3(k).folder, F3(k).name);
    tag = "";
    pName = ""; cName = ""; sName = "";

    try
        %% --- load the pair ---
        S3 = load(fp3,'stage3');
        assert(isfield(S3,'stage3'), 'Stage 3 file has no "stage3" variable.');
        stage3 = S3.stage3;

        [~, base] = fileparts(F3(k).name);
        base = strrep(base, '_stage3_EMGnorm', '');
        fp4  = fullfile(F3(k).folder, sprintf('%s_stage4_centroid_mm.mat', base));
        assert(isfile(fp4), 'No matching Stage 4 file for %s', base);

        S4 = load(fp4,'stage4');
        assert(isfield(S4,'stage4'), 'Stage 4 file has no "stage4" variable.');
        stage4 = S4.stage4;

        if isfield(stage3,'meta')
            pName = string(stage3.meta.participantName);
            cName = string(stage3.meta.conditionName);
            sName = string(stage3.meta.setName);
        else
            pName = "?"; cName = "?"; sName = "?";
        end
        tag = pName + "_" + cName + "_" + sName;

        nReps3 = numel(stage3.rep);
        nReps4 = numel(stage4.rep);
        assert(nReps3 == nReps4, 'Rep counts differ: Stage 3 %d, Stage 4 %d.', nReps3, nReps4);

        nT = stage3.params.timeNorm_points;

        % Channel count taken from the map rather than assumed, so the checks
        % below stay valid if the grid ever changes.
        if isfield(stage3,'map') && isfield(stage3.map,'bipolar_rowcol')
            nB = size(stage3.map.bipolar_rowcol,1);
        else
            nB = 28;
        end

        if isfield(stage3,'groupReady') && isfield(stage3.groupReady,'x_norm_101')
            xAxis = stage3.groupReady.x_norm_101;
        else
            xAxis = linspace(0,100,nT);
        end

        if isfield(stage3,'groupReady') && isfield(stage3.groupReady,'meanAngle_101')
            meanAngle_101 = stage3.groupReady.meanAngle_101;
        else
            meanAngle_101 = [];
        end

        %% --- valid repetitions, early and late selection ---
        valid3 = arrayfun(@(s) logical(s.is_valid), stage3.rep);
        valid4 = arrayfun(@(s) logical(s.is_valid), stage4.rep);
        repValid = find(valid3 & valid4);
        assert(~isempty(repValid), 'No repetitions valid in both stages.');

        [repEarly, repLate] = pick_early_late(repValid, CFG.nEarly, CFG.nLate, ...
                                              CFG.modeIfTooFew);

        %% --- structure ---
        setAgg = struct();
        if isfield(stage3,'meta'), setAgg.meta = stage3.meta; else, setAgg.meta = struct(); end

        setAgg.params = struct( ...
            'nTime', nT, 'nBipolar', nB, ...
            'nEarly', CFG.nEarly, 'nLate', CFG.nLate, ...
            'collapse_method', string(CFG.collapse_method), ...
            'modeIfTooFew', string(CFG.modeIfTooFew), ...
            'time_axis_label', "%CV (0 to 100)", ...
            'time_axis_values', xAxis, ...
            'meanAngle_101', meanAngle_101, ...
            'createdOn', string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'definition', ["collapses repetition-level outputs into set-level " ...
                           "all, early and late curves"]);

        setAgg.reps = struct('valid', repValid(:), 'early', repEarly(:), 'late', repLate(:));

        if isfield(stage3,'map'),    setAgg.mapInfo = stage3.map;    else, setAgg.mapInfo = struct(); end
        if isfield(stage4,'coords'), setAgg.coords  = stage4.coords; else, setAgg.coords  = struct(); end
        if isfield(stage4,'params'), setAgg.centroid_params = stage4.params; else, setAgg.centroid_params = struct(); end

        %% --- collapse each bin ---
        bins = {'all', repValid; 'early', repEarly; 'late', repLate};

        for b = 1:size(bins,1)
            binName = bins{b,1};
            reps    = bins{b,2};

            if isempty(reps)
                setAgg.curve.(binName)   = struct();
                setAgg.bipolar.(binName) = struct();
                setAgg.map.(binName)     = struct();
                continue
            end

            setAgg.curve.(binName).torque_101      = collapse_rows(stage3.rep, reps, 'torque_101',                 CFG.collapse_method);
            setAgg.curve.(binName).angle_101       = collapse_rows(stage3.rep, reps, 'angle_101',                  CFG.collapse_method);
            setAgg.curve.(binName).emgMeanNorm_101 = collapse_rows(stage3.rep, reps, 'emgRMSmeanNorm_101',         CFG.collapse_method);
            setAgg.curve.(binName).centroidXc_101  = collapse_rows(stage4.rep, reps, 'centroidX_centered_mm_101',  CFG.collapse_method);
            setAgg.curve.(binName).centroidYc_101  = collapse_rows(stage4.rep, reps, 'centroidY_centered_mm_101',  CFG.collapse_method);

            B = cat(3, stage3.rep(reps).emgRMS_bipNorm_101);     % [nB x nT x nRepsInBin]
            setAgg.bipolar.(binName).emgRMS_bipNorm_101 = collapse_dim(B, CFG.collapse_method, 3);

            % Time-averaged amplitude map, one value per channel
            setAgg.map.(binName).emgRMS_bipNorm_map = ...
                mean(setAgg.bipolar.(binName).emgRMS_bipNorm_101, 2, 'omitnan');
        end

        %% --- size checks ---
        for binName = ["all","early"]
            assert(numel(setAgg.curve.(binName).torque_101) == nT, ...
                '%s torque curve is not %d points.', binName, nT);
            assert(isequal(size(setAgg.bipolar.(binName).emgRMS_bipNorm_101), [nB nT]), ...
                '%s bipolar matrix is not %dx%d.', binName, nB, nT);
        end
        if ~isempty(repLate)
            assert(numel(setAgg.curve.late.torque_101) == nT, ...
                'Late curve is not %d points.', nT);
            assert(isequal(size(setAgg.bipolar.late.emgRMS_bipNorm_101), [nB nT]), ...
                'Late bipolar matrix is not %dx%d.', nB, nT);
        end

        %% --- summary values over the analysis window ---
        inCVP = xAxis >= CFG.cvp_lo & xAxis <= CFG.cvp_hi;

        emgE = mean(setAgg.curve.early.emgMeanNorm_101(inCVP),'omitnan');
        cyE  = mean(setAgg.curve.early.centroidYc_101(inCVP),'omitnan');

        if ~isempty(repLate)
            emgL = mean(setAgg.curve.late.emgMeanNorm_101(inCVP),'omitnan');
            cyL  = mean(setAgg.curve.late.centroidYc_101(inCVP),'omitnan');
        else
            emgL = NaN; cyL = NaN;
        end

        %% --- QC figure ---
        if CFG.saveFigures
            plot_stage5(setAgg, xAxis, CFG.cvp_lo, CFG.cvp_hi, tag, qcDir);
        end

        %% --- save ---
                outName = fullfile(F3(k).folder, sprintf('%s_stage5A_setAgg.mat', base));
        save(outName, 'setAgg', '-v7.3');

        res(k,:) = {pName, cName, sName, numel(repValid), ...
                    string(mat2str(repEarly)), string(mat2str(repLate)), ...
                    emgE, emgL, emgL-emgE, cyL-cyE, "ok"};

        fprintf(['[%3d/%3d] %-26s %d reps | early %-7s late %-7s | ' ...
                 'EMG %5.1f -> %5.1f (%+5.1f) | dY %+5.2f mm\n'], ...
            k, nF, tag, numel(repValid), mat2str(repEarly), mat2str(repLate), ...
            emgE, emgL, emgL-emgE, cyL-cyE);

    catch ME
        res(k,:) = {pName, cName, sName, NaN, "", "", ...
                    NaN, NaN, NaN, NaN, string(ME.message)};
        fprintf('[%3d/%3d] %-26s FAILED: %s\n', k, nF, tag, ME.message);
    end
end

%% ---------------------------
%  (4) Summary
%% ---------------------------
ok = res.status == "ok";
fprintf('\n===== STAGE 5 SUMMARY =====\n');
fprintf('Processed %d | ok %d | failed %d\n', nF, nnz(ok), nnz(~ok));

if any(ok)
    fprintf('\nEarly to late change over the %d-%d%% CVP, by condition:\n', ...
        CFG.cvp_lo, CFG.cvp_hi);
    fprintf('  %-8s %16s %16s\n', 'cond', 'dEMG (%MVA)', 'dCentroidY (mm)');
    for cc = ["CON_75","ECC_75","CON_90","ECC_90"]
        sel = ok & res.condition == cc;
        if any(sel)
            fprintf('  %-8s %8.2f +/- %5.2f %8.2f +/- %5.2f  (n = %d)\n', cc, ...
                mean(res.dEMG_pct(sel),'omitnan'), std(res.dEMG_pct(sel),'omitnan'), ...
                mean(res.dCentroidY_mm(sel),'omitnan'), std(res.dCentroidY_mm(sel),'omitnan'), ...
                nnz(sel));
        end
    end

    fprintf('\nRepetitions per set:\n');
    for n = unique(res.nValidReps(ok)).'
        fprintf('  %d reps : %d file(s)\n', n, nnz(ok & res.nValidReps == n));
    end
    fprintf(['\nAt 8 repetitions the early and late bins are separated by four\n' ...
             'contractions; at 4 they are adjacent. The two intensities are\n' ...
             'therefore not directly comparable on the early to late contrast.\n']);
end

if any(~ok)
    fprintf('\nFailures:\n');
    disp(res(~ok, {'participant','condition','set','status'}));
end

writetable(res, fullfile(ROOT,'stage5_summary.csv'));
fprintf('\nSaved summary: %s\n', fullfile(ROOT,'stage5_summary.csv'));
if CFG.saveFigures, fprintf('QC figures   : %s\n', qcDir); end

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

function y = collapse_rows(repStruct, reps, field, method)
% Stacks one field across the selected repetitions and collapses to one curve.
    C = arrayfun(@(r) reshape(repStruct(r).(field), 1, []), reps, 'UniformOutput', false);
    y = collapse_dim(cat(1, C{:}), method, 1);
end

function Y = collapse_dim(X, method, dim)
    if isempty(X), Y = []; return; end
    switch lower(string(method))
        case "mean",   Y = mean(X, dim, 'omitnan');
        case "median", Y = median(X, dim, 'omitnan');
        otherwise, error('Unknown collapse method: %s', method);
    end
end

function [earlyReps, lateReps] = pick_early_late(validReps, nEarly, nLate, modeIfTooFew)
% Takes the first nEarly and the last nLate valid repetitions. The two bins
% never overlap when the set holds enough repetitions; what happens when it
% does not is governed by modeIfTooFew.
    validReps = validReps(:).';
    nV = numel(validReps);

    if nV >= (nEarly + nLate)
        earlyReps = validReps(1:nEarly);
        lateReps  = validReps(end-nLate+1:end);
        return
    end

    switch lower(string(modeIfTooFew))
        case "error"
            error('Not enough valid repetitions (%d) for nEarly = %d and nLate = %d.', ...
                nV, nEarly, nLate);
        case "shrink"
            nE = min(nEarly, nV);
            nL = min(nLate, nV - nE);
            earlyReps = validReps(1:nE);
            lateReps  = validReps(max(1, nV-nL+1):nV);
        case "skiplate"
            earlyReps = validReps(1:min(nEarly, nV));
            lateReps  = [];
        otherwise
            error('Unknown modeIfTooFew: %s', modeIfTooFew);
    end
end

function plot_stage5(setAgg, xAxis, lo, hi, tag, qcDir)
% All five set-level curves, with the whole set against its early and late
% halves, and the analysis window shaded.
    hasLate = isfield(setAgg.curve,'late') && isfield(setAgg.curve.late,'torque_101') ...
              && ~isempty(setAgg.curve.late.torque_101);

    f = figure('Color','w','Visible','off','Position',[100 100 1100 800]);
    tiledlayout(f,3,2,'TileSpacing','compact','Padding','compact');

    panels = {'torque_101','Torque (Nm)'; 'angle_101','Angle (deg)'; ...
              'emgMeanNorm_101','EMG (%MVA)'; ...
              'centroidXc_101','Centroid X (mm)'; ...
              'centroidYc_101','Centroid Y (mm)'};

    for i = 1:size(panels,1)
        ax = nexttile; hold(ax,'on');
        plot(ax, xAxis, setAgg.curve.all.(panels{i,1}),   'k',  'LineWidth',2);
        plot(ax, xAxis, setAgg.curve.early.(panels{i,1}), 'Color',[0 .45 .74], 'LineWidth',1.5);
        if hasLate
            plot(ax, xAxis, setAgg.curve.late.(panels{i,1}), 'Color',[.85 .33 .1], 'LineWidth',1.5);
        end
        if contains(panels{i,1},'centroid'), yline(ax, 0, '--'); end
        xregion(ax, lo, hi, 'FaceAlpha',0.07);
        ylabel(ax, panels{i,2}, 'FontWeight','bold');
        if i >= 4, xlabel(ax,'%CV','FontWeight','bold'); end
        set(ax,'Box','off','TickDir','out');
        if i == 1
            if hasLate, legend(ax, {'All','Early','Late'}, 'Box','off','Location','best');
            else,       legend(ax, {'All','Early'},        'Box','off','Location','best');
            end
        end
    end

    nexttile; axis off;    % sixth tile of the 3x2 layout, deliberately empty
    sgtitle(f, tag + " | set-level aggregation", 'FontWeight','bold','Interpreter','none');

    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + ".png"), 'Resolution',110);
    close(f);
end
