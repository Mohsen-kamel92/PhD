%% ===== Stage 04: Amplitude-weighted centroid in physical coordinates =====
%
% PURPOSE
%   Reduce the bipolar amplitude map at each instant to a single point, the
%   amplitude-weighted centroid, expressed in millimetres on the electrode
%   grid.
%
% SUPPORTS
%   The centroid definition in the Analysis of Spatial Activity section of the
%   manuscript.
%
% INPUT
%   <set folder>/<tag>_stage3_EMGnorm.mat                 (Stage 03)
%
% OUTPUT
%   <set folder>/<tag>_stage4_centroid_mm.mat
%   <ROOT>/QC_Stage4/<tag>.png
%   <ROOT>/stage4_summary.csv
%
% COORDINATE SYSTEM
%   The grid is 8 monopolar rows by 4 columns at 10 mm spacing. Bipolar
%   derivation along the proximal to distal axis leaves 7 bipolar rows by 4
%   columns, each bipolar channel sitting at the midpoint between its two
%   monopolar electrodes:
%       x(col) = (col-1)*IED           gives  0, 10, 20, 30 mm, lateral to medial
%       y(row) = (row-1)*IED + IED/2   gives  5, 15, ..., 65 mm, proximal to distal
%   The geometric centre of the array is therefore (15, 35) mm, and the centred
%   coordinates express displacement relative to that point.
%
% WHY THE WEIGHTS ARE VALID HERE
%   The centroid is amplitude-weighted, so it depends on the relative values
%   between channels. Stage 3 normalizes with a single scalar per participant,
%   which cancels in the ratio and leaves the spatial pattern intact. The
%   centroid computed from %MVA is therefore identical to one computed from the
%   raw bipolar RMS.
%
% ON NON-FINITE WEIGHTS
%   A channel carrying no value at a given instant is treated as contributing
%   nothing, so it drops out of both the numerator and the denominator. The
%   number of such samples is counted and reported, and is expected to be zero
%   for a clean recording.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped. A file
%   that errors is logged and the batch continues.
%
% DEPENDENCIES
%   None beyond base MATLAB.

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'IED_mm',      10);     % inter-electrode distance
CFG = set_default(CFG, 'eps_denom',   1e-12);  % guard for the weight sum
CFG = set_default(CFG, 'cvp_lo',      20);     % analysis window, % of the CV phase
CFG = set_default(CFG, 'cvp_hi',      80);
CFG = set_default(CFG, 'saveFigures', true);

%% ---------------------------
%  (2) Data root and file list
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 4: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

qcDir = fullfile(ROOT,'QC_Stage4');
if CFG.saveFigures && ~exist(qcDir,'dir'), mkdir(qcDir); end

F = dir(fullfile(ROOT, '**', '*_stage3_EMGnorm.mat'));
assert(~isempty(F), 'Stage 4: no Stage 3 files found under %s', ROOT);

fprintf('\n===== Stage 4: weighted centroid =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Found %d Stage 3 files.\n\n', numel(F));

%% ---------------------------
%  (3) Batch loop
%% ---------------------------
nF  = numel(F);
res = table('Size',[nF 12], ...
    'VariableTypes',{'string','string','string','double','double','double', ...
                     'double','double','double','double','double','string'}, ...
    'VariableNames',{'participant','condition','set','nReps','nValid', ...
                     'meanX_mm','meanY_mm','rangeX_mm','rangeY_mm', ...
                     'sdAcrossReps_mm','nNonFiniteWeights','status'});

for k = 1:nF
    fp  = fullfile(F(k).folder, F(k).name);
    tag = "";
    pName = ""; cName = ""; sName = "";

    try
        L = load(fp);
        assert(isfield(L,'stage3'), 'File missing stage3 struct.');
        stage3 = L.stage3;

        assert(isfield(stage3,'rep') && isfield(stage3,'map') && isfield(stage3,'params'), ...
            'stage3 missing rep, map or params.');
        assert(isfield(stage3.map,'bipolar_rowcol'), 'stage3.map missing bipolar_rowcol.');

        if isfield(stage3,'meta')
            pName = string(stage3.meta.participantName);
            cName = string(stage3.meta.conditionName);
            sName = string(stage3.meta.setName);
        else
            [condPath, sTmp] = fileparts(F(k).folder);
            [~,        cTmp] = fileparts(condPath);
            pName = "?"; cName = string(cTmp); sName = string(sTmp);
        end
        tag = pName + "_" + cName + "_" + sName;

        rc    = stage3.map.bipolar_rowcol;
        nB    = size(rc,1);
        nReps = numel(stage3.rep);
        nT    = stage3.params.timeNorm_points;

        if isfield(stage3,'groupReady') && isfield(stage3.groupReady,'x_norm_101')
            x_norm_101 = stage3.groupReady.x_norm_101;
        else
            x_norm_101 = linspace(0,100,nT);
        end

        if isfield(stage3,'groupReady') && isfield(stage3.groupReady,'meanAngle_101')
            meanAngle_101 = stage3.groupReady.meanAngle_101;
        else
            meanAngle_101 = [];
        end

        %% --- physical coordinates ---
        x_i_mm = (rc(:,2) - 1) * CFG.IED_mm;
        y_i_mm = ((rc(:,1) - 1) * CFG.IED_mm) + (CFG.IED_mm/2);

        x_center_mm = mean([min(x_i_mm) max(x_i_mm)]);
        y_center_mm = mean([min(y_i_mm) max(y_i_mm)]);

        %% --- output struct ---
        stage4 = struct();
        stage4.meta = stage3.meta;
        stage4.map  = stage3.map;

        stage4.params = struct( ...
            'IED_mm', CFG.IED_mm, ...
            'coord_x', "lateral_to_medial_mm", ...
            'coord_y', "proximal_to_distal_mm", ...
            'origin_note', "x = 0 at column 1 (lateral), y = 5 mm at bipolar row 1", ...
            'x_center_mm', x_center_mm, ...
            'y_center_mm', y_center_mm, ...
            'weights', "emgRMS_bipNorm_101", ...
            'weights_note', ["global scalar normalization at Stage 3, so the " ...
                             "centroid is identical to one computed from raw RMS"], ...
            'nonfinite_note', "non-finite weights treated as zero, see the header", ...
            'x_norm_101', x_norm_101, ...
            'meanAngle_101', meanAngle_101, ...
            'eps_denom', CFG.eps_denom, ...
            'createdOn', string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));

        stage4.coords = struct('x_i_mm', x_i_mm, 'y_i_mm', y_i_mm, ...
            'x_center_mm', x_center_mm, 'y_center_mm', y_center_mm);

        stage4.rep = repmat(struct( ...
            'is_valid', false, ...
            'centroidX_mm_101', nan(1,nT), ...
            'centroidY_mm_101', nan(1,nT), ...
            'centroidX_centered_mm_101', nan(1,nT), ...
            'centroidY_centered_mm_101', nan(1,nT), ...
            'sumW_101', nan(1,nT)), nReps, 1);

        %% --- centroid per repetition ---
        nNonFinite = 0;

        for r = 1:nReps
            stage4.rep(r).is_valid = logical(stage3.rep(r).is_valid);
            if ~stage4.rep(r).is_valid, continue; end

            W = stage3.rep(r).emgRMS_bipNorm_101;
            if isempty(W) || ~ismatrix(W) || size(W,1) ~= nB || size(W,2) ~= nT
                warning('%s rep %d: weight matrix is %s, expected %dx%d. Skipped.', ...
                    tag, r, mat2str(size(W)), nB, nT);
                stage4.rep(r).is_valid = false;
                continue
            end

            % Zero rather than NaN, so a missing channel drops out of both the
            % numerator and the denominator instead of voiding the whole time
            % point. See the header note.
            bad = ~isfinite(W);
            nNonFinite = nNonFinite + nnz(bad);
            W(bad)   = 0;
            W(W < 0) = 0;                      % an RMS cannot be negative

            sumW      = sum(W, 1);
            sumW_safe = max(sumW, CFG.eps_denom);

            cx = (x_i_mm.' * W) ./ sumW_safe;
            cy = (y_i_mm.' * W) ./ sumW_safe;

            stage4.rep(r).centroidX_mm_101 = cx;
            stage4.rep(r).centroidY_mm_101 = cy;
            stage4.rep(r).centroidX_centered_mm_101 = cx - x_center_mm;
            stage4.rep(r).centroidY_centered_mm_101 = cy - y_center_mm;
            stage4.rep(r).sumW_101 = sumW;
        end

        stage4.params.nNonFiniteWeights = nNonFinite;

        %% --- summary values over the analysis window ---
        inCVP  = x_norm_101 >= CFG.cvp_lo & x_norm_101 <= CFG.cvp_hi;
        valid  = find(arrayfun(@(s) s.is_valid, stage4.rep));
        nValid = numel(valid);

        if nValid == 0
            mX = NaN; mY = NaN; rX = NaN; rY = NaN; sdReps = NaN;
        else
            CX = cell2mat(arrayfun(@(r) stage4.rep(r).centroidX_mm_101, valid(:), ...
                'UniformOutput', false));
            CY = cell2mat(arrayfun(@(r) stage4.rep(r).centroidY_mm_101, valid(:), ...
                'UniformOutput', false));

            mX = mean(mean(CX(:,inCVP),2,'omitnan'),'omitnan');
            mY = mean(mean(CY(:,inCVP),2,'omitnan'),'omitnan');

            % How far the centroid travels across the analysis window
            rX = mean(max(CX(:,inCVP),[],2) - min(CX(:,inCVP),[],2), 'omitnan');
            rY = mean(max(CY(:,inCVP),[],2) - min(CY(:,inCVP),[],2), 'omitnan');

            % Consistency of the mean centroid between repetitions
            sdReps = mean([std(mean(CX(:,inCVP),2,'omitnan'),'omitnan'), ...
                           std(mean(CY(:,inCVP),2,'omitnan'),'omitnan')]);
        end

        %% --- QC figure ---
        if CFG.saveFigures
            plot_stage4(stage4, x_norm_101, CFG.cvp_lo, CFG.cvp_hi, nB, tag, qcDir);
        end

        %% --- save ---
        % stage3 is not duplicated here, since it sits next to this file.
        [~, baseName] = fileparts(F(k).name);
        baseName = strrep(baseName, '_stage3_EMGnorm', '');
        outName = fullfile(F(k).folder, sprintf('%s_stage4_centroid_mm.mat', baseName));
        save(outName, 'stage4', '-v7.3');

        res(k,:) = {pName, cName, sName, nReps, nValid, mX, mY, rX, rY, ...
                    sdReps, nNonFinite, "ok"};

        fprintf(['[%3d/%3d] %-26s %d/%d reps | mean X %5.2f Y %5.2f mm | ' ...
                 'travel X %4.2f Y %4.2f | SD between reps %4.2f\n'], ...
            k, nF, tag, nValid, nReps, mX, mY, rX, rY, sdReps);

    catch ME
        res(k,:) = {pName, cName, sName, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                    NaN, string(ME.message)};
        fprintf('[%3d/%3d] %-26s FAILED: %s\n', k, nF, tag, ME.message);
    end
end

%% ---------------------------
%  (4) Summary
%% ---------------------------
ok = res.status == "ok";
fprintf('\n===== STAGE 4 SUMMARY =====\n');
fprintf('Processed %d | ok %d | failed %d\n', nF, nnz(ok), nnz(~ok));

if any(ok)
    fprintf('\nMean centroid over the %d-%d%% CVP, by condition (grid centre = 15.0, 35.0 mm):\n', ...
        CFG.cvp_lo, CFG.cvp_hi);
    fprintf('  %-8s %14s %14s\n', 'cond', 'X (mm)', 'Y (mm)');
    for cc = ["CON_75","ECC_75","CON_90","ECC_90"]
        sel = ok & res.condition == cc;
        if any(sel)
            fprintf('  %-8s %7.2f +/- %4.2f %7.2f +/- %4.2f   (n = %d)\n', cc, ...
                mean(res.meanX_mm(sel)), std(res.meanX_mm(sel)), ...
                mean(res.meanY_mm(sel)), std(res.meanY_mm(sel)), nnz(sel));
        end
    end

    fprintf('\nCON minus ECC difference in the mean centroid:\n');
    for I = ["75","90"]
        cs = ok & res.condition == "CON_"+I;
        es = ok & res.condition == "ECC_"+I;
        if any(cs) && any(es)
            fprintf('  %s%%: X %+6.3f mm | Y %+6.3f mm\n', I, ...
                mean(res.meanX_mm(cs)) - mean(res.meanX_mm(es)), ...
                mean(res.meanY_mm(cs)) - mean(res.meanY_mm(es)));
        end
    end

    fprintf('\nCentroid travel across the window : X %.2f mm | Y %.2f mm (median)\n', ...
        median(res.rangeX_mm(ok)), median(res.rangeY_mm(ok)));
    fprintf('SD of mean centroid between reps  : %.2f mm (median), max %.2f mm\n', ...
        median(res.sdAcrossReps_mm(ok)), max(res.sdAcrossReps_mm(ok)));

    % Expected to be zero. Anything else means channels dropped out during the
    % analysis window and the affected files should be inspected.
    totalBad = sum(res.nNonFiniteWeights(ok));
    if totalBad == 0
        fprintf('No non-finite weights encountered.\n');
    else
        fprintf('\nWARNING: %d non-finite weight samples treated as zero.\n', totalBad);
        badF = res(ok & res.nNonFiniteWeights > 0, :);
        disp(badF(:, {'participant','condition','set','nNonFiniteWeights'}));
    end

    miss = res(ok & res.nValid < res.nReps, :);
    if ~isempty(miss)
        fprintf('\nFiles with repetitions skipped:\n');
        disp(miss(:, {'participant','condition','set','nValid','nReps'}));
    end
end

if any(~ok)
    fprintf('\nFailures:\n');
    disp(res(~ok, {'participant','condition','set','status'}));
end

writetable(res, fullfile(ROOT,'stage4_summary.csv'));
fprintf('\nSaved summary: %s\n', fullfile(ROOT,'stage4_summary.csv'));
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

function plot_stage4(stage4, xAxis, lo, hi, nB, tag, qcDir)
% Centroid X and Y relative to the grid centre, with the mean weight per
% channel below so that a centroid computed from very low amplitude is visible.
    valid = find(arrayfun(@(s) s.is_valid, stage4.rep));
    if isempty(valid), return; end

    CX = cell2mat(arrayfun(@(r) stage4.rep(r).centroidX_centered_mm_101, valid(:), ...
        'UniformOutput', false));
    CY = cell2mat(arrayfun(@(r) stage4.rep(r).centroidY_centered_mm_101, valid(:), ...
        'UniformOutput', false));
    SW = cell2mat(arrayfun(@(r) stage4.rep(r).sumW_101, valid(:), ...
        'UniformOutput', false)) / nB;

    f = figure('Color','w','Visible','off','Position',[100 100 1000 800]);
    tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');

    a1 = nexttile; hold(a1,'on');
    plot(a1, xAxis, CX.', 'Color',[.75 .75 .95]);
    plot(a1, xAxis, mean(CX,1,'omitnan'), 'b','LineWidth',2);
    yline(a1, 0, '--'); xregion(a1, lo, hi, 'FaceAlpha',0.07);
    ylabel(a1,'Centroid X (mm)','FontWeight','bold');
    title(a1, tag + "  |  lateral (-) to medial (+), relative to the grid centre", ...
        'Interpreter','none','FontWeight','bold');
    set(a1,'XTickLabel',[]);

    a2 = nexttile; hold(a2,'on');
    plot(a2, xAxis, CY.', 'Color',[.95 .8 .7]);
    plot(a2, xAxis, mean(CY,1,'omitnan'), 'Color',[.85 .33 .1],'LineWidth',2);
    yline(a2, 0, '--'); xregion(a2, lo, hi, 'FaceAlpha',0.07);
    ylabel(a2,'Centroid Y (mm)','FontWeight','bold');
    set(a2,'XTickLabel',[]);

    a3 = nexttile; hold(a3,'on');
    plot(a3, xAxis, SW.', 'Color',[.8 .8 .8]);
    plot(a3, xAxis, mean(SW,1,'omitnan'), 'k','LineWidth',2);
    xregion(a3, lo, hi, 'FaceAlpha',0.07);
    xlabel(a3,'%CV','FontWeight','bold');
    ylabel(a3,'Mean %MVA per channel','FontWeight','bold');

    set([a1 a2 a3],'Box','off','TickDir','out','FontWeight','bold');
    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + ".png"), ...
        'Resolution',110);
    close(f);
end
