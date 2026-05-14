%% ===== Stage 6B: FULL torque scalar figure + MID(20-80%) full analysis =====
% INPUT:
%   groupData_stage6A.mat
%
% DESIGN:
%   Factor 1 = contractionType : CON vs ECC
%   Factor 2 = intensity       : 75 vs 90
%
% CONDITION ORDER (fixed):
%   1 = CON_75
%   2 = ECC_75
%   3 = CON_90
%   4 = ECC_90
%
% MAIN BIN:
%   all
%
% THIS VERSION:
%   - FULL TRACE: torque scalar figure only
%   - MID 20-80%: full scalar analysis + RM-ANOVA + post hoc + figures
%
% OUTPUT:
%   stage6B_stats_mid20to80.mat

clear; clc; close all;

%% ---------------------------
%  (1) Load groupData
%% ---------------------------
[f,p] = uigetfile('groupData_stage6A.mat', 'Select Stage6A groupData file');
if isequal(f,0)
    error('No groupData file selected.');
end

L = load(fullfile(p,f));
assert(isfield(L,'groupData'), 'Selected file does not contain "groupData".');
groupData = L.groupData;

%% ---------------------------
%  (2) Settings
%% ---------------------------
mainBin = 'all';
doPlots = true;
saveFigures = false;
alphaLevel = 0.05;
flipECCforPlots = true;

midWin = [20 80];

condNames = string(groupData.conditions(:));
expectedOrder = ["CON_75"; "ECC_75"; "CON_90"; "ECC_90"];
assert(all(condNames == expectedOrder), ...
    'Condition order mismatch. Expected [CON_75, ECC_75, CON_90, ECC_90].');

nSub = groupData.nSub;
nCond = groupData.nCond;
nT = groupData.nTime;
xAxis = groupData.xAxis(:).';

fprintf('\n===== Stage 6B loaded =====\n');
fprintf('Participants: %d\n', nSub);
fprintf('Conditions: %s\n', strjoin(cellstr(condNames.'), ', '));
fprintf('Main bin: %s\n', mainBin);
fprintf('Mid analysis window: %d-%d %%CV\n', midWin(1), midWin(2));

%% ---------------------------
%  (3) Build stats output struct
%% ---------------------------
stats6B = struct();
stats6B.meta = struct();
stats6B.meta.createdOn = datestr(now);
stats6B.meta.script = mfilename;
stats6B.meta.sourceFile = fullfile(p,f);
stats6B.meta.mainBin = mainBin;
stats6B.meta.nSub = nSub;
stats6B.meta.conditions = condNames;
stats6B.meta.alpha = alphaLevel;
stats6B.meta.midWin_percentCV = midWin;

stats6B.design = struct();
stats6B.design.factor1_name = 'contractionType';
stats6B.design.factor1_levels = {'CON','ECC'};
stats6B.design.factor2_name = 'intensity';
stats6B.design.factor2_levels = {'75','90'};

%% ---------------------------
%  (4A) Extract FULL scalar torque
%% ---------------------------
Y_torque_full = groupData.scalar.(mainBin).meanTorque;
assert(all(size(Y_torque_full)==[nSub nCond]), ...
    'groupData.scalar.(mainBin).meanTorque has unexpected size.');

stats6B.full = struct();
stats6B.full.meanTorque = Y_torque_full;


%% ---------------------------
%  (5) Extract time-series data
%% ---------------------------
TS = groupData.timeseries.(mainBin);

requiredFields = {'torque_101','emgMeanNorm_101','centroidXc_101','centroidYc_101','angle_101'};
for k = 1:numel(requiredFields)
    assert(isfield(TS, requiredFields{k}), ...
        'groupData.timeseries.%s.%s is missing.', mainBin, requiredFields{k});
end

T_torque = TS.torque_101;        % [nSub x nCond x nT]
T_emg    = TS.emgMeanNorm_101;   % [nSub x nCond x nT]
T_cx     = TS.centroidXc_101;    % [nSub x nCond x nT]
T_cy     = TS.centroidYc_101;    % [nSub x nCond x nT]
T_angle  = TS.angle_101;         % [nSub x nCond x nT]

assert(all(size(T_torque)==[nSub nCond nT]), 'Unexpected size for torque_101.');
assert(all(size(T_emg)==[nSub nCond nT]),    'Unexpected size for emgMeanNorm_101.');
assert(all(size(T_cx)==[nSub nCond nT]),     'Unexpected size for centroidXc_101.');
assert(all(size(T_cy)==[nSub nCond nT]),     'Unexpected size for centroidYc_101.');
assert(all(size(T_angle)==[nSub nCond nT]),  'Unexpected size for angle_101.');

%% ---------------------------
%  (6) Build MID 20-80% scalars from time-series
%% ---------------------------
midMask = xAxis >= midWin(1) & xAxis <= midWin(2);
assert(any(midMask), 'No samples found inside the requested mid window.');

Y_torque_mid = mean(T_torque(:,:,midMask), 3, 'omitnan');
Y_emg_mid    = mean(T_emg(:,:,midMask),    3, 'omitnan');
Y_cx_mid     = mean(T_cx(:,:,midMask),     3, 'omitnan');
Y_cy_mid     = mean(T_cy(:,:,midMask),     3, 'omitnan');

assert(all(size(Y_torque_mid)==[nSub nCond]));
assert(all(size(Y_emg_mid)==[nSub nCond]));
assert(all(size(Y_cx_mid)==[nSub nCond]));
assert(all(size(Y_cy_mid)==[nSub nCond]));

stats6B.mid20to80 = struct();
stats6B.mid20to80.mask = midMask;
stats6B.mid20to80.xAxis = xAxis(midMask);
stats6B.mid20to80.Y_torque = Y_torque_mid;
stats6B.mid20to80.Y_emg    = Y_emg_mid;
stats6B.mid20to80.Y_cx     = Y_cx_mid;
stats6B.mid20to80.Y_cy     = Y_cy_mid;




%% ---------------------------
%  (6B) MID 20-80% TORQUE MISMATCH + EMG NORMALIZED TO TORQUE
%% ---------------------------

fprintf('\n===== MID 20-80%% TORQUE MISMATCH + EMG/TORQUE ANALYSIS =====\n');

% ---------------------------------------
% (A) TORQUE MISMATCH (ECC - CON) in MID 20-80%
% ---------------------------------------
d75 = Y_torque_mid(:,2) - Y_torque_mid(:,1);
d90 = Y_torque_mid(:,4) - Y_torque_mid(:,3);

% Relative (%), referenced to CON
rel75 = (d75 ./ Y_torque_mid(:,1)) * 100;
rel90 = (d90 ./ Y_torque_mid(:,3)) * 100;

fprintf('\n--- Torque mismatch (MID 20-80%%) ---\n');
fprintf('75%%: %.3f ± %.3f Nm (%.2f ± %.2f %%)\n', ...
    mean(d75,'omitnan'), std(d75,0,'omitnan'), ...
    mean(rel75,'omitnan'), std(rel75,0,'omitnan'));

fprintf('90%%: %.3f ± %.3f Nm (%.2f ± %.2f %%)\n', ...
    mean(d90,'omitnan'), std(d90,0,'omitnan'), ...
    mean(rel90,'omitnan'), std(rel90,0,'omitnan'));

% Effect size (paired dz)
dz75 = mean(d75,'omitnan') / std(d75,0,'omitnan');
dz90 = mean(d90,'omitnan') / std(d90,0,'omitnan');

fprintf('Effect size dz (75%%): %.3f\n', dz75);
fprintf('Effect size dz (90%%): %.3f\n', dz90);

% 95% CI
tcrit = tinv(0.975, nSub-1);
ci75 = mean(d75,'omitnan') + [-1 1] * tcrit * std(d75,0,'omitnan') / sqrt(nSub);
ci90 = mean(d90,'omitnan') + [-1 1] * tcrit * std(d90,0,'omitnan') / sqrt(nSub);

fprintf('95%% CI 75%%: [%.3f  %.3f]\n', ci75(1), ci75(2));
fprintf('95%% CI 90%%: [%.3f  %.3f]\n', ci90(1), ci90(2));

% Save
stats6B.mid20to80.mismatch = struct();
stats6B.mid20to80.mismatch.d75 = d75;
stats6B.mid20to80.mismatch.d90 = d90;
stats6B.mid20to80.mismatch.rel75 = rel75;
stats6B.mid20to80.mismatch.rel90 = rel90;
stats6B.mid20to80.mismatch.dz75 = dz75;
stats6B.mid20to80.mismatch.dz90 = dz90;
stats6B.mid20to80.mismatch.ci75 = ci75;
stats6B.mid20to80.mismatch.ci90 = ci90;

% ---------------------------------------
% (B) EMG NORMALIZED TO TORQUE in MID 20-80%
% ---------------------------------------
EMG_perTorque_mid = Y_emg_mid ./ Y_torque_mid;

fprintf('\n--- EMG per Torque (MID 20-80%%) ---\n');
disp(mean(EMG_perTorque_mid,1,'omitnan'))

% Save
stats6B.mid20to80.emgPerTorque = struct();
stats6B.mid20to80.emgPerTorque.values = EMG_perTorque_mid;

% ---------------------------------------
% (C) RM-ANOVA ON EMG/TORQUE
% ---------------------------------------
stats6B.mid20to80.scalar.meanEMG_perTorque = run_rm2x2( ...
    EMG_perTorque_mid, groupData.participants, 'meanEMG_perTorque_mid20to80');

fprintf('\n--- EMG per Torque RM-ANOVA (MID 20-80%%) ---\n');
print_rm_result(stats6B.mid20to80.scalar.meanEMG_perTorque, 'EMG per Torque | mid 20-80%');

% ---------------------------------------
% (D) Pairwise follow-up tests for EMG/Torque
% ---------------------------------------
stats6B.mid20to80.posthoc.meanEMG_perTorque.CONvsECC_75 = ...
    paired_compare(EMG_perTorque_mid(:,1), EMG_perTorque_mid(:,2), ...
    'EMG/Torque mid20to80: CON75 vs ECC75');

stats6B.mid20to80.posthoc.meanEMG_perTorque.CONvsECC_90 = ...
    paired_compare(EMG_perTorque_mid(:,3), EMG_perTorque_mid(:,4), ...
    'EMG/Torque mid20to80: CON90 vs ECC90');

fprintf('\n--- EMG/Torque means (MID 20-80%%) ---\n');
fprintf('75%%: CON = %.3f, ECC = %.3f\n', ...
    mean(EMG_perTorque_mid(:,1),'omitnan'), mean(EMG_perTorque_mid(:,2),'omitnan'));

fprintf('90%%: CON = %.3f, ECC = %.3f\n', ...
    mean(EMG_perTorque_mid(:,3),'omitnan'), mean(EMG_perTorque_mid(:,4),'omitnan'));
%% ---------------------------
%  (7) MID 20-80% Scalar 2x2 RM-ANOVA
%% ---------------------------
stats6B.mid20to80.scalar = struct();

stats6B.mid20to80.scalar.meanTorque    = run_rm2x2(Y_torque_mid, groupData.participants, 'meanTorque_mid20to80');
stats6B.mid20to80.scalar.meanEMG       = run_rm2x2(Y_emg_mid,    groupData.participants, 'meanEMG_mid20to80');
stats6B.mid20to80.scalar.meanCentroidX = run_rm2x2(Y_cx_mid,     groupData.participants, 'meanCentroidX_mid20to80');
stats6B.mid20to80.scalar.meanCentroidY = run_rm2x2(Y_cy_mid,     groupData.participants, 'meanCentroidY_mid20to80');

%% ---------------------------
%  (8) MID 20-80% Pairwise follow-up tests
%% ---------------------------
stats6B.mid20to80.posthoc = struct();

stats6B.mid20to80.posthoc.meanTorque.CONvsECC_75    = paired_compare(Y_torque_mid(:,1), Y_torque_mid(:,2), 'Torque mid20to80: CON75 vs ECC75');
stats6B.mid20to80.posthoc.meanTorque.CONvsECC_90    = paired_compare(Y_torque_mid(:,3), Y_torque_mid(:,4), 'Torque mid20to80: CON90 vs ECC90');

stats6B.mid20to80.posthoc.meanEMG.CONvsECC_75       = paired_compare(Y_emg_mid(:,1), Y_emg_mid(:,2), 'EMG mid20to80: CON75 vs ECC75');
stats6B.mid20to80.posthoc.meanEMG.CONvsECC_90       = paired_compare(Y_emg_mid(:,3), Y_emg_mid(:,4), 'EMG mid20to80: CON90 vs ECC90');

stats6B.mid20to80.posthoc.meanCentroidX.CONvsECC_75 = paired_compare(Y_cx_mid(:,1), Y_cx_mid(:,2), 'CentroidX mid20to80: CON75 vs ECC75');
stats6B.mid20to80.posthoc.meanCentroidX.CONvsECC_90 = paired_compare(Y_cx_mid(:,3), Y_cx_mid(:,4), 'CentroidX mid20to80: CON90 vs ECC90');

stats6B.mid20to80.posthoc.meanCentroidY.CONvsECC_75 = paired_compare(Y_cy_mid(:,1), Y_cy_mid(:,2), 'CentroidY mid20to80: CON75 vs ECC75');
stats6B.mid20to80.posthoc.meanCentroidY.CONvsECC_90 = paired_compare(Y_cy_mid(:,3), Y_cy_mid(:,4), 'CentroidY mid20to80: CON90 vs ECC90');

%% ---------------------------
%  (9) Print results
%% ---------------------------
fprintf('\n===== FULL TRACE TORQUE SUMMARY =====\n');
fprintf('Condition means [CON75 ECC75 CON90 ECC90]:\n');
disp(mean(Y_torque_full,1,'omitnan'))
fprintf('Condition SDs   [CON75 ECC75 CON90 ECC90]:\n');
disp(std(Y_torque_full,0,1,'omitnan'))

fprintf('\n===== MID 20-80%% SCALAR RESULTS =====\n');
print_rm_result(stats6B.mid20to80.scalar.meanTorque,    'Mean Torque | mid 20-80%');
print_rm_result(stats6B.mid20to80.scalar.meanEMG,       'Mean EMG | mid 20-80%');
print_rm_result(stats6B.mid20to80.scalar.meanCentroidX, 'Mean Centroid X | mid 20-80%');
print_rm_result(stats6B.mid20to80.scalar.meanCentroidY, 'Mean Centroid Y | mid 20-80%');

%% ---------------------------
%  (10) Figures
%% ---------------------------
if doPlots
    % =========================
    % FULL TRACE torque scalar figure only
    % =========================
    figTorqueFull = figure('Color','w','Name','Stage6B FULL Scalar | Mean Torque');
    axTF = axes(figTorqueFull);
    plot_scalar_points_panel(axTF, Y_torque_full, 'Matched torque | FULL trace', 'Torque (Nm)');

    % =========================
    % MID 20-80% scalar figures
    % =========================
    figTorqueMid = figure('Color','w','Name','Stage6B MID Scalar | Mean Torque');
    axTM = axes(figTorqueMid);
    plot_scalar_points_panel(axTM, Y_torque_mid, 'Matched torque | MID 20-80% CV', 'Torque (Nm)');

    figEMGMid = figure('Color','w','Name','Stage6B MID Scalar | Mean EMG');
    axEM = axes(figEMGMid);
    plot_scalar_points_panel(axEM, Y_emg_mid, 'Mean EMG-RMS | MID 20-80% CV', 'EMG-RMS (%MVC)');

    figCXMid = figure('Color','w','Name','Stage6B MID Scalar | Mean Centroid X');
    axCXM = axes(figCXMid);
    plot_scalar_points_panel(axCXM, Y_cx_mid, 'Mean Centroid (Lateral-Medial) | MID 20-80% CV', 'Distance (mm)');

    figCYMid = figure('Color','w','Name','Stage6B MID Scalar | Mean Centroid Y');
    axCYM = axes(figCYMid);
    plot_scalar_points_panel(axCYM, Y_cy_mid, 'Mean Centroid (Proximal-Distal) | MID 20-80% CV', 'Distance (mm)');

    % =========================
    % MID 20-80% Combined CV figures
    % =========================
    fig75 = figure('Color','w','Name','Stage6B MID Combined CV Figure | 75%');
    tl75 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    plot_condition_pair(nexttile(tl75), xAxis, ...
        squeeze(T_torque(:,1,:)), maybe_flip_timeseries(squeeze(T_torque(:,2,:)), flipECCforPlots), ...
        'CON','ECC','Torque (Nm)','',false,false,midWin);

    plot_condition_pair(nexttile(tl75), xAxis, ...
        squeeze(T_angle(:,1,:)), maybe_flip_timeseries(squeeze(T_angle(:,2,:)), flipECCforPlots), ...
        'CON','ECC','Angle (deg)','',false,false,midWin);

    plot_condition_pair(nexttile(tl75), xAxis, ...
        squeeze(T_emg(:,1,:)), maybe_flip_timeseries(squeeze(T_emg(:,2,:)), flipECCforPlots), ...
        'CON','ECC','Mean EMG (%MVC)','',false,true,midWin);

    sgtitle(sprintf('CON vs ECC at 75%% MVC | MID scalar window %d-%d%% CV', midWin(1), midWin(2)), ...
        'FontWeight','bold');

    fig90 = figure('Color','w','Name','Stage6B MID Combined CV Figure | 90%');
    tl90 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    plot_condition_pair(nexttile(tl90), xAxis, ...
        squeeze(T_torque(:,3,:)), maybe_flip_timeseries(squeeze(T_torque(:,4,:)), flipECCforPlots), ...
        'CON','ECC','Torque (Nm)','',false,false,midWin);

    plot_condition_pair(nexttile(tl90), xAxis, ...
        squeeze(T_angle(:,3,:)), maybe_flip_timeseries(squeeze(T_angle(:,4,:)), flipECCforPlots), ...
        'CON','ECC','Angle (deg)','',false,false,midWin);

    plot_condition_pair(nexttile(tl90), xAxis, ...
        squeeze(T_emg(:,3,:)), maybe_flip_timeseries(squeeze(T_emg(:,4,:)), flipECCforPlots), ...
        'CON','ECC','Mean EMG (%MVC)','',false,true,midWin);

    sgtitle(sprintf('CON vs ECC at 90%% MVC | MID scalar window %d-%d%% CV', midWin(1), midWin(2)), ...
        'FontWeight','bold');

    % =========================
    % MID 20-80% Centroid time series
    % =========================
    figCX = figure('Color','w','Name','Stage6B MID Centroid X Time Series');
    tlCX = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    plot_condition_pair(nexttile(tlCX), xAxis, ...
        squeeze(T_cx(:,1,:)), maybe_flip_timeseries(squeeze(T_cx(:,2,:)), flipECCforPlots), ...
        'CON 75', 'ECC 75', 'Centroid X (mm)', 'Centroid X | 75%', true, true, midWin);

    plot_condition_pair(nexttile(tlCX), xAxis, ...
        squeeze(T_cx(:,3,:)), maybe_flip_timeseries(squeeze(T_cx(:,4,:)), flipECCforPlots), ...
        'CON 90', 'ECC 90', 'Centroid X (mm)', 'Centroid X | 90%', true, true, midWin);

    sgtitle(sprintf('Centroid Lateral-Medial | MID scalar window %d-%d%% CV', midWin(1), midWin(2)), ...
        'FontWeight','bold');

    figCY = figure('Color','w','Name','Stage6B MID Centroid Y Time Series');
    tlCY = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    plot_condition_pair(nexttile(tlCY), xAxis, ...
        squeeze(T_cy(:,1,:)), maybe_flip_timeseries(squeeze(T_cy(:,2,:)), flipECCforPlots), ...
        'CON 75', 'ECC 75', 'Centroid Y (mm)', 'Centroid Y | 75%', true, true, midWin);

    plot_condition_pair(nexttile(tlCY), xAxis, ...
        squeeze(T_cy(:,3,:)), maybe_flip_timeseries(squeeze(T_cy(:,4,:)), flipECCforPlots), ...
        'CON 90', 'ECC 90', 'Centroid Y (mm)', 'Centroid Y | 90%', true, true, midWin);

    sgtitle(sprintf('Centroid Proximal-Distal | MID scalar window %d-%d%% CV', midWin(1), midWin(2)), ...
        'FontWeight','bold');



    % =========================
    % EMG per Torque figure
    % =========================
    figEMGnorm = figure('Color','w','Name','Stage6B Scalar | EMG per Torque');
    axN = axes(figEMGnorm);
    plot_scalar_points_panel(axN, EMG_perTorque_mid, ...
    'EMG normalized to torque | MID 20-80% CV', 'EMG / Torque');
    % =========================
    % Optional export
    % =========================
    if saveFigures
        outFigDir = uigetdir(pwd, 'Select folder to save Stage6B figures'); %#ok<*UNRCH>
        if ~isequal(outFigDir,0)
            exportgraphics(figTorqueFull, fullfile(outFigDir, 'Stage6B_FULL_scalar_meanTorque.png'), 'Resolution', 300);

            exportgraphics(figTorqueMid, fullfile(outFigDir, 'Stage6B_MID20to80_scalar_meanTorque.png'), 'Resolution', 300);
            exportgraphics(figEMGMid,   fullfile(outFigDir, 'Stage6B_MID20to80_scalar_meanEMG.png'), 'Resolution', 300);
            exportgraphics(figCXMid,    fullfile(outFigDir, 'Stage6B_MID20to80_scalar_meanCentroidX.png'), 'Resolution', 300);
            exportgraphics(figCYMid,    fullfile(outFigDir, 'Stage6B_MID20to80_scalar_meanCentroidY.png'), 'Resolution', 300);

            exportgraphics(fig75, fullfile(outFigDir, 'Stage6B_MID20to80_CV_75.png'), 'Resolution', 300);
            exportgraphics(fig90, fullfile(outFigDir, 'Stage6B_MID20to80_CV_90.png'), 'Resolution', 300);

            exportgraphics(figCX, fullfile(outFigDir, 'Stage6B_MID20to80_centroidX_timeseries.png'), 'Resolution', 300);
            exportgraphics(figCY, fullfile(outFigDir, 'Stage6B_MID20to80_centroidY_timeseries.png'), 'Resolution', 300);
        end
    end
end

%% ---------------------------
%  (11) Save stats
%% ---------------------------
[outFile, outPath] = uiputfile('stage6B_stats_mid20to80.mat', 'Save Stage6B stats as');
if isequal(outFile,0)
    outPath = pwd;
    outFile = 'stage6B_stats_mid20to80.mat';
end

save(fullfile(outPath, outFile), 'stats6B', '-v7.3');
fprintf('\n✅ Saved Stage 6B stats:\n%s\n', fullfile(outPath, outFile));

%% =========================================================
%  Local functions
%% =========================================================
function S = run_rm2x2(Y, participants, varLabel)
% Y is [nSub x 4]
% Order:
% 1=CON_75, 2=ECC_75, 3=CON_90, 4=ECC_90

    T = table( ...
        string(participants(:)), ...
        Y(:,1), Y(:,2), Y(:,3), Y(:,4), ...
        'VariableNames', {'Subject','CON_75','ECC_75','CON_90','ECC_90'});

    within = table( ...
        categorical({'CON';'ECC';'CON';'ECC'}), ...
        categorical({'75';'75';'90';'90'}), ...
        'VariableNames', {'Contraction','Intensity'});

    rm = fitrm(T, 'CON_75-ECC_90 ~ 1', 'WithinDesign', within);
    ranovatbl = ranova(rm, 'WithinModel', 'Contraction*Intensity');

    S = struct();
    S.label = varLabel;
    S.data = Y;
    S.subjects = string(participants(:));
    S.table = T;
    S.within = within;
    S.rm = rm;
    S.ranova = ranovatbl;

    S.p_contraction = get_p_from_ranova(ranovatbl, 'Contraction');
    S.p_intensity   = get_p_from_ranova(ranovatbl, 'Intensity');
    S.p_interaction = get_p_from_ranova(ranovatbl, 'Contraction:Intensity');

    S.mean = mean(Y,1,'omitnan');
    S.sd   = std(Y,0,1,'omitnan');
end

function p = get_p_from_ranova(tbl, effectLabel)
    p = NaN;
    try
        rn = string(tbl.Properties.RowNames);
        idx = find(contains(rn, effectLabel), 1, 'first');

        if isempty(idx)
            return
        end

        if any(strcmp(tbl.Properties.VariableNames, 'pValue')) && isfinite(tbl.pValue(idx))
            p = tbl.pValue(idx);
        elseif any(strcmp(tbl.Properties.VariableNames, 'pValueGG')) && isfinite(tbl.pValueGG(idx))
            p = tbl.pValueGG(idx);
        end
    catch
        p = NaN;
    end
end

function print_rm_result(S, labelText)
    fprintf('\n--- %s ---\n', labelText);
    fprintf('Condition means [CON75 ECC75 CON90 ECC90]:\n');
    disp(S.mean)
    fprintf('Condition SDs   [CON75 ECC75 CON90 ECC90]:\n');
    disp(S.sd)

    fprintf('p(Contraction) = %.4f\n', S.p_contraction);
    fprintf('p(Intensity)   = %.4f\n', S.p_intensity);
    fprintf('p(Interaction) = %.4f\n', S.p_interaction);
end

function H = paired_compare(y1, y2, labelText)
    [~, p, ~, stats] = ttest(y1, y2);

    H = struct();
    H.label = labelText;
    H.mean1 = mean(y1,'omitnan');
    H.mean2 = mean(y2,'omitnan');
    H.sd1   = std(y1,0,'omitnan');
    H.sd2   = std(y2,0,'omitnan');
    H.p     = p;
    H.stats = stats;
end

function plot_scalar_points_panel(ax, Y, ttl, ylab)
    hold(ax,'on');

    xPos = [1 2 4 5];

    yMin0 = min(Y(:),[],'omitnan');
    yMax0 = max(Y(:),[],'omitnan');
    yRange0 = yMax0 - yMin0;
    if yRange0 == 0
        yRange0 = 1;
    end
    pad = 0.12 * yRange0;
    yLo = yMin0 - pad;
    yHi = yMax0 + pad;

    nSub = size(Y,1);
    subjColors = distinguishable_colors(nSub);
    subjColors = subjColors(1:nSub,:);

    lineColor    = [0.75 0.75 0.75];
    lineWidthInd = 0.5;
    markerSize   = 28;

    for i = 1:nSub
        yi = Y(i,:);
        thisColor = subjColors(i,:);

        if ~isnan(yi(1)) && ~isnan(yi(2))
            plot(ax, xPos(1:2), yi(1:2), '-', 'Color', lineColor, 'LineWidth', lineWidthInd);
        end

        if ~isnan(yi(3)) && ~isnan(yi(4))
            plot(ax, xPos(3:4), yi(3:4), '-', 'Color', lineColor, 'LineWidth', lineWidthInd);
        end

        for c = 1:4
            if ~isnan(yi(c))
                scatter(ax, xPos(c), yi(c), markerSize, ...
                    'MarkerFaceColor', thisColor, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.4);
            end
        end
    end

    m = mean(Y,1,'omitnan');
    s = std(Y,0,1,'omitnan');

    for c = 1:4
        if ~isnan(m(c)) && ~isnan(s(c))
            plot(ax, [xPos(c) xPos(c)], [m(c)-s(c), m(c)+s(c)], ...
                '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0);

            capW = 0.10;
            plot(ax, [xPos(c)-capW xPos(c)+capW], [m(c)-s(c) m(c)-s(c)], ...
                '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0);
            plot(ax, [xPos(c)-capW xPos(c)+capW], [m(c)+s(c) m(c)+s(c)], ...
                '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0);
        end
    end

    plot(ax, xPos(1:2), m(1:2), '--', 'Color', 'k', 'LineWidth', 1.2);
    plot(ax, xPos(3:4), m(3:4), '--', 'Color', 'k', 'LineWidth', 1.2);

    scatter(ax, xPos, m, 42, ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5);

    xlim(ax, [0.5 5.5]);
    ylim(ax, [yLo yHi]);

    set(ax, ...
        'XTick', xPos, ...
        'XTickLabel', {'CON75','ECC75','CON90','ECC90'}, ...
        'FontSize', 10, ...
        'FontWeight','bold', ...
        'LineWidth',1.5, ...
        'TickDir','out');

    ylabel(ax, ylab,'FontWeight','bold','FontSize',10);
    title(ax, ttl,'Interpreter','none','FontWeight','bold','FontSize',10);
    box(ax, 'off');

    if yLo < 0 && yHi > 0
        yline(ax, 0, ':k', 'LineWidth', 1.0);
    end
end

function Yout = maybe_flip_timeseries(Yin, doFlip)
    Yout = Yin;
    if doFlip
        Yout = fliplr(Yin);
    end
end

function plot_condition_pair(ax, xAxis, Y1, Y2, name1, name2, ylab, ttl, addZeroLine, showXlabel, winCV)

    if nargin < 9 || isempty(addZeroLine)
        addZeroLine = false;
    end

    if nargin < 10 || isempty(showXlabel)
        showXlabel = true;
    end

    if nargin < 11 || isempty(winCV)
        winCV = [];
    end

    hold(ax,'on');

    m1 = mean(Y1,1,'omitnan');
    s1 = std(Y1,0,1,'omitnan');

    m2 = mean(Y2,1,'omitnan');
    s2 = std(Y2,0,1,'omitnan');

    colCON = [0.00 0.4470 0.7410];
    colECC = [0.4660 0.6740 0.1880];

    fill(ax, [xAxis fliplr(xAxis)], [m1-s1 fliplr(m1+s1)], ...
        colCON, 'FaceAlpha',0.10, 'EdgeColor','none', 'HandleVisibility','off');

    fill(ax, [xAxis fliplr(xAxis)], [m2-s2 fliplr(m2+s2)], ...
        colECC, 'FaceAlpha',0.10, 'EdgeColor','none', 'HandleVisibility','off');

    hCON = plot(ax, xAxis, m1, 'Color', colCON, 'LineWidth', 2.2);
    hECC = plot(ax, xAxis, m2, 'Color', colECC, 'LineWidth', 2.2);

    if ~isempty(winCV) && numel(winCV)==2
        xline(ax, winCV(1), '--k', 'LineWidth', 1.0, 'HandleVisibility','off');
        xline(ax, winCV(2), '--k', 'LineWidth', 1.0, 'HandleVisibility','off');
    end

    if addZeroLine
        yline(ax, 0, '--k', 'LineWidth', 0.9);
    end

    if showXlabel
        xlabel(ax, 'Contraction Cycle (CV%)', 'FontWeight','bold');
    else
        set(ax,'XTickLabel',[])
    end

    ylabel(ax, ylab, 'FontWeight','bold');

    if ~isempty(ttl)
        title(ax, ttl, 'Interpreter','none', 'FontWeight','bold');
    end

    legend(ax, [hCON hECC], {char(name1), char(name2)}, ...
        'Box','off','Location','best','Interpreter','none');

    box(ax,'off');
    set(ax, 'TickDir','out', 'LineWidth',1.4, 'FontWeight','bold');

    if contains(ylab, 'Centroid Y', 'IgnoreCase', true)
        ylim(ax, [-5 5]);
        yticks(ax, -5:1:5);
    end
end