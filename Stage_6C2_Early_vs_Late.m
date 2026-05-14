%% ===== Stage 6C2: Early vs Late MID 20-80% Analysis (3-way RM ANOVA + post hoc paired SPM) =====
% PURPOSE:
%   Analyze EARLY vs LATE effects within the MID 20-80% CV phase for:
%       - torque
%       - mean EMG
%
% WITHIN-SUBJECT FACTORS:
%   A) Contraction type : CON vs ECC
%   B) Intensity        : 75% vs 90%
%   C) Time-bin         : Early vs Late
%
% INPUT:
%   groupData_stage6A.mat
%
% OUTPUT:
%   stage6C2_earlyLate_mid20to80.mat
%
% NOTES:
%   - Restricts scalar and time-series analyses to 20-80% CV
%   - Uses same ECC flipping convention as Stage 6C1
%   - Uses corrected cluster coordinate mapping for cropped domain


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
alphaLevel = 0.05;
twoTailed = true;
interpClusters = true;
flipECCforAnalysis = true;
doPlots = true;
saveFigures = false;

midWin = [20 80];

condNames = string(groupData.conditions(:));
expectedOrder = ["CON_75"; "ECC_75"; "CON_90"; "ECC_90"];
assert(all(condNames == expectedOrder), ...
    'Condition order mismatch. Expected [CON_75, ECC_75, CON_90, ECC_90].');

nSub  = groupData.nSub;
nCond = groupData.nCond;
xAxis_full = groupData.xAxis(:).';

midMask = xAxis_full >= midWin(1) & xAxis_full <= midWin(2);
assert(any(midMask), 'No samples found inside the requested mid window.');

xAxis = xAxis_full(midMask);
nT = numel(xAxis);

fprintf('\n===== Stage 6C2 MID loaded =====\n');
fprintf('Participants: %d\n', nSub);
fprintf('Conditions: %s\n', strjoin(cellstr(condNames.'), ', '));
fprintf('Alpha: %.3f\n', alphaLevel);
fprintf('Flip ECC for analysis: %d\n', flipECCforAnalysis);
fprintf('Analysis window: %d-%d %%CV\n', midWin(1), midWin(2));

%% ---------------------------
%  (3) Check bins + spm1d
%% ---------------------------
reqBins = {'early','late'};
for i = 1:numel(reqBins)
    assert(isfield(groupData.scalar, reqBins{i}), 'Missing scalar.%s', reqBins{i});
    assert(isfield(groupData.timeseries, reqBins{i}), 'Missing timeseries.%s', reqBins{i});
end

spmFun = which('spm1d.stats.anova3rm');
assert(~isempty(spmFun), ...
    ['spm1d toolbox not found on MATLAB path. ', ...
     'Add spm1d first, e.g. addpath(genpath(''path_to_spm1d''));']);

fprintf('SPM1D found at:\n%s\n', spmFun);

%% ---------------------------
%  (4) Extract scalar outcomes from MID 20-80%
%% ---------------------------
% Build from timeseries so scalar window matches Stage 6B / 6C1 exactly

T_torque_early_full = groupData.timeseries.early.torque_101;       % [nSub x 4 x 101]
T_torque_late_full  = groupData.timeseries.late.torque_101;

T_emg_early_full    = groupData.timeseries.early.emgMeanNorm_101;  % [nSub x 4 x 101]
T_emg_late_full     = groupData.timeseries.late.emgMeanNorm_101;

assert(all(size(T_torque_early_full)==[nSub nCond numel(xAxis_full)]));
assert(all(size(T_torque_late_full) ==[nSub nCond numel(xAxis_full)]));
assert(all(size(T_emg_early_full)   ==[nSub nCond numel(xAxis_full)]));
assert(all(size(T_emg_late_full)    ==[nSub nCond numel(xAxis_full)]));

Y_torque_early = mean(T_torque_early_full(:,:,midMask), 3, 'omitnan');   % [nSub x 4]
Y_torque_late  = mean(T_torque_late_full(:,:,midMask),  3, 'omitnan');

Y_emg_early    = mean(T_emg_early_full(:,:,midMask),    3, 'omitnan');
Y_emg_late     = mean(T_emg_late_full(:,:,midMask),     3, 'omitnan');

assert(all(size(Y_torque_early)==[nSub nCond]));
assert(all(size(Y_torque_late) ==[nSub nCond]));
assert(all(size(Y_emg_early)   ==[nSub nCond]));
assert(all(size(Y_emg_late)    ==[nSub nCond]));

%% ---------------------------
%  (5) Extract MID 20-80% time series
%% ---------------------------
T_torque_early = T_torque_early_full(:,:,midMask);   % [nSub x 4 x nT]
T_torque_late  = T_torque_late_full(:,:,midMask);

T_emg_early    = T_emg_early_full(:,:,midMask);
T_emg_late     = T_emg_late_full(:,:,midMask);

assert(all(size(T_torque_early)==[nSub nCond nT]));
assert(all(size(T_torque_late) ==[nSub nCond nT]));
assert(all(size(T_emg_early)   ==[nSub nCond nT]));
assert(all(size(T_emg_late)    ==[nSub nCond nT]));

%% ---------------------------
%  (5b) Optionally flip ECC conditions
%% ---------------------------
if flipECCforAnalysis
    % Torque
    T_torque_early(:,2,:) = flip(T_torque_early(:,2,:), 3);
    T_torque_early(:,4,:) = flip(T_torque_early(:,4,:), 3);
    T_torque_late(:,2,:)  = flip(T_torque_late(:,2,:),  3);
    T_torque_late(:,4,:)  = flip(T_torque_late(:,4,:),  3);

    % EMG
    T_emg_early(:,2,:) = flip(T_emg_early(:,2,:), 3);
    T_emg_early(:,4,:) = flip(T_emg_early(:,4,:), 3);
    T_emg_late(:,2,:)  = flip(T_emg_late(:,2,:),  3);
    T_emg_late(:,4,:)  = flip(T_emg_late(:,4,:),  3);
end

%% ---------------------------
%  (6) Build output struct
%% ---------------------------
fatigue6C2 = struct();
fatigue6C2.meta = struct();
fatigue6C2.meta.createdOn = datestr(now);
fatigue6C2.meta.script = mfilename;
fatigue6C2.meta.sourceFile = fullfile(p,f);
fatigue6C2.meta.nSub = nSub;
fatigue6C2.meta.nCond = nCond;
fatigue6C2.meta.nTime = nT;
fatigue6C2.meta.xAxis = xAxis;
fatigue6C2.meta.alpha = alphaLevel;
fatigue6C2.meta.conditions = condNames;
fatigue6C2.meta.flipECCforAnalysis = flipECCforAnalysis;
fatigue6C2.meta.analysisWindow = midWin;
fatigue6C2.meta.design = '2x2x2 repeated-measures ANOVA';
fatigue6C2.meta.factorA = 'Contraction';
fatigue6C2.meta.factorALevels = {'CON','ECC'};
fatigue6C2.meta.factorB = 'Intensity';
fatigue6C2.meta.factorBLevels = {'75','90'};
fatigue6C2.meta.factorC = 'Time-bin';
fatigue6C2.meta.factorCLevels = {'Early','Late'};

%% ---------------------------
%  (7) Scalar paired tests
%% ---------------------------
fatigue6C2.scalar = struct();

for c = 1:nCond
    condLabel = char(condNames(c));

    fatigue6C2.scalar.torque.(condLabel) = paired_scalar_earlyLate( ...
        Y_torque_early(:,c), Y_torque_late(:,c), ...
        sprintf('Torque | %s | Early vs Late | mid20to80', condLabel));

    fatigue6C2.scalar.emg.(condLabel) = paired_scalar_earlyLate( ...
        Y_emg_early(:,c), Y_emg_late(:,c), ...
        sprintf('Mean EMG | %s | Early vs Late | mid20to80', condLabel));
end

%% ---------------------------
%  (8) Omnibus SPM RM-ANOVA (3-way)
%% ---------------------------
fatigue6C2.spmRM = struct();

fatigue6C2.spmRM.torque = run_spm_anova3rm_earlyLate( ...
    T_torque_early, T_torque_late, xAxis, alphaLevel, interpClusters, ...
    'Torque', 'Nm');

fatigue6C2.spmRM.emg = run_spm_anova3rm_earlyLate( ...
    T_emg_early, T_emg_late, xAxis, alphaLevel, interpClusters, ...
    'Mean EMG', '%MVC');

%% ---------------------------
%  (9) Post hoc paired SPM tests: Early vs Late within each condition
%% ---------------------------
fatigue6C2.posthoc = struct();

for c = 1:nCond
    condLabel = char(condNames(c));

    fatigue6C2.posthoc.torque.(condLabel) = run_spm_pair_earlyLate( ...
        squeeze(T_torque_early(:,c,:)), ...
        squeeze(T_torque_late(:,c,:)), ...
        xAxis, alphaLevel, twoTailed, interpClusters, ...
        'Torque', condLabel, 'Nm');

    fatigue6C2.posthoc.emg.(condLabel) = run_spm_pair_earlyLate( ...
        squeeze(T_emg_early(:,c,:)), ...
        squeeze(T_emg_late(:,c,:)), ...
        xAxis, alphaLevel, twoTailed, interpClusters, ...
        'Mean EMG', condLabel, '%MVC');
end

%% ---------------------------
%  (10) Print summary
%% ---------------------------
fprintf('\n===== Stage 6C2 MID SCALAR SUMMARY =====\n');

for c = 1:nCond
    condLabel = char(condNames(c));

    S1 = fatigue6C2.scalar.torque.(condLabel);
    fprintf('\n--- Torque | %s | Early vs Late | MID 20-80%% ---\n', condLabel);
    fprintf('Early mean ± SD: %.3f ± %.3f\n', S1.meanEarly, S1.sdEarly);
    fprintf('Late  mean ± SD: %.3f ± %.3f\n', S1.meanLate,  S1.sdLate);
    fprintf('Paired t-test p = %.4f\n', S1.p);

    S2 = fatigue6C2.scalar.emg.(condLabel);
    fprintf('\n--- Mean EMG | %s | Early vs Late | MID 20-80%% ---\n', condLabel);
    fprintf('Early mean ± SD: %.3f ± %.3f\n', S2.meanEarly, S2.sdEarly);
    fprintf('Late  mean ± SD: %.3f ± %.3f\n', S2.meanLate,  S2.sdLate);
    fprintf('Paired t-test p = %.4f\n', S2.p);
end

fprintf('\n===== Stage 6C2 MID SPM RM-ANOVA SUMMARY =====\n');

print_spm_summary(fatigue6C2.spmRM.torque.effects.contraction);
print_spm_summary(fatigue6C2.spmRM.torque.effects.intensity);
print_spm_summary(fatigue6C2.spmRM.torque.effects.earlyLate);
print_spm_summary(fatigue6C2.spmRM.torque.effects.AxB);
print_spm_summary(fatigue6C2.spmRM.torque.effects.AxC);
print_spm_summary(fatigue6C2.spmRM.torque.effects.BxC);
print_spm_summary(fatigue6C2.spmRM.torque.effects.AxBxC);

print_spm_summary(fatigue6C2.spmRM.emg.effects.contraction);
print_spm_summary(fatigue6C2.spmRM.emg.effects.intensity);
print_spm_summary(fatigue6C2.spmRM.emg.effects.earlyLate);
print_spm_summary(fatigue6C2.spmRM.emg.effects.AxB);
print_spm_summary(fatigue6C2.spmRM.emg.effects.AxC);
print_spm_summary(fatigue6C2.spmRM.emg.effects.BxC);
print_spm_summary(fatigue6C2.spmRM.emg.effects.AxBxC);

fprintf('\n===== Stage 6C2 MID POST HOC PAIRED SPM SUMMARY =====\n');

for c = 1:nCond
    condLabel = char(condNames(c));
    print_spm_summary_paired(fatigue6C2.posthoc.torque.(condLabel));
    print_spm_summary_paired(fatigue6C2.posthoc.emg.(condLabel));
end

%% ---------------------------
%  (11) Figures
%% ---------------------------
if doPlots
    figHandles = struct();

    figHandles.scalarTorque = figure('Color','w','Name','Stage6C2 MID Scalar | Torque Early vs Late');
    ax1 = axes(figHandles.scalarTorque);
    plot_scalar_earlyLate_panel(ax1, Y_torque_early, Y_torque_late, condNames, ...
        'Mean Torque | Early vs Late | MID 20-80%', 'Torque (Nm)');

    figHandles.scalarEMG = figure('Color','w','Name','Stage6C2 MID Scalar | Mean EMG Early vs Late');
    ax2 = axes(figHandles.scalarEMG);
    plot_scalar_earlyLate_panel(ax2, Y_emg_early, Y_emg_late, condNames, ...
        'Mean EMG | Early vs Late | MID 20-80%', 'Mean EMG (%MVC)');

    figHandles.torque_contraction = plot_spm_result_rm(fatigue6C2.spmRM.torque.effects.contraction, fatigue6C2.spmRM.torque.descriptives);
    figHandles.torque_intensity   = plot_spm_result_rm(fatigue6C2.spmRM.torque.effects.intensity,   fatigue6C2.spmRM.torque.descriptives);
    figHandles.torque_earlyLate   = plot_spm_result_rm(fatigue6C2.spmRM.torque.effects.earlyLate,   fatigue6C2.spmRM.torque.descriptives);
    figHandles.torque_AxB         = plot_spm_result_rm(fatigue6C2.spmRM.torque.effects.AxB,         fatigue6C2.spmRM.torque.descriptives);
    figHandles.torque_AxC         = plot_spm_result_rm(fatigue6C2.spmRM.torque.effects.AxC,         fatigue6C2.spmRM.torque.descriptives);
    figHandles.torque_BxC         = plot_spm_result_rm(fatigue6C2.spmRM.torque.effects.BxC,         fatigue6C2.spmRM.torque.descriptives);
    figHandles.torque_AxBxC       = plot_spm_result_rm(fatigue6C2.spmRM.torque.effects.AxBxC,       fatigue6C2.spmRM.torque.descriptives);

    figHandles.emg_contraction = plot_spm_result_rm(fatigue6C2.spmRM.emg.effects.contraction, fatigue6C2.spmRM.emg.descriptives);
    figHandles.emg_intensity   = plot_spm_result_rm(fatigue6C2.spmRM.emg.effects.intensity,   fatigue6C2.spmRM.emg.descriptives);
    figHandles.emg_earlyLate   = plot_spm_result_rm(fatigue6C2.spmRM.emg.effects.earlyLate,   fatigue6C2.spmRM.emg.descriptives);
    figHandles.emg_AxB         = plot_spm_result_rm(fatigue6C2.spmRM.emg.effects.AxB,         fatigue6C2.spmRM.emg.descriptives);
    figHandles.emg_AxC         = plot_spm_result_rm(fatigue6C2.spmRM.emg.effects.AxC,         fatigue6C2.spmRM.emg.descriptives);
    figHandles.emg_BxC         = plot_spm_result_rm(fatigue6C2.spmRM.emg.effects.BxC,         fatigue6C2.spmRM.emg.descriptives);
    figHandles.emg_AxBxC       = plot_spm_result_rm(fatigue6C2.spmRM.emg.effects.AxBxC,       fatigue6C2.spmRM.emg.descriptives);

    for c = 1:nCond
        condLabel = char(condNames(c));
        figHandles.(['posthoc_torque_' condLabel]) = plot_spm_result_paired(fatigue6C2.posthoc.torque.(condLabel));
        figHandles.(['posthoc_emg_' condLabel])    = plot_spm_result_paired(fatigue6C2.posthoc.emg.(condLabel));
    end

    if saveFigures
        outFigDir = uigetdir(pwd, 'Select folder to save Stage6C2 figures');
        if ~isequal(outFigDir,0)
            exportgraphics(figHandles.scalarTorque, fullfile(outFigDir, 'Stage6C2_MID_scalar_torque_earlyLate.png'), 'Resolution', 300);
            exportgraphics(figHandles.scalarEMG,    fullfile(outFigDir, 'Stage6C2_MID_scalar_emg_earlyLate.png'), 'Resolution', 300);

            exportgraphics(figHandles.torque_contraction, fullfile(outFigDir, 'Stage6C2_MID_SPMRM_torque_mainContraction.png'), 'Resolution', 300);
            exportgraphics(figHandles.torque_intensity,   fullfile(outFigDir, 'Stage6C2_MID_SPMRM_torque_mainIntensity.png'),   'Resolution', 300);
            exportgraphics(figHandles.torque_earlyLate,   fullfile(outFigDir, 'Stage6C2_MID_SPMRM_torque_mainEarlyLate.png'),   'Resolution', 300);
            exportgraphics(figHandles.torque_AxB,         fullfile(outFigDir, 'Stage6C2_MID_SPMRM_torque_AxB.png'),             'Resolution', 300);
            exportgraphics(figHandles.torque_AxC,         fullfile(outFigDir, 'Stage6C2_MID_SPMRM_torque_AxC.png'),             'Resolution', 300);
            exportgraphics(figHandles.torque_BxC,         fullfile(outFigDir, 'Stage6C2_MID_SPMRM_torque_BxC.png'),             'Resolution', 300);
            exportgraphics(figHandles.torque_AxBxC,       fullfile(outFigDir, 'Stage6C2_MID_SPMRM_torque_AxBxC.png'),           'Resolution', 300);

            exportgraphics(figHandles.emg_contraction, fullfile(outFigDir, 'Stage6C2_MID_SPMRM_emg_mainContraction.png'), 'Resolution', 300);
            exportgraphics(figHandles.emg_intensity,   fullfile(outFigDir, 'Stage6C2_MID_SPMRM_emg_mainIntensity.png'),   'Resolution', 300);
            exportgraphics(figHandles.emg_earlyLate,   fullfile(outFigDir, 'Stage6C2_MID_SPMRM_emg_mainEarlyLate.png'),   'Resolution', 300);
            exportgraphics(figHandles.emg_AxB,         fullfile(outFigDir, 'Stage6C2_MID_SPMRM_emg_AxB.png'),             'Resolution', 300);
            exportgraphics(figHandles.emg_AxC,         fullfile(outFigDir, 'Stage6C2_MID_SPMRM_emg_AxC.png'),             'Resolution', 300);
            exportgraphics(figHandles.emg_BxC,         fullfile(outFigDir, 'Stage6C2_MID_SPMRM_emg_BxC.png'),             'Resolution', 300);
            exportgraphics(figHandles.emg_AxBxC,       fullfile(outFigDir, 'Stage6C2_MID_SPMRM_emg_AxBxC.png'),           'Resolution', 300);

            for c = 1:nCond
                condLabel = char(condNames(c));
                exportgraphics(figHandles.(['posthoc_torque_' condLabel]), ...
                    fullfile(outFigDir, sprintf('Stage6C2_MID_posthoc_torque_%s.png', condLabel)), 'Resolution', 300);
                exportgraphics(figHandles.(['posthoc_emg_' condLabel]), ...
                    fullfile(outFigDir, sprintf('Stage6C2_MID_posthoc_emg_%s.png', condLabel)), 'Resolution', 300);
            end
        end
    end
end

%% ---------------------------
%  (12) Save results
%% ---------------------------
[outFile, outPath] = uiputfile('stage6C2_earlyLate_mid20to80.mat', 'Save Stage6C2 early-late results as');
if isequal(outFile,0)
    outPath = pwd;
    outFile = 'stage6C2_earlyLate_mid20to80.mat';
end

save(fullfile(outPath, outFile), 'fatigue6C2', '-v7.3');
fprintf('\n✅ Saved Stage 6C2 MID early-late results:\n%s\n', fullfile(outPath, outFile));

%% =========================================================
%  Local functions
%% =========================================================
function H = paired_scalar_earlyLate(yEarly, yLate, labelText)
    [~, p, ~, stats] = ttest(yEarly, yLate);

    H = struct();
    H.label = labelText;
    H.meanEarly = mean(yEarly,'omitnan');
    H.meanLate  = mean(yLate,'omitnan');
    H.sdEarly   = std(yEarly,0,'omitnan');
    H.sdLate    = std(yLate,0,'omitnan');
    H.p = p;
    H.stats = stats;
end

function R = run_spm_anova3rm_earlyLate(Tearly, Tlate, xAxis, alphaLevel, interpClusters, varLabel, units)
% Tearly, Tlate are [nSub x 4 x nT]

    [nSub, nCond, nT] = size(Tearly);
    assert(nCond == 4, 'Expected 4 conditions in dim 2.');
    assert(all(size(Tlate) == [nSub 4 nT]), 'Tlate size mismatch.');

    Y    = zeros(nSub*8, nT);
    A    = zeros(nSub*8, 1);
    B    = zeros(nSub*8, 1);
    C    = zeros(nSub*8, 1);
    SUBJ = zeros(nSub*8, 1);

    row = 0;
    for s = 1:nSub
        for cond = 1:4
            row = row + 1;
            Y(row,:) = squeeze(Tearly(s,cond,:)).';
            [A(row), B(row)] = map_cond_to_AB(cond);
            C(row) = 1;
            SUBJ(row) = s;

            row = row + 1;
            Y(row,:) = squeeze(Tlate(s,cond,:)).';
            [A(row), B(row)] = map_cond_to_AB(cond);
            C(row) = 2;
            SUBJ(row) = s;
        end
    end

    spmList = spm1d.stats.anova3rm(Y, A, B, C, SUBJ);

    spmA   = spmList('A');
    spmB   = spmList('B');
    spmC   = spmList('C');
    spmAB  = spmList('AB');
    spmAC  = spmList('AC');
    spmBC  = spmList('BC');
    spmABC = spmList('ABC');

    spmiA   = local_infer_spm(spmA,   alphaLevel, interpClusters);
    spmiB   = local_infer_spm(spmB,   alphaLevel, interpClusters);
    spmiC   = local_infer_spm(spmC,   alphaLevel, interpClusters);
    spmiAB  = local_infer_spm(spmAB,  alphaLevel, interpClusters);
    spmiAC  = local_infer_spm(spmAC,  alphaLevel, interpClusters);
    spmiBC  = local_infer_spm(spmBC,  alphaLevel, interpClusters);
    spmiABC = local_infer_spm(spmABC, alphaLevel, interpClusters);

    D = struct();
    D.xAxis = xAxis;
    D.varLabel = varLabel;
    D.units = units;
    D.condNames8 = { ...
        'CON_75_Early', 'CON_75_Late', ...
        'ECC_75_Early', 'ECC_75_Late', ...
        'CON_90_Early', 'CON_90_Late', ...
        'ECC_90_Early', 'ECC_90_Late'};

    D.Y = cell(1,8);
    D.mean = cell(1,8);
    D.sd = cell(1,8);

    D.Y{1} = squeeze(Tearly(:,1,:));  D.Y{2} = squeeze(Tlate(:,1,:));
    D.Y{3} = squeeze(Tearly(:,2,:));  D.Y{4} = squeeze(Tlate(:,2,:));
    D.Y{5} = squeeze(Tearly(:,3,:));  D.Y{6} = squeeze(Tlate(:,3,:));
    D.Y{7} = squeeze(Tearly(:,4,:));  D.Y{8} = squeeze(Tlate(:,4,:));

    for i = 1:8
        D.mean{i} = mean(D.Y{i},1,'omitnan');
        D.sd{i}   = std(D.Y{i},0,1,'omitnan');
    end

    R = struct();
    R.varLabel = varLabel;
    R.units = units;
    R.xAxis = xAxis;
    R.descriptives = D;

    R.effects = struct();
    R.effects.contraction = build_effect_struct(spmA,   spmiA,   xAxis, varLabel, units, 'Main effect: Contraction', 'Contraction');
    R.effects.intensity   = build_effect_struct(spmB,   spmiB,   xAxis, varLabel, units, 'Main effect: Intensity', 'Intensity');
    R.effects.earlyLate   = build_effect_struct(spmC,   spmiC,   xAxis, varLabel, units, 'Main effect: Early vs Late', 'EarlyLate');
    R.effects.AxB         = build_effect_struct(spmAB,  spmiAB,  xAxis, varLabel, units, 'Interaction: Contraction × Intensity', 'AxB');
    R.effects.AxC         = build_effect_struct(spmAC,  spmiAC,  xAxis, varLabel, units, 'Interaction: Contraction × Early/Late', 'AxC');
    R.effects.BxC         = build_effect_struct(spmBC,  spmiBC,  xAxis, varLabel, units, 'Interaction: Intensity × Early/Late', 'BxC');
    R.effects.AxBxC       = build_effect_struct(spmABC, spmiABC, xAxis, varLabel, units, 'Interaction: Contraction × Intensity × Early/Late', 'AxBxC');
end

function [a,b] = map_cond_to_AB(cond)
    switch cond
        case 1
            a = 1; b = 1;
        case 2
            a = 2; b = 1;
        case 3
            a = 1; b = 2;
        case 4
            a = 2; b = 2;
        otherwise
            error('Unknown condition index.');
    end
end

function R = run_spm_pair_earlyLate(Yearly, Ylate, xAxis, alphaLevel, twoTailed, interpClusters, varLabel, condLabel, units)
    assert(ismatrix(Yearly) && ismatrix(Ylate), 'Input time series must be 2D matrices.');
    assert(all(size(Yearly) == size(Ylate)), 'Early and Late time series size mismatch.');

    spm = spm1d.stats.ttest_paired(Yearly, Ylate);
    spmi = spm.inference(alphaLevel, 'two_tailed', twoTailed, 'interp', interpClusters);

    R = struct();
    R.varLabel = sprintf('%s | %s | Early vs Late', varLabel, condLabel);
    R.condLabel = condLabel;
    R.units = units;
    R.xAxis = xAxis;

    R.Yearly = Yearly;
    R.Ylate  = Ylate;

    R.meanEarly = mean(Yearly,1,'omitnan');
    R.meanLate  = mean(Ylate,1,'omitnan');
    R.sdEarly   = std(Yearly,0,1,'omitnan');
    R.sdLate    = std(Ylate,0,1,'omitnan');

    R.spm = spm;
    R.spmi = spmi;

    R.z = try_get_field(spmi, 'z');
    R.zstar = try_get_field(spmi, 'zstar');
    R.h0reject = try_get_field(spmi, 'h0reject');
    R.df = try_get_field(spmi, 'df');

    R.clusters = summarize_spm_clusters(spmi, xAxis);
end

function spmi = local_infer_spm(spmObj, alphaLevel, interpClusters)
    try
        spmi = spmObj.inference(alphaLevel, 'interp', interpClusters);
    catch
        spmi = spmObj.inference(alphaLevel);
    end
end

function E = build_effect_struct(spm, spmi, xAxis, varLabel, units, effectLabel, effectName)
    E = struct();
    E.varLabel = varLabel;
    E.units = units;
    E.effectLabel = effectLabel;
    E.effectName = effectName;
    E.xAxis = xAxis;

    E.spm = spm;
    E.spmi = spmi;

    E.z = try_get_field(spmi, 'z');
    E.zstar = try_get_field(spmi, 'zstar');
    E.h0reject = try_get_field(spmi, 'h0reject');
    E.df = try_get_field(spmi, 'df');

    E.clusters = summarize_spm_clusters(spmi, xAxis);
end

function val = try_get_field(S, fieldName)
    val = [];
    try
        if isprop(S, fieldName)
            val = S.(fieldName);
        elseif isfield(S, fieldName)
            val = S.(fieldName);
        end
    catch
        val = [];
    end
end

function C = summarize_spm_clusters(spmi, xAxis)
    C = struct('startPct', {}, 'endPct', {}, 'extentPct', {}, 'P', {});

    try
        if isempty(spmi.clusters)
            return
        end

        nT = numel(xAxis);
        nodeAxis = 0:(nT-1);

        for k = 1:numel(spmi.clusters)
            try
                cl = spmi.clusters{k};
            catch
                cl = spmi.clusters(k);
            end

            startPct = NaN;
            endPct   = NaN;
            extentPct = NaN;
            pval = NaN;

            try
                xy = cl.endpoints;
                if all(isfinite(xy))
                    startPct = interp1(nodeAxis, xAxis, xy(1), 'linear', 'extrap');
                    endPct   = interp1(nodeAxis, xAxis, xy(2), 'linear', 'extrap');
                    extentPct = endPct - startPct;
                end
            catch
                try
                    inds = cl.indices;
                    inds = inds(:);
                    inds = max(min(round(inds), nT-1), 0);
                    startPct = xAxis(min(inds)+1);
                    endPct   = xAxis(max(inds)+1);
                    extentPct = endPct - startPct;
                catch
                end
            end

            try
                pval = cl.P;
            catch
            end

            C(k).startPct = startPct;
            C(k).endPct = endPct;
            C(k).extentPct = extentPct;
            C(k).P = pval;
        end
    catch
    end
end

function print_spm_summary(R)
    fprintf('\n--- %s | %s ---\n', R.varLabel, R.effectLabel);

    if ~isempty(R.zstar)
        fprintf('Threshold F*: %.4f\n', R.zstar);
    else
        fprintf('Threshold F*: unavailable\n');
    end

    if ~isempty(R.h0reject)
        fprintf('Reject H0: %d\n', logical(R.h0reject));
    else
        fprintf('Reject H0: unavailable\n');
    end

    if isempty(R.clusters)
        fprintf('Significant clusters: none\n');
    else
        fprintf('Significant clusters:\n');
        for k = 1:numel(R.clusters)
            fprintf('  Cluster %d: %.1f%% to %.1f%% | extent = %.1f%% | p = %.4f\n', ...
                k, R.clusters(k).startPct, R.clusters(k).endPct, ...
                R.clusters(k).extentPct, R.clusters(k).P);
        end
    end
end

function print_spm_summary_paired(R)
    fprintf('\n--- %s ---\n', R.varLabel);

    if ~isempty(R.zstar)
        fprintf('Threshold z*: %.4f\n', R.zstar);
    else
        fprintf('Threshold z*: unavailable\n');
    end

    if ~isempty(R.h0reject)
        fprintf('Reject H0: %d\n', logical(R.h0reject));
    else
        fprintf('Reject H0: unavailable\n');
    end

    if isempty(R.clusters)
        fprintf('Significant clusters: none\n');
    else
        fprintf('Significant clusters:\n');
        for k = 1:numel(R.clusters)
            fprintf('  Cluster %d: %.1f%% to %.1f%% | extent = %.1f%% | p = %.4f\n', ...
                k, R.clusters(k).startPct, R.clusters(k).endPct, ...
                R.clusters(k).extentPct, R.clusters(k).P);
        end
    end
end

function fig = plot_spm_result_rm(R, D)
    fig = figure('Color','w', 'Name', sprintf('SPM | %s | %s', R.varLabel, R.effectLabel));
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl,1); hold(ax1,'on');

    cols = [...
        0.00 0.4470 0.7410;
        0.3010 0.7450 0.9330;
        0.4660 0.6740 0.1880;
        0.7290 0.8730 0.5410;
        0.8500 0.3250 0.0980;
        0.9290 0.6940 0.1250;
        0.4940 0.1840 0.5560;
        0.6350 0.0780 0.1840];

    h = gobjects(1,8);
    for c = 1:8
        fill(ax1, [D.xAxis fliplr(D.xAxis)], ...
            [D.mean{c}-D.sd{c} fliplr(D.mean{c}+D.sd{c})], ...
            cols(c,:), 'FaceAlpha',0.08, 'EdgeColor','none', 'HandleVisibility','off');

        h(c) = plot(ax1, D.xAxis, D.mean{c}, 'Color', cols(c,:), 'LineWidth',1.8);
    end

    ylabel(ax1, sprintf('%s (%s)', R.varLabel, R.units), 'FontWeight','bold');
    title(ax1, sprintf('%s | %s', R.varLabel, R.effectLabel), ...
        'Interpreter','none', 'FontWeight','bold');

    legend(ax1, h, strrep(D.condNames8,'_',' '), 'Box','off', 'Location','eastoutside');
    box(ax1,'off');
    set(ax1,'TickDir','out', 'LineWidth',1.2, 'FontWeight','bold');
    xlim(ax1, [D.xAxis(1) D.xAxis(end)]);

    ax2 = nexttile(tl,2); hold(ax2,'on');

    if ~isempty(R.z)
        hSPM = plot(ax2, R.xAxis, R.z, 'k', 'LineWidth',1.8, 'DisplayName','SPM{F}');
    else
        hSPM = gobjects(0);
    end

    if ~isempty(R.zstar)
        yline(ax2, R.zstar, '--r', 'LineWidth',1.2, 'HandleVisibility','off');
    end

    if ~isempty(R.clusters)
        yl = ylim(ax2);
        for k = 1:numel(R.clusters)
            x1 = R.clusters(k).startPct;
            x2 = R.clusters(k).endPct;
            patch(ax2, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.75 0.75 0.75], 'FaceAlpha',0.25, ...
                'EdgeColor','none', 'HandleVisibility','off');
        end
        ylim(ax2, yl);
        uistack(findobj(ax2,'Type','line'),'top');
    end

    xlabel(ax2, 'Normalized movement (%)', 'FontWeight','bold');
    ylabel(ax2, 'SPM{F}', 'FontWeight','bold');
    title(ax2, sprintf('SPM RM ANOVA | alpha = %.3f', try_get_field(R.spmi,'alpha')), ...
        'Interpreter','none', 'FontWeight','bold');

    if ~isempty(hSPM)
        legend(ax2, hSPM, 'SPM{F}', 'Box','off', 'Location','best');
    end

    box(ax2,'off');
    set(ax2,'TickDir','out', 'LineWidth',1.2, 'FontWeight','bold');
    xlim(ax2, [R.xAxis(1) R.xAxis(end)]);
end

function fig = plot_spm_result_paired(R)
    fig = figure('Color','w', 'Name', sprintf('SPM | %s', R.varLabel));
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl,1); hold(ax1,'on');

    col1 = [0.00 0.4470 0.7410];
    col2 = [0.8500 0.3250 0.0980];

    fill(ax1, [R.xAxis fliplr(R.xAxis)], [R.meanEarly-R.sdEarly fliplr(R.meanEarly+R.sdEarly)], ...
        col1, 'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
    fill(ax1, [R.xAxis fliplr(R.xAxis)], [R.meanLate-R.sdLate fliplr(R.meanLate+R.sdLate)], ...
        col2, 'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');

    h1 = plot(ax1, R.xAxis, R.meanEarly, 'Color', col1, 'LineWidth',2.2);
    h2 = plot(ax1, R.xAxis, R.meanLate,  'Color', col2, 'LineWidth',2.2);

    ylabel(ax1, sprintf('%s (%s)', erase(R.varLabel, [' | ' R.condLabel ' | Early vs Late']), R.units), ...
        'FontWeight','bold');
    title(ax1, strrep(R.varLabel,'_',' '), 'Interpreter','none', 'FontWeight','bold');
    legend(ax1, [h1 h2], {'Early','Late'}, 'Box','off', 'Location','best');
    box(ax1,'off');
    set(ax1,'TickDir','out', 'LineWidth',1.2, 'FontWeight','bold');
    xlim(ax1, [R.xAxis(1) R.xAxis(end)]);

    ax2 = nexttile(tl,2); hold(ax2,'on');

    if ~isempty(R.z)
        hSPM = plot(ax2, R.xAxis, R.z, 'k', 'LineWidth',1.8, 'DisplayName','SPM{t}');
    else
        hSPM = gobjects(0);
    end

    if ~isempty(R.zstar)
        yline(ax2, R.zstar, '--r', 'LineWidth',1.2, 'HandleVisibility','off');
        if ~isempty(R.spmi) && try_get_field(R.spmi, 'two_tailed')
            yline(ax2, -R.zstar, '--r', 'LineWidth',1.2, 'HandleVisibility','off');
        end
    end

    if ~isempty(R.clusters)
        yl = ylim(ax2);
        for k = 1:numel(R.clusters)
            x1 = R.clusters(k).startPct;
            x2 = R.clusters(k).endPct;
            patch(ax2, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.75 0.75 0.75], 'FaceAlpha',0.25, ...
                'EdgeColor','none', 'HandleVisibility','off');
        end
        ylim(ax2, yl);
        uistack(findobj(ax2,'Type','line'),'top');
    end

    xlabel(ax2, 'Normalized movement (%)', 'FontWeight','bold');
    ylabel(ax2, 'SPM{t}', 'FontWeight','bold');
    title(ax2, sprintf('SPM paired t-test | alpha = %.3f', try_get_field(R.spmi,'alpha')), ...
        'Interpreter','none', 'FontWeight','bold');

    if ~isempty(hSPM)
        legend(ax2, hSPM, 'SPM{t}', 'Box','off', 'Location','best');
    end

    box(ax2,'off');
    set(ax2,'TickDir','out', 'LineWidth',1.2, 'FontWeight','bold');
    xlim(ax2, [R.xAxis(1) R.xAxis(end)]);
end

function plot_scalar_earlyLate_panel(ax, Yearly, Ylate, condNames, ttl, ylab)
    hold(ax,'on');

    xPos = [1 2 4 5 7 8 10 11];
    labels = {'E','L','E','L','E','L','E','L'};

    nSub = size(Yearly,1);
    subjColors = lines(nSub);

    Y = nan(nSub, 8);
    Y(:,1) = Yearly(:,1); Y(:,2) = Ylate(:,1);
    Y(:,3) = Yearly(:,2); Y(:,4) = Ylate(:,2);
    Y(:,5) = Yearly(:,3); Y(:,6) = Ylate(:,3);
    Y(:,7) = Yearly(:,4); Y(:,8) = Ylate(:,4);

    for i = 1:nSub
        yi = Y(i,:);
        col = subjColors(i,:);

        for k = 1:2:8
            plot(ax, xPos(k:k+1), yi(k:k+1), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
            scatter(ax, xPos(k),   yi(k),   28, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
            scatter(ax, xPos(k+1), yi(k+1), 28, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
        end
    end

    m = mean(Y,1,'omitnan');
    s = std(Y,0,1,'omitnan');

    for k = 1:2:8
        plot(ax, xPos(k:k+1), m(k:k+1), '--', 'Color','k', 'LineWidth',1.4);
    end

    for c = 1:8
        plot(ax, [xPos(c) xPos(c)], [m(c)-s(c), m(c)+s(c)], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
        capW = 0.10;
        plot(ax, [xPos(c)-capW xPos(c)+capW], [m(c)-s(c) m(c)-s(c)], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
        plot(ax, [xPos(c)-capW xPos(c)+capW], [m(c)+s(c) m(c)+s(c)], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
    end

    scatter(ax, xPos, m, 42, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

    yTop = max(Y(:),[],'omitnan');
    yBot = min(Y(:),[],'omitnan');
    yRange = yTop - yBot;
    if yRange == 0, yRange = 1; end
    ylim(ax, [yBot - 0.12*yRange, yTop + 0.18*yRange]);

    pairCenters = [1.5 4.5 7.5 10.5];
    for i = 1:numel(pairCenters)
        text(ax, pairCenters(i), yTop + 0.10*yRange, strrep(condNames(i),'_',' '), ...
            'HorizontalAlignment','center', 'FontWeight','bold');
    end

    set(ax, 'XTick', xPos, 'XTickLabel', labels, ...
        'FontSize',10, 'FontWeight','bold', 'LineWidth',1.4, 'TickDir','out');

    ylabel(ax, ylab, 'FontWeight','bold');
    title(ax, ttl, 'Interpreter','none', 'FontWeight','bold');
    box(ax,'off');

    if yBot < 0 && yTop > 0
        yline(ax, 0, ':k', 'LineWidth', 1.0);
    end
end