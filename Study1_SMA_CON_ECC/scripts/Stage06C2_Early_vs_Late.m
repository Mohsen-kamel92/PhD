%% ===== Stage 06C2: Early against late repetitions within the analysis window =====
%
% PURPOSE
%   Test whether the outcomes changed between the first and the last
%   repetitions of a set, and whether that change differed between contraction
%   modes or intensities.
%
% SUPPORTS
%   The within-trial modulation reported in the Results, and the two
%   supplementary figures covering repetition effects.
%
% INPUT
%   <ROOT>/groupData_stage6A.mat                      (Stage 06A)
%
% OUTPUT
%   <ROOT>/stage6C2_earlyLate_mid20to80.mat
%   <ROOT>/Figures_Stage6C2/*.png
%   <ROOT>/Figures_Stage6C2/Stage6C2_SFig_emg_repetition.{png,pdf}
%   <ROOT>/Figures_Stage6C2/Stage6C2_SFig_centroid_repetition.{png,pdf}
%   <ROOT>/stage6C2_clusters.csv
%   <ROOT>/stage6C2_scalar.csv
%
% DESIGN
%   2 x 2 x 2 repeated measures SPM ANOVA on mean EMG, centroid X and centroid
%   Y. Factor A is contraction type (CON, ECC), B is intensity (75, 90), and C
%   is repetition bin (early, late). Paired SPM t-tests compare early against
%   late within each of the four conditions.
%
%   Torque is tested only as a scalar. The claim under test is that torque
%   stayed stable across repetitions while activity changed, so torque is a
%   manipulation check rather than an outcome.
%
% ON THE EARLY TO LATE CONTRAST
%   Sets at 75% MVC hold 8 repetitions, so early (1 and 2) and late (7 and 8)
%   are separated by four intervening contractions. Sets at 90% hold 4, so
%   early (1 and 2) and late (3 and 4) are adjacent. The contrast therefore
%   spans a different amount of preceding work at the two intensities, which
%   matters when interpreting any difference between them in the early to late
%   change.
%
% WHY ECC IS REVERSED
%   Same convention as Stage 6C1. The two modes traverse the range of motion in
%   opposite directions, so ECC is reversed in time to align them by joint angle
%   rather than by elapsed time. The scalar means are unaffected by the
%   reversal.
%
% WHICH CONDITIONS EACH PANEL SHOWS
%   The upper panel shows only the grouping the test below it compares, so a
%   main effect of repetition bin shows early and late collapsed across the
%   other factors rather than all eight conditions.
%
% ON THE SHADED BANDS
%   Centroid position differs by up to 11 mm between participants, because the
%   grid was positioned relative to each participant's own muscle-tendon
%   junction and their anatomy differed. For centroid outcomes the bands are
%   therefore computed after removing each participant's own mean across the
%   groups shown, so they reflect within-subject variability. EMG keeps
%   conventional between-subject SDs.
%
% ON APPROXIMATE RESIDUALS
%   With one observation per participant per condition, spm1d computes
%   residuals approximately and prints a warning to that effect. This is
%   expected for this design and is noted in the manuscript.
%
% SUPPLEMENTARY COMPOSITES
%   Section (11) assembles two figures: EMG across repetitions in four panels,
%   and centroid across repetitions in two. Styling matches Stages 06C1 and 06C4.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped.
%
% DEPENDENCIES
%   spm1d for MATLAB (https://spm1d.org)
%   distinguishable_colors.m  participant colours  (MATLAB File Exchange)

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'midWin',             [20 80]);   % analysis window, % of CV
CFG = set_default(CFG, 'alphaLevel',         0.05);
CFG = set_default(CFG, 'twoTailed',          true);
CFG = set_default(CFG, 'interpClusters',     true);
CFG = set_default(CFG, 'flipECCforAnalysis', true);
CFG = set_default(CFG, 'doPlots',            true);
CFG = set_default(CFG, 'saveFigures',        true);
CFG = set_default(CFG, 'showTitles',         true);

% Composite supplementary figures, section (11). The SPM limits were chosen
% from the ranges printed at the end of that section on an earlier run.
CFG = set_default(CFG, 'comp_emgYLim',  [20 100]);   yStepEMG = 20;
CFG = set_default(CFG, 'comp_cenYLim',  [-2 2]);     yStepCen = 1;
CFG = set_default(CFG, 'comp_spmYLimF', [-8 100]);   % F panels, EMG composite
CFG = set_default(CFG, 'comp_spmStepF', 20);
CFG = set_default(CFG, 'comp_spmYLimT', [-10 6]);    % t panels, EMG composite
CFG = set_default(CFG, 'comp_spmStepT', 5);
CFG = set_default(CFG, 'comp_spmYLimC', [-1.5 20]);  % F panels, centroid composite
CFG = set_default(CFG, 'comp_spmStepC', 5);

% Styling, matched to Stages 6B, 6C1 and 6C4
STY = struct('FontSize',12, 'AxLineWidth',2.5, 'TraceWidth',2.2, ...
             'SPMWidth',1.8, 'FaceAlpha',0.15, 'TickLen',0.008);

midWin             = CFG.midWin;
alphaLevel         = CFG.alphaLevel;
twoTailed          = CFG.twoTailed;
interpClusters     = CFG.interpClusters;
flipECCforAnalysis = CFG.flipECCforAnalysis;

%% ---------------------------
%  (2) Data root and group data
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 6C2: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

spmFun = which('spm1d.stats.anova3rm');
assert(~isempty(spmFun), ...
    ['Stage 6C2: spm1d not found on the MATLAB path. ' ...
     'Add it with addpath(genpath(''path_to_spm1d'')).']);
assert(~isempty(which('distinguishable_colors')), ...
    'Stage 6C2: distinguishable_colors not found. Add the /functions folder.');

groupFile = fullfile(ROOT,'groupData_stage6A.mat');
assert(isfile(groupFile), ...
    'Stage 6C2: %s not found. Run Stage 6A first.', groupFile);

L = load(groupFile);
assert(isfield(L,'groupData'), 'File does not contain "groupData".');
groupData = L.groupData;

condNames     = string(groupData.conditions(:));
expectedOrder = ["CON_75"; "ECC_75"; "CON_90"; "ECC_90"];
assert(all(condNames == expectedOrder), ...
    'Condition order mismatch. Expected [CON_75, ECC_75, CON_90, ECC_90].');

for b = {'early','late'}
    assert(isfield(groupData.timeseries, b{1}), 'Missing timeseries.%s', b{1});
end

nSub       = groupData.nSub;
nCond      = groupData.nCond;
xAxis_full = groupData.xAxis(:).';
midMask    = xAxis_full >= midWin(1) & xAxis_full <= midWin(2);
assert(any(midMask), 'No samples inside the requested window.');

xAxis = xAxis_full(midMask);
nT    = numel(xAxis);

figDir = fullfile(ROOT,'Figures_Stage6C2');
if CFG.saveFigures && ~exist(figDir,'dir'), mkdir(figDir); end

fprintf('\n===== Stage 6C2: early against late =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('spm1d found at %s\n', spmFun);
fprintf('Participants %d | alpha %.3f | window %d-%d %%CV | ECC reversed %d\n', ...
    nSub, alphaLevel, midWin(1), midWin(2), flipECCforAnalysis);

%% ---------------------------
%  (3) Extract and align
%% ---------------------------
getTS = @(bin, fld) groupData.timeseries.(bin).(fld)(:,:,midMask);

T = struct();
T.torque.early = getTS('early','torque_101');
T.torque.late  = getTS('late', 'torque_101');
T.emg.early    = getTS('early','emgMeanNorm_101');
T.emg.late     = getTS('late', 'emgMeanNorm_101');
T.cx.early     = getTS('early','centroidXc_101');
T.cx.late      = getTS('late', 'centroidXc_101');
T.cy.early     = getTS('early','centroidYc_101');
T.cy.late      = getTS('late', 'centroidYc_101');

fn = fieldnames(T);
for i = 1:numel(fn)
    assert(isequal(size(T.(fn{i}).early), [nSub nCond nT]), 'Size mismatch: %s early.', fn{i});
    assert(isequal(size(T.(fn{i}).late),  [nSub nCond nT]), 'Size mismatch: %s late.',  fn{i});
end

if flipECCforAnalysis
    for i = 1:numel(fn)
        for c = [2 4]
            T.(fn{i}).early(:,c,:) = flip(T.(fn{i}).early(:,c,:), 3);
            T.(fn{i}).late(:,c,:)  = flip(T.(fn{i}).late(:,c,:),  3);
        end
    end
end

% Window means, for the scalar tests
Y = struct();
for i = 1:numel(fn)
    Y.(fn{i}).early = mean(T.(fn{i}).early, 3, 'omitnan');
    Y.(fn{i}).late  = mean(T.(fn{i}).late,  3, 'omitnan');
end

%% ---------------------------
%  (4) Output structure
%% ---------------------------
fatigue6C2 = struct();
fatigue6C2.meta = struct( ...
    'createdOn', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'script', mfilename, 'sourceFile', groupFile, ...
    'nSub', nSub, 'nCond', nCond, 'nTime', nT, ...
    'xAxis', xAxis, 'xAxis_label', '%CV', 'alpha', alphaLevel, ...
    'conditions', condNames, 'flipECCforAnalysis', flipECCforAnalysis, ...
    'analysisWindow', midWin, 'design', '2x2x2 repeated-measures ANOVA', ...
    'spmOutcomes', {{'emg','cx','cy'}}, ...
    'factorA', 'Contraction', 'factorALevels', {{'CON','ECC'}}, ...
    'factorB', 'Intensity',   'factorBLevels', {{'75','90'}}, ...
    'factorC', 'Repetition',  'factorCLevels', {{'Early','Late'}}, ...
    'torque_note', "torque is tested only as a scalar, as a stability check", ...
    'binning_note', ["early and late are the first two and last two valid " ...
                     "repetitions; at 75% these are separated by four " ...
                     "contractions, at 90% they are adjacent"], ...
    'band_note', ["centroid bands are within-subject: each participant's own " ...
                  "mean across the groups shown is removed first"]);

%% ---------------------------
%  (5) Scalar tests, Holm corrected
%% ---------------------------
fatigue6C2.scalar = struct();
praw = []; ref = {};

for i = 1:numel(fn)
    nm = fn{i};
    for c = 1:nCond
        cl = char(condNames(c));
        H = paired_scalar(Y.(nm).early(:,c), Y.(nm).late(:,c), ...
            sprintf('%s | %s | early vs late', nm, cl));
        fatigue6C2.scalar.(nm).(cl) = H;

        praw(end+1) = H.p;      %#ok<SAGROW>
        ref{end+1}  = {nm, cl}; %#ok<SAGROW>
    end
end

% All sixteen tests enter one Holm family, including the torque stability
% check, so no comparison is corrected more leniently than the others.
pHolm = holm(praw);
for i = 1:numel(ref)
    fatigue6C2.scalar.(ref{i}{1}).(ref{i}{2}).p_holm = pHolm(i);
end

fprintf('\n===== SCALAR: EARLY vs LATE (Holm corrected across %d tests) =====\n', numel(praw));
fprintf('%-32s %15s %15s %9s %8s %8s\n','test','early','late','change','p','p_holm');
scalarRows = {};
for i = 1:numel(ref)
    H = fatigue6C2.scalar.(ref{i}{1}).(ref{i}{2});
    fprintf('%-32s %7.2f +/-%5.2f %7.2f +/-%5.2f %+9.2f %8.4f %8.4f  dz = %+.2f\n', ...
        H.label, H.meanEarly, H.sdEarly, H.meanLate, H.sdLate, ...
        H.change, H.p, H.p_holm, H.dz);
    scalarRows(end+1,:) = {string(ref{i}{1}), string(ref{i}{2}), ...
        H.meanEarly, H.sdEarly, H.meanLate, H.sdLate, H.change, ...
        H.stats.tstat, H.stats.df, H.p, H.p_holm, H.dz}; %#ok<SAGROW>
end

writetable(cell2table(scalarRows, 'VariableNames', ...
    {'outcome','condition','early','early_sd','late','late_sd','change','t','df','p','p_holm','dz'}), ...
    fullfile(ROOT,'stage6C2_scalar.csv'));

%% ---------------------------
%  (6) SPM three-way RM-ANOVA
%% ---------------------------
spmOutcomes = {'emg','cx','cy'};
spmLabels   = {'RMS EMG','Centroid X','Centroid Y'};
spmUnits    = {'% MVA','mm','mm'};
isCentroid  = [false true true];

fatigue6C2.spmRM = struct();
for i = 1:numel(spmOutcomes)
    nm = spmOutcomes{i};
    fatigue6C2.spmRM.(nm) = run_anova3rm(T.(nm).early, T.(nm).late, xAxis, ...
        alphaLevel, interpClusters, spmLabels{i}, spmUnits{i}, isCentroid(i));
end

fprintf('\n===== SPM 3-WAY RM-ANOVA =====\n');
effNames = {'contraction','intensity','repetition','AxB','AxC','BxC','AxBxC'};
for i = 1:numel(spmOutcomes)
    R = fatigue6C2.spmRM.(spmOutcomes{i});
    for e = 1:numel(effNames)
        print_summary(R.effects.(effNames{e}));
    end
end

%% ---------------------------
%  (7) Paired SPM t-tests, early against late within condition
%% ---------------------------
fatigue6C2.posthoc = struct();
fprintf('\n===== PAIRED SPM t-TESTS (early vs late within condition) =====\n');

for i = 1:numel(spmOutcomes)
    nm = spmOutcomes{i};
    for c = 1:nCond
        cl = char(condNames(c));
        spmTi = infer(spm1d.stats.ttest_paired( ...
            squeeze(T.(nm).early(:,c,:)), squeeze(T.(nm).late(:,c,:))), ...
            alphaLevel, interpClusters, twoTailed);

        E = struct('varLabel', spmLabels{i}, 'units', spmUnits{i}, ...
            'effectLabel', sprintf('early vs late, %s', strrep(cl,'_',' ')), ...
            'effectName', "posthoc_"+cl, 'xAxis', xAxis, 'spmi', spmTi, ...
            'z', getf(spmTi,'z'), 'zstar', getf(spmTi,'zstar'), ...
            'h0reject', getf(spmTi,'h0reject'), 'df', getf(spmTi,'df'), ...
            'clusters', clusters_of(spmTi, xAxis), 'statLabel', 'SPM{t}', ...
            'groups', {{2*c-1, 2*c}}, 'groupNames', {{'Early','Late'}});

        fatigue6C2.posthoc.(nm).(cl) = E;
        print_summary(E);
    end
end

%% ---------------------------
%  (8) Cluster table
%% ---------------------------
rows = {};
for i = 1:numel(spmOutcomes)
    nm = spmOutcomes{i};
    R  = fatigue6C2.spmRM.(nm);
    eff = cellfun(@(e) R.effects.(e), effNames, 'UniformOutput', false);
    for c = 1:nCond
        eff{end+1} = fatigue6C2.posthoc.(nm).(char(condNames(c))); %#ok<SAGROW>
    end

    for e = 1:numel(eff)
        E = eff{e};
        if isempty(E.clusters)
            rows(end+1,:) = {string(E.varLabel), string(E.effectLabel), 0, ...
                             NaN, NaN, NaN, NaN}; %#ok<SAGROW>
        else
            for k = 1:numel(E.clusters)
                rows(end+1,:) = {string(E.varLabel), string(E.effectLabel), k, ...
                    E.clusters(k).startPct, E.clusters(k).endPct, ...
                    E.clusters(k).extentPct, E.clusters(k).P}; %#ok<SAGROW>
            end
        end
    end
end

writetable(cell2table(rows, 'VariableNames', ...
    {'outcome','effect','cluster','startPct','endPct','extentPct','p'}), ...
    fullfile(ROOT,'stage6C2_clusters.csv'));

%% ---------------------------
%  (9) Figures, one per effect
%% ---------------------------
% Every effect is drawn and saved. Which of them reach the supplement is decided
% in section (11); keeping the full set here means a result can be produced
% later without rerunning the analysis.

if CFG.doPlots
    figs = struct();

    for i = 1:numel(fn)
        nm = fn{i};
        figs.(['scalar_' nm]) = scalar_figure(Y.(nm).early, Y.(nm).late, ...
            condNames, nm, CFG.showTitles, STY);
    end

    for i = 1:numel(spmOutcomes)
        nm = spmOutcomes{i};
        R  = fatigue6C2.spmRM.(nm);
        for e = 1:numel(effNames)
            figs.([nm '_' effNames{e}]) = ...
                plot_result(R.effects.(effNames{e}), R.descriptives, CFG.showTitles, STY);
        end
        for c = 1:nCond
            cl = char(condNames(c));
            figs.([nm '_posthoc_' cl]) = ...
                plot_result(fatigue6C2.posthoc.(nm).(cl), R.descriptives, CFG.showTitles, STY);
        end
    end

    if CFG.saveFigures
        names = fieldnames(figs);
        for i = 1:numel(names)
            exportgraphics(figs.(names{i}), ...
                fullfile(figDir, sprintf('Stage6C2_%s.png', names{i})), 'Resolution', 300);
        end
        fprintf('\nFigures saved to %s\n', figDir);
    end
end

%% ---------------------------
%  (10) Composite supplementary figures
%% ---------------------------
%  EMG across repetitions
%    A  main effect of repetition bin
%    B  contraction by repetition interaction
%    C  early against late within CON 75
%    D  early against late within CON 90
%
%  ECC 75 produced no clusters and ECC 90 produced two brief ones, both stated
%  in the text, so neither needs a panel.
%
%  Centroid across repetitions
%    A  centroid X, main effect of repetition bin
%    B  centroid Y, main effect of repetition bin
%
%  The SPM ranges printed at the end of this section are how the shared
%  y-limits in section (1) were chosen.

if CFG.doPlots
    Remg = fatigue6C2.spmRM.emg;
    Rcx  = fatigue6C2.spmRM.cx;
    Rcy  = fatigue6C2.spmRM.cy;

    % ---- EMG composite ----
    figA = figure('Color','w','Name','Supplementary: EMG across repetitions', ...
                  'Position',[100 100 1250 950]);
    outerA = tiledlayout(figA,2,2,'TileSpacing','compact','Padding','compact');

    %          effect                          descriptives       letter  yLim              spmYLim            spmStep            xLab   legend xTick  yTick
    panelsA = { Remg.effects.repetition,       Remg.descriptives, 'A', CFG.comp_emgYLim, CFG.comp_spmYLimF, CFG.comp_spmStepF, false, true,  false, true;
                Remg.effects.AxC,              Remg.descriptives, 'B', CFG.comp_emgYLim, CFG.comp_spmYLimF, CFG.comp_spmStepF, false, true,  false, false;
                fatigue6C2.posthoc.emg.CON_75, Remg.descriptives, 'C', CFG.comp_emgYLim, CFG.comp_spmYLimT, CFG.comp_spmStepT, true,  true,  true,  true;
                fatigue6C2.posthoc.emg.CON_90, Remg.descriptives, 'D', CFG.comp_emgYLim, CFG.comp_spmYLimT, CFG.comp_spmStepT, true,  true,  true,  false };

    for k = 1:size(panelsA,1)
        innerA = tiledlayout(outerA,2,1,'TileSpacing','compact','Padding','tight');
        innerA.Layout.Tile = k;

        optsA = struct('panelLetter',  panelsA{k,3}, ...
                       'yLim',         panelsA{k,4}, ...
                       'yTickStep',    yStepEMG, ...
                       'spmYLim',      panelsA{k,5}, ...
                       'spmYTickStep', panelsA{k,6}, ...
                       'showXLabel',   panelsA{k,7}, ...
                       'showLegend',   panelsA{k,8}, ...
                       'showXTicks',   panelsA{k,9}, ...
                       'showYTicks',   panelsA{k,10});

        draw_spm_pair(innerA, panelsA{k,1}, panelsA{k,2}, false, STY, optsA);
    end

    % ---- centroid composite ----
    figB = figure('Color','w','Name','Supplementary: centroid across repetitions', ...
                  'Position',[100 100 900 950]);
    outerB = tiledlayout(figB,2,1,'TileSpacing','compact','Padding','compact');

    panelsB = { Rcx.effects.repetition, Rcx.descriptives, 'A', false, false;
                Rcy.effects.repetition, Rcy.descriptives, 'B', true,  true };

    for k = 1:size(panelsB,1)
        innerB = tiledlayout(outerB,2,1,'TileSpacing','compact','Padding','tight');
        innerB.Layout.Tile = k;

        optsB = struct('panelLetter',  panelsB{k,3}, ...
                       'yLim',         CFG.comp_cenYLim, ...
                       'yTickStep',    yStepCen, ...
                       'spmYLim',      CFG.comp_spmYLimC, ...
                       'spmYTickStep', CFG.comp_spmStepC, ...
                       'showXLabel',   panelsB{k,4}, ...
                       'showLegend',   true, ...
                       'showXTicks',   panelsB{k,5}, ...
                       'showYTicks',   true);

        draw_spm_pair(innerB, panelsB{k,1}, panelsB{k,2}, false, STY, optsB);
    end

    fprintf('\n===== SPM RANGES FOR THE COMPOSITE PANELS =====\n');
    report_range('A  EMG repetition', Remg.effects.repetition);
    report_range('B  EMG AxC',        Remg.effects.AxC);
    report_range('C  EMG CON 75',     fatigue6C2.posthoc.emg.CON_75);
    report_range('D  EMG CON 90',     fatigue6C2.posthoc.emg.CON_90);
    report_range('A  centroid X rep', Rcx.effects.repetition);
    report_range('B  centroid Y rep', Rcy.effects.repetition);

    if CFG.saveFigures
        exportgraphics(figA, fullfile(figDir,'Stage6C2_SFig_emg_repetition.png'), ...
            'Resolution', 300);
        exportgraphics(figA, fullfile(figDir,'Stage6C2_SFig_emg_repetition.pdf'), ...
            'ContentType','vector');
        exportgraphics(figB, fullfile(figDir,'Stage6C2_SFig_centroid_repetition.png'), ...
            'Resolution', 300);
        exportgraphics(figB, fullfile(figDir,'Stage6C2_SFig_centroid_repetition.pdf'), ...
            'ContentType','vector');
        fprintf('\nComposite figures saved to %s\n', figDir);
    end
end

%% ---------------------------
%  (11) Save
%% ---------------------------
outName = fullfile(ROOT, 'stage6C2_earlyLate_mid20to80.mat');
save(outName, 'fatigue6C2', '-v7.3');
fprintf('Saved results : %s\n', outName);
fprintf('Saved scalars : %s\n', fullfile(ROOT,'stage6C2_scalar.csv'));
fprintf('Saved clusters: %s\n', fullfile(ROOT,'stage6C2_clusters.csv'));

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

    assert(hasParticipant, ...
        ['The selected folder does not look like the expected data root.\n' ...
         'Expected <ROOT>/<participant>/<CON_75|CON_90|ECC_75|ECC_90>/<Set_N>/\n' ...
         'Selected: %s'], ROOT);
end

function pAdj = holm(p)
    [ps, idx] = sort(p(:));
    m = numel(ps);
    adj = min(cummax(ps .* (m - (1:m)' + 1)), 1);
    pAdj = nan(size(ps));
    pAdj(idx) = adj;
    pAdj = reshape(pAdj, size(p));
end

function H = paired_scalar(yE, yL, labelText)
% Change is late minus early, so a positive value means the outcome rose across
% the set.
    [~, pv, ci, st] = ttest(yE, yL);
    d = yL - yE;
    H = struct('label', labelText, ...
        'meanEarly', mean(yE,'omitnan'), 'meanLate', mean(yL,'omitnan'), ...
        'sdEarly', std(yE,0,'omitnan'),  'sdLate', std(yL,0,'omitnan'), ...
        'change', mean(d,'omitnan'), 'sdChange', std(d,0,'omitnan'), ...
        'dz', mean(d,'omitnan')/std(d,0,'omitnan'), ...
        'ci', ci(:).', 'p', pv, 'p_holm', NaN, 'stats', st);
end

function spmi = infer(spmObj, alphaLevel, interpClusters, twoTailed)
% Older spm1d versions do not accept the interp or two_tailed options, so the
% call falls back to the plain form.
    if nargin < 4, twoTailed = []; end
    try
        if isempty(twoTailed)
            spmi = spmObj.inference(alphaLevel, 'interp', interpClusters);
        else
            spmi = spmObj.inference(alphaLevel, 'two_tailed', twoTailed, ...
                'interp', interpClusters);
        end
    catch
        spmi = spmObj.inference(alphaLevel);
    end
end

function R = run_anova3rm(Tearly, Tlate, xAxis, alphaLevel, interpClusters, ...
                          varLabel, units, isCentroid)
% Reshapes the two [nSub x 4 x nT] arrays into the long form spm1d expects,
% runs the three-way repeated measures ANOVA, and packages each effect with the
% conditions its panel should show.
    [nSub, nCond, nT] = size(Tearly);
    assert(nCond == 4, 'Expected 4 conditions.');

    Y = zeros(nSub*8, nT);
    A = zeros(nSub*8,1); B = zeros(nSub*8,1);
    C = zeros(nSub*8,1); S = zeros(nSub*8,1);

    row = 0;
    for s = 1:nSub
        for c = 1:4
            [a,b] = ab_of(c);
            row = row + 1;
            Y(row,:) = squeeze(Tearly(s,c,:)).';
            A(row)=a; B(row)=b; C(row)=1; S(row)=s;

            row = row + 1;
            Y(row,:) = squeeze(Tlate(s,c,:)).';
            A(row)=a; B(row)=b; C(row)=2; S(row)=s;
        end
    end

    sl = spm1d.stats.anova3rm(Y, A, B, C, S);

    % Descriptives, indexed 1 to 8 as
    % 1 CON75E  2 CON75L  3 ECC75E  4 ECC75L
    % 5 CON90E  6 CON90L  7 ECC90E  8 ECC90L
    D = struct('xAxis', xAxis, 'varLabel', varLabel, 'units', units, ...
        'isCentroid', isCentroid, ...
        'condNames', {{'CON 75 E','CON 75 L','ECC 75 E','ECC 75 L', ...
                       'CON 90 E','CON 90 L','ECC 90 E','ECC 90 L'}}, ...
        'Y', {cell(1,8)});
    for c = 1:4
        D.Y{2*c-1} = squeeze(Tearly(:,c,:));
        D.Y{2*c}   = squeeze(Tlate(:,c,:));
    end

    R = struct('varLabel', varLabel, 'units', units, 'xAxis', xAxis, 'descriptives', D);

    mk = @(key, lbl, nm, grp, gnm) build(infer(sl(key), alphaLevel, interpClusters), ...
        xAxis, varLabel, units, lbl, nm, grp, gnm);

    R.effects.contraction = mk('A',  'Main effect: contraction type', 'Contraction', ...
        {[1 2 5 6],[3 4 7 8]}, {'CON','ECC'});
    R.effects.intensity   = mk('B',  'Main effect: intensity', 'Intensity', ...
        {[1 2 3 4],[5 6 7 8]}, {'75% MVC','90% MVC'});
    R.effects.repetition  = mk('C',  'Main effect: repetition bin', 'Repetition', ...
        {[1 3 5 7],[2 4 6 8]}, {'Early','Late'});
    R.effects.AxB         = mk('AB', 'Interaction: contraction x intensity', 'AxB', ...
        {[1 2],[3 4],[5 6],[7 8]}, {'CON 75','ECC 75','CON 90','ECC 90'});
    R.effects.AxC         = mk('AC', 'Interaction: contraction x repetition', 'AxC', ...
        {[1 5],[2 6],[3 7],[4 8]}, {'CON early','CON late','ECC early','ECC late'});
    R.effects.BxC         = mk('BC', 'Interaction: intensity x repetition', 'BxC', ...
        {[1 3],[2 4],[5 7],[6 8]}, {'75 early','75 late','90 early','90 late'});
    R.effects.AxBxC       = mk('ABC','Interaction: contraction x intensity x repetition', ...
        'AxBxC', num2cell(1:8), D.condNames);
end

function [a,b] = ab_of(cond)
    switch cond
        case 1, a = 1; b = 1;   % CON 75
        case 2, a = 2; b = 1;   % ECC 75
        case 3, a = 1; b = 2;   % CON 90
        case 4, a = 2; b = 2;   % ECC 90
        otherwise, error('Unknown condition index.');
    end
end

function E = build(spmi, xAxis, varLabel, units, effectLabel, effectName, groups, groupNames)
    E = struct('varLabel', varLabel, 'units', units, ...
        'effectLabel', effectLabel, 'effectName', effectName, 'xAxis', xAxis, ...
        'spmi', spmi, 'z', getf(spmi,'z'), 'zstar', getf(spmi,'zstar'), ...
        'h0reject', getf(spmi,'h0reject'), 'df', getf(spmi,'df'), ...
        'clusters', clusters_of(spmi, xAxis), 'statLabel', 'SPM{F}', ...
        'groups', {groups}, 'groupNames', {groupNames});
end

function val = getf(S, fieldName)
% spm1d exposes some quantities as properties and some as fields, depending on
% the version, so both are tried.
    val = [];
    try
        if isprop(S, fieldName) || isfield(S, fieldName)
            val = S.(fieldName);
        end
    catch
    end
end

function C = clusters_of(spmi, xAxis)
% Converts each suprathreshold cluster from spm1d node coordinates into percent
% of the analysis window.
%
% The whole body is wrapped in a try block so that a version mismatch returns an
% empty cluster list rather than aborting the run. The cost is that a genuine
% error here would also look like "no clusters found".
    C = struct('startPct',{}, 'endPct',{}, 'extentPct',{}, 'P',{});
    try
        if isempty(spmi.clusters), return; end
        nT = numel(xAxis);
        nodeAxis = 0:(nT-1);

        for k = 1:numel(spmi.clusters)
            % Clusters come back as a cell array in some versions and a struct
            % array in others.
            try
                cl = spmi.clusters{k};
            catch
                cl = spmi.clusters(k);
            end

            s1 = NaN; s2 = NaN; ex = NaN; pv = NaN;

            try
                xy = cl.endpoints;
                if all(isfinite(xy))
                    s1 = interp1(nodeAxis, xAxis, xy(1), 'linear','extrap');
                    s2 = interp1(nodeAxis, xAxis, xy(2), 'linear','extrap');
                    ex = s2 - s1;
                end
            catch
                % No interpolated endpoints, so fall back to sample indices.
                try
                    ind = max(min(round(cl.indices(:)), nT-1), 0);
                    s1 = xAxis(min(ind)+1); s2 = xAxis(max(ind)+1); ex = s2 - s1;
                catch
                end
            end

            % Older versions do not expose a per-cluster P value.
            try
                pv = cl.P;
            catch
            end

            C(k).startPct = s1; C(k).endPct = s2;
            C(k).extentPct = ex; C(k).P = pv;
        end
    catch
    end
end

function print_summary(R)
    fprintf('\n--- %s | %s ---\n', R.varLabel, R.effectLabel);
    if ~isempty(R.zstar),    fprintf('  critical threshold = %.3f\n', R.zstar); end
    if ~isempty(R.h0reject), fprintf('  reject H0: %d\n', logical(R.h0reject)); end
    if isempty(R.clusters)
        fprintf('  no suprathreshold clusters\n');
    else
        for k = 1:numel(R.clusters)
            fprintf('  cluster %d: %.1f%% to %.1f%% (extent %.1f%%), p = %.4f\n', ...
                k, R.clusters(k).startPct, R.clusters(k).endPct, ...
                R.clusters(k).extentPct, R.clusters(k).P);
        end
    end
end

function report_range(tag, R)
% Prints the range the SPM trace actually occupies, so shared y-limits can be
% chosen for the composite panels rather than guessed.
    if isempty(R.z)
        fprintf('%-20s no trace\n', tag);
        return;
    end
    zs = NaN;
    if ~isempty(R.zstar), zs = R.zstar; end
    fprintf('%-20s trace %.1f to %.1f | threshold %.2f\n', ...
        tag, min(R.z), max(R.z), zs);
end

function cols = palette_for(groupNames)
% Colour assignment matched to Stages 6C1 and 6C4, so a factor always looks the
% same across every figure in the set. CON blue and ECC green; 75% purple and
% 90% orange; early light blue and late dark blue. The four-group early and late
% case uses two hues in dark and light pairs, so contraction mode reads as hue
% and repetition bin as lightness.
    base = [0.00 0.4470 0.7410;   % 1 blue
            0.4660 0.6740 0.1880; % 2 green
            0.8500 0.3250 0.0980; % 3 orange
            0.4940 0.1840 0.5560; % 4 purple
            0.3010 0.7450 0.9330; % 5 light blue
            0.7290 0.8730 0.5410; % 6 light green
            0.9290 0.6940 0.1250; % 7 gold
            0.6350 0.0780 0.1840];% 8 dark red

    nG = numel(groupNames);
    g  = string(groupNames(:)).';

    if nG == 2 && all(g == ["CON","ECC"])
        cols = base([1 2],:);
    elseif nG == 2 && all(g == ["75% MVC","90% MVC"])
        cols = base([4 3],:);
    elseif nG == 2 && all(g == ["Early","Late"])
        cols = [0.40 0.65 0.85; 0.05 0.25 0.45];
    elseif nG == 4 && all(contains(g, ["early","late"]))
        cols = [0.20 0.20 0.20;   % CON early, near black
                0.60 0.60 0.60;   % CON late,  grey
                0.65 0.16 0.16;   % ECC early, dark red
                0.93 0.55 0.38];  % ECC late,  salmon
    else
        cols = base(1:min(nG,8),:);
        if nG > 8, cols = [cols; base(1:(nG-8),:)]; end
    end
end

function [M, S] = traces_for(D, groups, isCentroid)
    nG = numel(groups);
    Yg = cell(1,nG);
    for g = 1:nG
        acc = 0;
        for c = groups{g}, acc = acc + D.Y{c}; end
        Yg{g} = acc / numel(groups{g});
    end

    if isCentroid
        ref = 0;
        for g = 1:nG, ref = ref + mean(Yg{g},2,'omitnan'); end
        ref = ref / nG;
        gm = mean(ref,'omitnan');
        for g = 1:nG, Yg{g} = Yg{g} - ref + gm; end
    end

    M = cell(1,nG); S = cell(1,nG);
    for g = 1:nG
        M{g} = mean(Yg{g},1,'omitnan');
        S{g} = std(Yg{g},0,1,'omitnan');
    end
end

function fig = plot_result(R, D, showTitles, STY)
% Standalone figure: one trace panel above one SPM panel.
    fig = figure('Color','w','Name', sprintf('%s | %s', R.varLabel, R.effectLabel));
    tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
    draw_spm_pair(tl, R, D, showTitles, STY, struct());
end

function draw_spm_pair(tl, R, D, showTitles, STY, opts)
% Draws the trace panel and the SPM panel into an existing 2x1 tiledlayout.
% opts fields, all optional:
%   panelLetter  : 'A', 'B', ... drawn inside the trace panel, top left
%   yLim         : [lo hi] forced limits for the trace panel
%   yTickStep    : tick spacing for the trace panel when yLim is given
%   spmYLim      : [lo hi] forced limits for the SPM panel
%   spmYTickStep : tick spacing for the SPM panel, ticks anchored at zero
%   showXLabel   : logical, default true
%   showLegend   : logical, default true
%   showXTicks   : logical, default true
%   showYTicks   : logical, default true; false also blanks the axis label

    if ~isfield(opts,'panelLetter'),  opts.panelLetter  = ''; end
    if ~isfield(opts,'yLim'),         opts.yLim         = []; end
    if ~isfield(opts,'yTickStep'),    opts.yTickStep    = []; end
    if ~isfield(opts,'spmYLim'),      opts.spmYLim      = []; end
    if ~isfield(opts,'spmYTickStep'), opts.spmYTickStep = []; end
    if ~isfield(opts,'showXLabel'),   opts.showXLabel   = true; end
    if ~isfield(opts,'showLegend'),   opts.showLegend   = true; end
    if ~isfield(opts,'showXTicks'),   opts.showXTicks   = true; end
    if ~isfield(opts,'showYTicks'),   opts.showYTicks   = true; end

    tickLen = [STY.TickLen STY.TickLen];

    [M, S] = traces_for(D, R.groups, D.isCentroid);
    cols   = palette_for(R.groupNames);

    % ---------- upper panel ----------
    ax1 = nexttile(tl,1); hold(ax1,'on');

    for g = 1:numel(M)
        fill(ax1, [D.xAxis fliplr(D.xAxis)], [M{g}-S{g} fliplr(M{g}+S{g})], ...
            cols(g,:), 'FaceAlpha',STY.FaceAlpha, 'EdgeColor','none', ...
            'HandleVisibility','off');
    end

    h = gobjects(1,numel(M));
    for g = 1:numel(M)
        h(g) = plot(ax1, D.xAxis, M{g}, 'Color', cols(g,:), 'LineWidth', STY.TraceWidth);
    end

    if D.isCentroid
        yline(ax1, 0, ':k', 'LineWidth',2, 'HandleVisibility','off');
    end

    ylabel(ax1, sprintf('%s [%s]', R.varLabel, R.units), ...
        'FontWeight','bold','FontSize',STY.FontSize);
    if showTitles
        title(ax1, sprintf('%s | %s', R.varLabel, R.effectLabel), ...
            'Interpreter','none','FontWeight','bold','FontSize',STY.FontSize);
    end
    if opts.showLegend
        legend(ax1, h, R.groupNames, 'Box','off','Location','southeast', ...
            'FontSize',STY.FontSize);
    end
    box(ax1,'off');
    set(ax1,'TickDir','out','LineWidth',STY.AxLineWidth, ...
        'TickLength',tickLen, ...
        'FontWeight','bold','FontSize',STY.FontSize,'XTick',[]);
    xlim(ax1,[D.xAxis(1) D.xAxis(end)]);

    if ~isempty(opts.yLim)
        ylim(ax1, opts.yLim);
        if ~isempty(opts.yTickStep)
            yticks(ax1, opts.yLim(1):opts.yTickStep:opts.yLim(2));
        end
    elseif D.isCentroid
        yl = ylim(ax1);
        Lm = ceil(max(abs(yl)));
        ylim(ax1, [-Lm Lm]);
        yticks(ax1, -Lm:1:Lm);
    end

    if ~opts.showYTicks
        set(ax1,'YTick',[]);
        ylabel(ax1,'');
    end

    if ~isempty(opts.panelLetter)
        text(ax1, 0.01, 0.97, opts.panelLetter, 'Units','normalized', ...
            'FontSize',STY.FontSize+4, 'FontWeight','bold', ...
            'HorizontalAlignment','left','VerticalAlignment','top');
    end

    % ---------- lower panel ----------
    ax2 = nexttile(tl,2); hold(ax2,'on');

    for k = 1:numel(R.clusters)
        x1 = R.clusters(k).startPct; x2 = R.clusters(k).endPct;
        if isfinite(x1) && isfinite(x2)
            xregion(ax2, x1, x2, 'FaceColor',[0.7 0.7 0.7], 'FaceAlpha',0.30);
        end
    end

    if ~isempty(R.z)
        plot(ax2, R.xAxis, R.z, 'k', 'LineWidth', STY.SPMWidth);
    end

    if ~isempty(R.zstar)
        yline(ax2, R.zstar, '--r', 'LineWidth',2, 'HandleVisibility','off');
        text(ax2, R.xAxis(end), R.zstar, '*', 'Color','k', ...
            'FontSize',18, 'FontWeight','bold', ...
            'HorizontalAlignment','left','VerticalAlignment','middle', ...
            'Clipping','off');
        if contains(R.statLabel,'t')
            yline(ax2, -R.zstar, '--r', 'LineWidth',2, 'HandleVisibility','off');
            text(ax2, R.xAxis(end), -R.zstar, '*', 'Color','k', ...
                'FontSize',18, 'FontWeight','bold', ...
                'HorizontalAlignment','left','VerticalAlignment','middle', ...
                'Clipping','off');
        end
    end

    if opts.showXLabel
        xlabel(ax2,'Constant velocity phase [%]', ...
            'FontWeight','bold','FontSize',STY.FontSize);
    end
    ylabel(ax2, R.statLabel, 'Interpreter','none', ...
        'FontWeight','bold','FontSize',STY.FontSize);
    box(ax2,'off');
    set(ax2,'TickDir','out','LineWidth',STY.AxLineWidth, ...
        'TickLength',tickLen, ...
        'FontWeight','bold','FontSize',STY.FontSize);
    xlim(ax2,[R.xAxis(1) R.xAxis(end)]);

    if ~isempty(opts.spmYLim)
        ylim(ax2, opts.spmYLim);
        if ~isempty(opts.spmYTickStep)
            % Ticks are anchored at zero rather than at the axis floor, which
            % sits slightly below zero so the trace clears the axis line.
            lo = opts.spmYLim(1); hi = opts.spmYLim(2);
            st = opts.spmYTickStep;
            if lo < 0
                tk = unique([fliplr(0:-st:lo) 0:st:hi]);
            else
                tk = 0:st:hi;
            end
            yticks(ax2, tk);
        end
    end

    if ~opts.showXTicks
        set(ax2,'XTick',[]);
    end
    if ~opts.showYTicks
        set(ax2,'YTick',[]);
        ylabel(ax2,'');
    end
end

function fig = scalar_figure(Yearly, Ylate, condNames, outcomeName, showTitles, STY)
% Early against late scalar panel, styled to match the Stage 6B scalar figures.
% Individual participants are connected within each condition, with the group
% mean and between-participant SD drawn over them.

    labs = struct('torque','Torque [Nm]', 'emg','RMS EMG [% MVA]', ...
                  'cx','Centroid X displacement [mm]', ...
                  'cy','Centroid Y displacement [mm]');
    ylab = labs.(outcomeName);

    fig = figure('Color','w','Name', sprintf('Early vs late | %s', outcomeName));
    ax = axes(fig); hold(ax,'on');

    xPos = [1 2 4 5 7 8 10 11];
    Y = nan(size(Yearly,1), 8);
    for c = 1:4
        Y(:,2*c-1) = Yearly(:,c);
        Y(:,2*c)   = Ylate(:,c);
    end

    nSub = size(Y,1);
    subjColors = distinguishable_colors(nSub, {'w','k'});
    subjColors = subjColors(1:nSub,:);

    lineColor    = [0.75 0.75 0.75];
    lineWidthInd = 0.5;
    markerSize   = 100;

    % --- individual participants ---
    for i = 1:nSub
        for k = 1:2:8
            yi = Y(i,k:k+1);
            if all(~isnan(yi))
                plot(ax, xPos(k:k+1), yi, '-', ...
                    'Color', lineColor, 'LineWidth', lineWidthInd);
            end
            for j = 0:1
                if ~isnan(Y(i,k+j))
                    scatter(ax, xPos(k+j), Y(i,k+j), markerSize, ...
                        'MarkerFaceColor', subjColors(i,:), ...
                        'MarkerEdgeColor','k', 'LineWidth',0.4);
                end
            end
        end
    end

    % --- group mean and SD ---
    m = mean(Y,1,'omitnan');
    s = std(Y,0,1,'omitnan');

    capW = 0.25;
    for c = 1:8
        if ~isnan(m(c)) && ~isnan(s(c))
            plot(ax, [xPos(c) xPos(c)], [m(c)-s(c), m(c)+s(c)], '-', ...
                'Color',[0.35 0.35 0.35], 'LineWidth',2.0);
            plot(ax, xPos(c)+[-capW capW], [m(c)-s(c) m(c)-s(c)], '-', ...
                'Color',[0.35 0.35 0.35], 'LineWidth',2.0);
            plot(ax, xPos(c)+[-capW capW], [m(c)+s(c) m(c)+s(c)], '-', ...
                'Color',[0.35 0.35 0.35], 'LineWidth',2.0);
        end
    end

    for k = 1:2:8
        plot(ax, xPos(k:k+1), m(k:k+1), '--', 'Color','k', 'LineWidth',1.2);
    end
    scatter(ax, xPos, m, 42, ...
        'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',0.5);

    % --- axis limits ---
    % Centroid panels share a fixed range with headroom for the condition
    % labels; the others scale to their own data.
    switch outcomeName
        case {'cx','cy'}
            ylim(ax, [-7 8]);
            yticks(ax, -6:2:6);
            yTop = 6; yBot = -7;
            labelY = 7.3;
        otherwise
            yTop = max(Y(:),[],'omitnan'); yBot = min(Y(:),[],'omitnan');
            rng  = max(yTop - yBot, eps);
            ylim(ax, [yBot - 0.12*rng, yTop + 0.18*rng]);
            labelY = yTop + 0.10*rng;
    end

    % --- condition labels above each early and late pair ---
    for i = 1:4
        text(ax, mean(xPos(2*i-1:2*i)), labelY, ...
            strrep(char(condNames(i)),'_',' '), ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold','FontSize',STY.FontSize);
    end

    xlim(ax, [0.5 11.5]);
    set(ax, 'XTick', xPos, 'XTickLabel', {'E','L','E','L','E','L','E','L'}, ...
        'FontSize',STY.FontSize, 'FontWeight','bold', ...
        'LineWidth',STY.AxLineWidth, 'TickDir','out');
    ylabel(ax, ylab, 'FontWeight','bold','FontSize',STY.FontSize);
    if showTitles
        title(ax, sprintf('%s | early vs late', ylab), ...
            'Interpreter','none','FontWeight','bold','FontSize',STY.FontSize);
    end
    box(ax,'off');

    if yBot < 0 && yTop > 0
        yline(ax, 0, ':k', 'LineWidth',2);
    end
end
