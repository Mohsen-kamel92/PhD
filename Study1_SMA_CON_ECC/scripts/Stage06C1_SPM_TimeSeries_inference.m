%% ===== Stage 06C1: Time-series inference with spm1d over the analysis window =====
%
% PURPOSE
%   Test whether the outcomes differ between contraction modes and intensities
%   at any point across the analysis window, rather than only on average.
%
% SUPPORTS
%   The time-series inference reported in the Results, Figure.3 of the
%   manuscript, and the centroid supplementary figures.
%
% INPUT
%   <ROOT>/groupData_stage6A.mat                      (Stage 06A)
%
% OUTPUT
%   <ROOT>/stage6C1_spmResults_mid20to80.mat
%   <ROOT>/Figures_Stage6C1/*.png
%   <ROOT>/Figures_Stage6C1/Stage6C1_SFig_centroid_composite.{png,pdf}
%   <ROOT>/stage6C1_clusters.csv
%
% DESIGN
%   2 x 2 repeated measures SPM ANOVA on mean EMG, centroid X and centroid Y.
%   Factor A is contraction type (CON, ECC), factor B is intensity (75, 90).
%   Paired SPM t-tests follow at each intensity.
%
%   Torque is not analysed here. It is a manipulation check rather than an
%   outcome, and the scalar comparison in Stage 06B already establishes the
%   degree of matching.
%
% WHY ECC IS REVERSED BEFORE INFERENCE
%   The two modes traverse the range of motion in opposite directions. CON runs
%   from about 37 to -3 deg and ECC from about -4 to 37 deg, so comparing them
%   sample by sample without reversal would contrast different joint angles at
%   every point. Reversing ECC in time aligns them by angle, so 20% of the
%   window corresponds to roughly 29 deg in both. The window is symmetric about
%   50%, so reversing after cropping is equivalent to reversing first.
%
% WHICH CONDITIONS EACH PANEL SHOWS
%   The upper panel shows only the conditions the test below it compares:
%     main effect of contraction : CON and ECC, collapsed across intensity
%     main effect of intensity   : 75 and 90, collapsed across contraction
%     interaction                : all four conditions
%     paired t-test              : the two conditions at that intensity
%
% ON THE SHADED BANDS
%   Centroid position differs by up to 11 mm between participants, because the
%   grid was positioned relative to each participant's own muscle-tendon
%   junction and their anatomy differed. For centroid outcomes the bands are
%   therefore computed after removing each participant's own mean across the
%   conditions shown, so they reflect within-subject variability. EMG keeps
%   conventional between-subject SDs.
%
% ON MULTIPLE COMPARISONS
%   SPM controls family-wise error across the time domain within each test
%   through random field theory. No further correction is applied across the
%   three outcomes, which are treated as separate hypotheses.
%
% ON APPROXIMATE RESIDUALS
%   With one observation per participant per condition, spm1d computes
%   residuals approximately and prints a warning to that effect. This is
%   expected for this design and is noted in the manuscript.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped.
%
% DEPENDENCIES
%   spm1d for MATLAB (https://spm1d.org)

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'mainBin',            'all');
CFG = set_default(CFG, 'midWin',             [20 80]);   % analysis window, % of CV
CFG = set_default(CFG, 'alphaLevel',         0.05);
CFG = set_default(CFG, 'interpClusters',     true);
CFG = set_default(CFG, 'flipECCforAnalysis', true);
CFG = set_default(CFG, 'doPlots',            true);
CFG = set_default(CFG, 'saveFigures',        true);
CFG = set_default(CFG, 'doPostHoc',          true);
CFG = set_default(CFG, 'doNormality',        true);
CFG = set_default(CFG, 'showTitles',         true);

% Composite supplementary figure, section (10b)
CFG = set_default(CFG, 'comp_cxYLim',  [-2 2]);      % centroid X trace panels
CFG = set_default(CFG, 'comp_cyYLim',  [-2 2]);      % centroid Y trace panels
CFG = set_default(CFG, 'comp_yStep',   1);           % mm
CFG = set_default(CFG, 'comp_spmYLim', [-3 60]);     % SPM panels
CFG = set_default(CFG, 'comp_spmStep', 20);

% Styling, matched to Stage 6B
STY = struct('FontSize',12, 'AxLineWidth',2.5, 'TraceWidth',2.2, ...
             'SPMWidth',1.8, 'FaceAlpha',0.15);

mainBin            = CFG.mainBin;
midWin             = CFG.midWin;
alphaLevel         = CFG.alphaLevel;
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
        error('Stage 6C1: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

spmPath = which('spm1d.stats.anova2rm');
assert(~isempty(spmPath), ...
    ['Stage 6C1: spm1d not found on the MATLAB path. ' ...
     'Add it with addpath(genpath(''path_to_spm1d'')).']);

groupFile = fullfile(ROOT,'groupData_stage6A.mat');
assert(isfile(groupFile), ...
    'Stage 6C1: %s not found. Run Stage 6A first.', groupFile);

L = load(groupFile);
assert(isfield(L,'groupData'), 'File does not contain "groupData".');
groupData = L.groupData;

condNames     = string(groupData.conditions(:));
expectedOrder = ["CON_75"; "ECC_75"; "CON_90"; "ECC_90"];
assert(all(condNames == expectedOrder), ...
    'Condition order mismatch. Expected [CON_75, ECC_75, CON_90, ECC_90].');

nSub       = groupData.nSub;
xAxis_full = groupData.xAxis(:).';
midMask    = xAxis_full >= midWin(1) & xAxis_full <= midWin(2);
assert(any(midMask), 'No samples inside the requested window.');

xAxis = xAxis_full(midMask);
nT    = numel(xAxis);

figDir = fullfile(ROOT,'Figures_Stage6C1');
if CFG.saveFigures && ~exist(figDir,'dir'), mkdir(figDir); end

fprintf('\n===== Stage 6C1: time-series inference =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('spm1d found at %s\n', spmPath);
fprintf('Participants %d | bin %s | alpha %.3f | window %d-%d %%CV | ECC reversed %d\n', ...
    nSub, mainBin, alphaLevel, midWin(1), midWin(2), flipECCforAnalysis);

%% ---------------------------
%  (3) Extract and align the time series
%% ---------------------------
T_emg   = groupData.timeseries.(mainBin).emgMeanNorm_101(:,:,midMask);
T_cx    = groupData.timeseries.(mainBin).centroidXc_101(:,:,midMask);
T_cy    = groupData.timeseries.(mainBin).centroidYc_101(:,:,midMask);
T_angle = groupData.timeseries.(mainBin).angle_101(:,:,midMask);

for V = {T_emg,T_cx,T_cy}
    assert(isequal(size(V{1}),[nSub 4 nT]), 'Unexpected time-series size.');
end

aCON = []; aECC = [];
if flipECCforAnalysis
    for c = [2 4]
        T_emg(:,c,:)   = flip(T_emg(:,c,:),   3);
        T_cx(:,c,:)    = flip(T_cx(:,c,:),    3);
        T_cy(:,c,:)    = flip(T_cy(:,c,:),    3);
        T_angle(:,c,:) = flip(T_angle(:,c,:), 3);
    end

    aCON = squeeze(mean(mean(T_angle(:,[1 3],:),1,'omitnan'),2,'omitnan')).';
    aECC = squeeze(mean(mean(T_angle(:,[2 4],:),1,'omitnan'),2,'omitnan')).';
    fprintf('Angle after alignment: CON %.1f to %.1f deg | ECC %.1f to %.1f deg\n', ...
        aCON(1), aCON(end), aECC(1), aECC(end));
    fprintf('Mean absolute angle mismatch across the window: %.2f deg\n', ...
        mean(abs(aCON - aECC)));
end

%% ---------------------------
%  (4) Output structure
%% ---------------------------
spmResults = struct();
spmResults.meta = struct( ...
    'createdOn', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'script', mfilename, 'sourceFile', groupFile, 'mainBin', mainBin, ...
    'nSub', nSub, 'nTime', nT, 'xAxis', xAxis, 'xAxis_label', '%CV', ...
    'alpha', alphaLevel, 'interpClusters', interpClusters, ...
    'flipECCforAnalysis', flipECCforAnalysis, 'analysisWindow', midWin, ...
    'design', '2x2 repeated-measures ANOVA', ...
    'outcomes', {{'emg','cx','cy'}}, ...
    'factorA', 'Contraction', 'factorALevels', {{'CON','ECC'}}, ...
    'factorB', 'Intensity',   'factorBLevels', {{'75','90'}}, ...
    'torque_note', ["torque is a manipulation check and is tested only as a " ...
                    "scalar in Stage 6B"], ...
    'alignment_note', ["ECC reversed in time so both modes are compared at " ...
                       "matched joint angles rather than matched elapsed time"], ...
    'mc_note', ["SPM controls family-wise error across the time domain within " ...
                "each test; no correction across the three outcomes"]);

if flipECCforAnalysis
    spmResults.meta.angleCON = aCON;
    spmResults.meta.angleECC = aECC;
end

%% ---------------------------
%  (5) Normality of the paired differences
%% ---------------------------
if CFG.doNormality
    fprintf('\n===== NORMALITY (spm1d, paired differences) =====\n');
    spmResults.normality = struct();
    outc = {'emg',T_emg; 'cx',T_cx; 'cy',T_cy};

    for i = 1:size(outc,1)
        T = outc{i,2};
        for I = 1:2
            if I == 1, a = 1; b = 2; lab = "75"; else, a = 3; b = 4; lab = "90"; end
            D = squeeze(T(:,a,:)) - squeeze(T(:,b,:));
            try
                spmNi = spm1d.stats.normality.ttest(D).inference(alphaLevel);
                rej = logical(spmNi.h0reject);
                fprintf('%-6s %s%%: normality rejected = %d\n', outc{i,1}, lab, rej);
                spmResults.normality.(outc{i,1}).("rej"+lab) = rej;
            catch ME
                fprintf('%-6s %s%%: normality test unavailable (%s)\n', ...
                    outc{i,1}, lab, ME.message);
            end
        end
    end
end

%% ---------------------------
%  (6) SPM RM-ANOVA
%% ---------------------------
spmResults.emg = run_spm_anova2rm(T_emg, xAxis, alphaLevel, interpClusters, 'RMS EMG',   '% MVA', false);
spmResults.cx  = run_spm_anova2rm(T_cx,  xAxis, alphaLevel, interpClusters, 'Centroid X','mm',    true);
spmResults.cy  = run_spm_anova2rm(T_cy,  xAxis, alphaLevel, interpClusters, 'Centroid Y','mm',    true);

fprintf('\n===== SPM RM-ANOVA SUMMARY =====\n');
outcomes = {'emg','cx','cy'};
for i = 1:numel(outcomes)
    R = spmResults.(outcomes{i});
    print_spm_summary(R.effects.contraction);
    print_spm_summary(R.effects.intensity);
    print_spm_summary(R.effects.interaction);
end

%% ---------------------------
%  (7) Paired SPM t-tests
%% ---------------------------
if CFG.doPostHoc
    fprintf('\n===== PAIRED SPM t-TESTS (CON vs ECC at each intensity) =====\n');
    data = struct('emg',T_emg, 'cx',T_cx, 'cy',T_cy);

    for i = 1:numel(outcomes)
        nm = outcomes{i};
        T  = data.(nm);
        R  = spmResults.(nm);

        spmResults.(nm).posthoc = struct();
        for I = 1:2
            if I == 1, a = 1; b = 2; lab = "75"; else, a = 3; b = 4; lab = "90"; end

            spmTi = local_infer(spm1d.stats.ttest_paired( ...
                squeeze(T(:,a,:)), squeeze(T(:,b,:))), alphaLevel, interpClusters);

            E = struct('varLabel', R.varLabel, 'units', R.units, ...
                'effectLabel', sprintf('CON vs ECC at %s%% MVC', lab), ...
                'effectName', "posthoc_"+lab, 'xAxis', xAxis, ...
                'spmi', spmTi, ...
                'z', getf(spmTi,'z'), 'zstar', getf(spmTi,'zstar'), ...
                'h0reject', getf(spmTi,'h0reject'), 'df', getf(spmTi,'df'), ...
                'clusters', summarize_clusters(spmTi, xAxis), ...
                'statLabel', 'SPM{t}', ...
                'showConds', [a b], 'collapse', 'none');

            spmResults.(nm).posthoc.("i"+lab) = E;
            print_spm_summary(E);
        end
    end
end

%% ---------------------------
%  (8) Cluster table
%% ---------------------------
rows = {};
for i = 1:numel(outcomes)
    R = spmResults.(outcomes{i});
    eff = {R.effects.contraction, R.effects.intensity, R.effects.interaction};
    if CFG.doPostHoc && isfield(R,'posthoc')
        pf = fieldnames(R.posthoc);
        for j = 1:numel(pf), eff{end+1} = R.posthoc.(pf{j}); end %#ok<SAGROW>
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

clusterTable = cell2table(rows, 'VariableNames', ...
    {'outcome','effect','cluster','startPct','endPct','extentPct','p'});
writetable(clusterTable, fullfile(ROOT,'stage6C1_clusters.csv'));

%% ---------------------------
%  (9) Figures, one per effect
%% ---------------------------
% Every effect is drawn and saved. Which of them reach the supplement is decided
% in section (10); keeping the full set here means a result can be produced
% later without rerunning the analysis.

if CFG.doPlots
    figs = struct();
    for i = 1:numel(outcomes)
        nm = outcomes{i};
        R  = spmResults.(nm);
        figs.([nm '_contraction']) = plot_spm_result(R.effects.contraction, R.descriptives, CFG.showTitles, STY);
        figs.([nm '_intensity'])   = plot_spm_result(R.effects.intensity,   R.descriptives, CFG.showTitles, STY);
        figs.([nm '_interaction']) = plot_spm_result(R.effects.interaction, R.descriptives, CFG.showTitles, STY);

        if CFG.doPostHoc && isfield(R,'posthoc')
            pf = fieldnames(R.posthoc);
            for j = 1:numel(pf)
                figs.([nm '_' char(pf{j})]) = ...
                    plot_spm_result(R.posthoc.(pf{j}), R.descriptives, CFG.showTitles, STY);
            end
        end
    end

    if CFG.saveFigures
        names = fieldnames(figs);
        for i = 1:numel(names)
            exportgraphics(figs.(names{i}), ...
                fullfile(figDir, sprintf('Stage6C1_%s.png', names{i})), 'Resolution', 300);
        end
        fprintf('\nFigures saved to %s\n', figDir);
    end
end

%% ---------------------------
%  (10) Composite supplementary figure: centroid, both main effects
%% ---------------------------
%  A  centroid X, main effect of contraction mode
%  B  centroid Y, main effect of contraction mode
%  C  centroid X, main effect of intensity
%  D  centroid Y, main effect of intensity
%
%  X panels share one y-range and Y panels share another, so A can be compared
%  against C and B against D.

if CFG.doPlots
    figC = figure('Color','w','Name','Supplementary: centroid SPM', ...
                  'Position',[100 100 1250 950]);
    outer = tiledlayout(figC,2,2,'TileSpacing','compact','Padding','compact');

    %            effect                            descriptives              letter  yLim   xLab  legend xTick yTick
    panels = { spmResults.cx.effects.contraction, spmResults.cx.descriptives, 'A', CFG.comp_cxYLim, false, true, false, true;
               spmResults.cy.effects.contraction, spmResults.cy.descriptives, 'B', CFG.comp_cyYLim, false, true, false, false;
               spmResults.cx.effects.intensity,   spmResults.cx.descriptives, 'C', CFG.comp_cxYLim, true,  true, true,  true;
               spmResults.cy.effects.intensity,   spmResults.cy.descriptives, 'D', CFG.comp_cyYLim, true,  true, true,  false };

    for k = 1:size(panels,1)
        inner = tiledlayout(outer,2,1,'TileSpacing','compact','Padding','tight');
        inner.Layout.Tile = k;

        opts = struct('panelLetter',  panels{k,3}, ...
                      'yLim',         panels{k,4}, ...
                      'yTickStep',    CFG.comp_yStep, ...
                      'spmYLim',      CFG.comp_spmYLim, ...
                      'spmYTickStep', CFG.comp_spmStep, ...
                      'showXLabel',   panels{k,5}, ...
                      'showLegend',   panels{k,6}, ...
                      'showXTicks',   panels{k,7}, ...
                      'showYTicks',   panels{k,8});

        draw_spm_pair(inner, panels{k,1}, panels{k,2}, false, STY, opts);
    end

    if CFG.saveFigures
        exportgraphics(figC, fullfile(figDir,'Stage6C1_SFig_centroid_composite.png'), ...
            'Resolution', 300);
        exportgraphics(figC, fullfile(figDir,'Stage6C1_SFig_centroid_composite.pdf'), ...
            'ContentType','vector');
        fprintf('Composite centroid figure saved to %s\n', figDir);
    end
end

%% ---------------------------
%  (11) Save
%% ---------------------------
outName = fullfile(ROOT, 'stage6C1_spmResults_mid20to80.mat');
save(outName, 'spmResults', '-v7.3');
fprintf('Saved results : %s\n', outName);
fprintf('Saved clusters: %s\n', fullfile(ROOT,'stage6C1_clusters.csv'));

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

function R = run_spm_anova2rm(T, xAxis, alphaLevel, interpClusters, varLabel, units, isCentroid)
% Reshapes [nSub x 4 x nT] into the long form spm1d expects, runs the two-way
% repeated measures ANOVA, and packages each effect with its descriptives.
    [nSub, nCond, nT] = size(T);
    assert(nCond == 4, 'T must hold 4 conditions in dimension 2.');

    Y = zeros(nSub*4, nT);
    A = zeros(nSub*4, 1);
    B = zeros(nSub*4, 1);
    SUBJ = zeros(nSub*4, 1);

    row = 0;
    for s = 1:nSub
        for c = 1:4
            row = row + 1;
            Y(row,:) = squeeze(T(s,c,:)).';
            switch c
                case 1, A(row) = 1; B(row) = 1;   % CON 75
                case 2, A(row) = 2; B(row) = 1;   % ECC 75
                case 3, A(row) = 1; B(row) = 2;   % CON 90
                case 4, A(row) = 2; B(row) = 2;   % ECC 90
            end
            SUBJ(row) = s;
        end
    end

    spmList = spm1d.stats.anova2rm(Y, A, B, SUBJ);
    spmiA  = local_infer(spmList('A'),  alphaLevel, interpClusters);
    spmiB  = local_infer(spmList('B'),  alphaLevel, interpClusters);
    spmiAB = local_infer(spmList('AB'), alphaLevel, interpClusters);

    D = struct('xAxis', xAxis, 'varLabel', varLabel, 'units', units, ...
        'isCentroid', isCentroid, ...
        'condNames', {{'CON 75','ECC 75','CON 90','ECC 90'}}, ...
        'Y', {cell(1,4)});
    for c = 1:4
        D.Y{c} = squeeze(T(:,c,:));
    end

    R = struct('varLabel', varLabel, 'units', units, 'xAxis', xAxis, 'descriptives', D);

    R.effects.contraction = build_effect(spmiA,  xAxis, varLabel, units, ...
        'Main effect: contraction type', 'Contraction', [1 2 3 4], 'intensity');
    R.effects.intensity   = build_effect(spmiB,  xAxis, varLabel, units, ...
        'Main effect: intensity', 'Intensity', [1 2 3 4], 'contraction');
    R.effects.interaction = build_effect(spmiAB, xAxis, varLabel, units, ...
        'Interaction: contraction x intensity', 'Interaction', [1 2 3 4], 'none');
end

function spmi = local_infer(spmObj, alphaLevel, interpClusters)
    try
        spmi = spmObj.inference(alphaLevel, 'interp', interpClusters);
    catch
        spmi = spmObj.inference(alphaLevel);
    end
end

function E = build_effect(spmi, xAxis, varLabel, units, effectLabel, effectName, showConds, collapse)
    E = struct('varLabel', varLabel, 'units', units, ...
        'effectLabel', effectLabel, 'effectName', effectName, 'xAxis', xAxis, ...
        'spmi', spmi, 'z', getf(spmi,'z'), 'zstar', getf(spmi,'zstar'), ...
        'h0reject', getf(spmi,'h0reject'), 'df', getf(spmi,'df'), ...
        'clusters', summarize_clusters(spmi, xAxis), ...
        'statLabel', 'SPM{F}', 'showConds', showConds, 'collapse', collapse);
end

function val = getf(S, fieldName)
    val = [];
    try
        if isprop(S, fieldName) || isfield(S, fieldName)
            val = S.(fieldName);
        end
    catch
    end
end

function C = summarize_clusters(spmi, xAxis)
% Converts each suprathreshold cluster from spm1d node coordinates into percent
% of the analysis window. The endpoints are interpolated where spm1d provides
% them, and fall back to sample indices where it does not.
    C = struct('startPct',{}, 'endPct',{}, 'extentPct',{}, 'P',{});
    try
        if isempty(spmi.clusters), return; end
        nT = numel(xAxis);
        nodeAxis = 0:(nT-1);

        for k = 1:numel(spmi.clusters)
              % spm1d returns clusters as a cell array in some versions and a
            % struct array in others.
            try
                cl = spmi.clusters{k};
            catch
                cl = spmi.clusters(k);
            end
            startPct = NaN; endPct = NaN; extentPct = NaN; pval = NaN;

            try
                xy = cl.endpoints;
                if all(isfinite(xy))
                    startPct  = interp1(nodeAxis, xAxis, xy(1), 'linear','extrap');
                    endPct    = interp1(nodeAxis, xAxis, xy(2), 'linear','extrap');
                    extentPct = endPct - startPct;
                end
            catch
                try
                    inds = max(min(round(cl.indices(:)), nT-1), 0);
                    startPct  = xAxis(min(inds)+1);
                    endPct    = xAxis(max(inds)+1);
                    extentPct = endPct - startPct;
                catch
                end
            end

             % Older spm1d versions do not expose a per-cluster P value, in
            % which case pval stays NaN.
            try
                pval = cl.P;
            catch
            end

            C(k).startPct = startPct;  C(k).endPct = endPct;
            C(k).extentPct = extentPct; C(k).P = pval;
        end
    catch
    end
end

function print_spm_summary(R)
    fprintf('\n--- %s | %s ---\n', R.varLabel, R.effectLabel);
    if ~isempty(R.zstar)
        fprintf('  critical threshold = %.3f\n', R.zstar);
    end
    if ~isempty(R.h0reject)
        fprintf('  reject H0: %d\n', logical(R.h0reject));
    end
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

function [M, S, names] = build_traces(D, showConds, collapse, isCentroid)
% Returns the mean and SD traces the panel should show, given which effect is
% being tested. Bands for centroid outcomes are computed after removing each
% participant's own mean across the conditions shown, so they reflect
% within-subject rather than grid-placement variability.

    switch lower(string(collapse))
        case "intensity"    % main effect of contraction: CON against ECC
            G = {[1 3], [2 4]};
            names = {'CON','ECC'};
        case "contraction"  % main effect of intensity: 75 against 90
            G = {[1 2], [3 4]};
            names = {'75% MVC','90% MVC'};
        otherwise
            G = num2cell(showConds);
            names = D.condNames(showConds);
    end

    nG = numel(G);
    Yg = cell(1,nG);
    for g = 1:nG
        acc = 0;
        for c = G{g}, acc = acc + D.Y{c}; end
        Yg{g} = acc / numel(G{g});
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

function fig = plot_spm_result(R, D, showTitles, STY)
% Standalone figure: one trace panel above one SPM panel.
    fig = figure('Color','w','Name', sprintf('SPM | %s | %s', R.varLabel, R.effectLabel));
    tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
    draw_spm_pair(tl, R, D, showTitles, STY, struct());
end

function draw_spm_pair(tl, R, D, showTitles, STY, opts)
% Draws the trace panel and the SPM panel into an existing 2x1 tiledlayout.
% opts fields, all optional:
%   panelLetter  : 'A', 'B', ... drawn inside the trace panel, top left
%   yLim         : [lo hi] forced limits for the trace panel
%   yTickStep    : tick spacing for the trace panel when yLim is given
%   spmYLim      : [lo hi] forced limits for the SPM panel
%   spmYTickStep : tick spacing for the SPM panel, ticks start at zero
%   showXLabel   : logical, default true
%   showLegend   : logical, default true
%   showXTicks   : logical, default true
%   showYTicks   : logical, default true

    if ~isfield(opts,'panelLetter'),  opts.panelLetter  = ''; end
    if ~isfield(opts,'yLim'),         opts.yLim         = []; end
    if ~isfield(opts,'yTickStep'),    opts.yTickStep    = []; end
    if ~isfield(opts,'spmYLim'),      opts.spmYLim      = []; end
    if ~isfield(opts,'spmYTickStep'), opts.spmYTickStep = []; end
    if ~isfield(opts,'showXLabel'),   opts.showXLabel   = true; end
    if ~isfield(opts,'showLegend'),   opts.showLegend   = true; end
    if ~isfield(opts,'showXTicks'),   opts.showXTicks   = true; end
    if ~isfield(opts,'showYTicks'),   opts.showYTicks   = true; end

    [M, S, names] = build_traces(D, R.showConds, R.collapse, D.isCentroid);

    % CON blue and ECC green; 75% purple and 90% orange, so a factor always
    % looks the same across every figure in the set.
    cols = [0.00 0.4470 0.7410;
            0.4660 0.6740 0.1880;
            0.8500 0.3250 0.0980;
            0.4940 0.1840 0.5560];
    if strcmpi(string(R.collapse),"contraction")   % main effect of intensity
        cols = cols([4 3],:);
    end

    % ---------- upper panel ----------
    ax1 = nexttile(tl,1); hold(ax1,'on');

    h = gobjects(1,numel(M));
    for g = 1:numel(M)
        fill(ax1, [D.xAxis fliplr(D.xAxis)], [M{g}-S{g} fliplr(M{g}+S{g})], ...
            cols(g,:), 'FaceAlpha',STY.FaceAlpha, 'EdgeColor','none', ...
            'HandleVisibility','off');
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
        legend(ax1, h, names, 'Box','off','Location','southeast','FontSize',STY.FontSize);
    end
    box(ax1,'off');
    set(ax1,'TickDir','out','LineWidth',STY.AxLineWidth, ...
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
    else
        ylim(ax1, [20 100]);
        yticks(ax1, 20:20:100);
    end

    if ~opts.showYTicks
        set(ax1,'YTick',[]);
    end

    if ~isempty(opts.panelLetter)
        text(ax1, 0.01, 0.97, opts.panelLetter, 'Units','normalized', ...
            'FontSize',STY.FontSize+4, 'FontWeight','bold', ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end

    % ---------- lower panel ----------
    ax2 = nexttile(tl,2); hold(ax2,'on');

    if ~isempty(R.clusters)
        for k = 1:numel(R.clusters)
            x1 = R.clusters(k).startPct; x2 = R.clusters(k).endPct;
            if isfinite(x1) && isfinite(x2)
                xregion(ax2, x1, x2, 'FaceColor',[0.7 0.7 0.7], 'FaceAlpha',0.30);
            end
        end
    end

    if ~isempty(R.z)
        plot(ax2, R.xAxis, R.z, 'k', 'LineWidth', STY.SPMWidth);
    end

    if ~isempty(R.zstar)
        yline(ax2, R.zstar, '--r', 'LineWidth',2, 'HandleVisibility','off');
        text(ax2, R.xAxis(end), R.zstar, '*', 'Color','k', ...
            'FontSize',18, 'FontWeight','bold', ...
            'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
            'Clipping','off');
        if contains(R.statLabel,'t')
            yline(ax2, -R.zstar, '--r', 'LineWidth',2, 'HandleVisibility','off');
            text(ax2, R.xAxis(end), -R.zstar, '*', 'Color','k', ...
                'FontSize',18, 'FontWeight','bold', ...
                'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
                'Clipping','off');
        end
    end

    if opts.showXLabel
        xlabel(ax2, 'Constant velocity phase [%]', ...
            'FontWeight','bold','FontSize',STY.FontSize);
    end
    ylabel(ax2, R.statLabel, 'Interpreter','none', ...
        'FontWeight','bold','FontSize',STY.FontSize);
    box(ax2,'off');
    set(ax2,'TickDir','out','LineWidth',STY.AxLineWidth, ...
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
    end
end
