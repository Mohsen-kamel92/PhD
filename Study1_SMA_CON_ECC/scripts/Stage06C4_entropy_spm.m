%% ===== Stage 06C4: Shannon spatial entropy over the analysis window =====
%
% PURPOSE
%   Quantify how evenly amplitude is spread across the electrode grid, and test
%   whether that evenness differs between contraction modes or intensities.
%
% SUPPORTS
%   The entropy results reported in the Results, Figure 5 of the manuscript,
%   and the entropy supplementary figure.
%
% INPUT
%   <ROOT>/groupData_stage6A.mat           close all;           (Stage 06A)
%
% OUTPUT
%   <ROOT>/stage6C4_entropy_mid20to80.mat
%   <ROOT>/Figures_Stage6C4/*.png
%   <ROOT>/Figures_Stage6C4/Stage6C4_SFig_entropy_composite.{png,pdf}
%   <ROOT>/stage6C4_results.csv
%   <ROOT>/stage6C4_clusters.csv
%
% WHAT ENTROPY MEASURES
%   How evenly amplitude is spread across the bipolar channels at each instant.
%   Values approach 100% when every channel contributes equally and fall as
%   activity concentrates in fewer channels. It is therefore independent of
%   overall amplitude and complements the centroid, which describes where
%   activity sits rather than how evenly it is spread.
%
%       E = -sum(p_i * ln p_i) / ln(N),   p_i = w_i / sum(w)
%
%   where N is the number of channels.
%
% ON SCALE INVARIANCE
%   Because p_i is a proportion, multiplying every channel by the same constant
%   leaves entropy unchanged. With the single scalar normalization used at Stage
%   3, entropy computed from %MVA is identical to entropy computed from raw RMS.
%   This would not hold under per-channel normalization, where each channel
%   receives a different multiplier and the proportions themselves are altered.
%
% ON TERMINOLOGY
%   Farina et al. (2008) and subsequent work describe this quantity as modified
%   entropy. It is reported here as normalized Shannon entropy, which is the
%   same measure expressed as a percentage of its maximum.
%
% ON THE SHADED BANDS
%   Entropy varies more between participants than between conditions, so a
%   conventional SD band would obscure the comparison. In the time-series
%   figures each participant's own mean across the conditions shown is removed
%   before computing the band, leaving within-subject variability. The scalar
%   figure keeps between-subject SD, since it plots individual participants and
%   that is the appropriate descriptor there.
%
% WHY ECC IS REVERSED
%   Same convention as Stages 6C1 and 6C2. The modes traverse the range of
%   motion in opposite directions, so ECC is reversed in time to compare them at
%   matched joint angles. Applied to both the analysis and the figures. The
%   scalar means are unaffected by the reversal.
%
% ON APPROXIMATE RESIDUALS
%   With one observation per participant per condition, spm1d computes residuals
%   approximately and prints a warning to that effect. This is expected for this
%   design and is noted in the manuscript.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped. The
%   scalar analysis runs without spm1d; only the time-series inference needs it.
%
% DEPENDENCIES
%   Statistics and Machine Learning Toolbox   fitrm, ranova, ttest
%   distinguishable_colors.m  participant colours  (MATLAB File Exchange)
%   spm1d for MATLAB (optional, https://spm1d.org)

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'binName',            'all');
CFG = set_default(CFG, 'midWin',             [20 80]);   % analysis window, % of CV
CFG = set_default(CFG, 'alphaLevel',         0.05);
CFG = set_default(CFG, 'interpClusters',     true);
CFG = set_default(CFG, 'flipECCforAnalysis', true);
CFG = set_default(CFG, 'doPlots',            true);
CFG = set_default(CFG, 'saveFigures',        true);
CFG = set_default(CFG, 'showTitles',         true);

% Scalar figure axis, fixed so the panel is comparable across runs
CFG = set_default(CFG, 'entYLim',   [88 100]);
CFG = set_default(CFG, 'entStep',   2);

% Composite supplementary figure, section (8)
CFG = set_default(CFG, 'comp_entYLim',  [88 100]);
CFG = set_default(CFG, 'comp_entStep',  4);
CFG = set_default(CFG, 'comp_spmYLim',  [-0.8 20]);
CFG = set_default(CFG, 'comp_spmStep',  10);
CFG = set_default(CFG, 'comp_spmStep',  10);

% Styling, matched to Stages 6B, 6C1 and 6C2
STY = struct('FontSize',12, 'AxLineWidth',2.5, 'TraceWidth',2.2, ...
             'SPMWidth',1.8, 'FaceAlpha',0.15, 'TickLen',0.008);

binName            = CFG.binName;
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
        error('Stage 6C4: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

assert(~isempty(which('distinguishable_colors')), ...
    'Stage 6C4: distinguishable_colors not found. Add the /functions folder.');

groupFile = fullfile(ROOT,'groupData_stage6A.mat');
assert(isfile(groupFile), ...
    'Stage 6C4: %s not found. Run Stage 6A first.', groupFile);

L = load(groupFile);
assert(isfield(L,'groupData'), 'File does not contain "groupData".');
groupData = L.groupData;

condNames     = string(groupData.conditions(:));
expectedOrder = ["CON_75"; "ECC_75"; "CON_90"; "ECC_90"];
assert(all(condNames == expectedOrder), 'Condition order mismatch.');

nSub  = groupData.nSub;
nCond = groupData.nCond;
nT    = groupData.nTime;
xAxis = groupData.xAxis(:).';

% Channel count taken from the data rather than assumed, since it sets the
% normalizing constant ln(N) in the entropy formula.
if isfield(groupData,'nBipolar')
    nCh = groupData.nBipolar;
elseif isfield(groupData,'mapInfo') && isfield(groupData.mapInfo,'bipolar_rowcol')
    nCh = size(groupData.mapInfo.bipolar_rowcol,1);
else
    nCh = 28;
end

figDir = fullfile(ROOT,'Figures_Stage6C4');
if CFG.saveFigures && ~exist(figDir,'dir'), mkdir(figDir); end

fprintf('\n===== Stage 6C4: Shannon spatial entropy =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Participants %d | %d channels | window %d-%d %%CV\n', ...
    nSub, nCh, midWin(1), midWin(2));

%% ---------------------------
%  (3) Entropy time series
%% ---------------------------
B = groupData.bipolar.(binName).emgRMS_bipNorm_101;
assert(isequal(size(B), [nSub nCond nCh nT]), 'Unexpected bipolar array size.');

rawH    = nan(nSub,nCond,nT);
entropy = nan(nSub,nCond,nT);   % fraction of maximum
entPct  = nan(nSub,nCond,nT);   % percent of maximum

for s = 1:nSub
    for c = 1:nCond
        W = squeeze(B(s,c,:,:));
        W(~isfinite(W)) = NaN;
        W(W < 0) = 0;                       % an RMS cannot be negative

        for t = 1:nT
            w = W(:,t);
            w = w(isfinite(w) & w >= 0);
            if isempty(w), continue; end

            tot = sum(w,'omitnan');
            if tot <= 0, continue; end

            pi_ = w ./ tot;
            pi_ = pi_(pi_ > 0 & isfinite(pi_));   % 0*ln(0) is defined as zero

            H = -sum(pi_ .* log(pi_), 'omitnan');
            rawH(s,c,t)    = H;
            entropy(s,c,t) = H / log(nCh);
            entPct(s,c,t)  = 100 * H / log(nCh);
        end
    end
end

nMissing = nnz(~isfinite(entPct));
if nMissing > 0
    fprintf('Note: %d of %d entropy samples could not be computed.\n', ...
        nMissing, numel(entPct));
end

%% ---------------------------
%  (4) Align ECC and take window means
%% ---------------------------
entPct_aligned = entPct;
if flipECCforAnalysis
    for c = [2 4]
        entPct_aligned(:,c,:) = flip(entPct_aligned(:,c,:), 3);
    end
end

midMask   = xAxis >= midWin(1) & xAxis <= midWin(2);
assert(any(midMask), 'No samples inside the requested window.');
xAxis_mid = xAxis(midMask);

Y = mean(entPct_aligned(:,:,midMask), 3, 'omitnan');

%% ---------------------------
%  (5) Scalar RM-ANOVA and follow-up tests
%% ---------------------------
entropyStats = struct();
entropyStats.meta = struct( ...
    'createdOn', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'script', mfilename, 'sourceFile', groupFile, 'binName', binName, ...
    'midWin', midWin, 'xAxis_mid_percentCV', xAxis_mid, ...
    'alpha', alphaLevel, 'nChannels', nCh, ...
    'flipECCforAnalysis', flipECCforAnalysis, ...
    'formula', 'E = -sum(p_i ln p_i) / ln(N), expressed as % of maximum', ...
    'invariance_note', ["entropy is invariant to a global scale factor, so " ...
                        "with the single scalar normalization used at Stage 3 " ...
                        "it is identical whether computed from %MVA or raw RMS"], ...
    'terminology_note', ["equivalent to the modified entropy of Farina et al. " ...
                         "(2008), expressed as a percentage of its maximum"], ...
    'band_note', ["time-series bands are within-subject: each participant's " ...
                  "own mean across the conditions shown is removed first"], ...
    'sphericity_note', ["two levels per factor, so sphericity holds by " ...
                        "definition and no correction is applicable"]);

entropyStats.timeseries = struct('rawH_101', rawH, ...
    'entropy_101', entropy, 'entropy_percent_101', entPct, ...
    'entropy_percent_aligned_101', entPct_aligned);
entropyStats.scalar.Y_entropy_mid = Y;

entropyStats.scalar.rm = run_rm2x2(Y, groupData.participants, 'ShannonEntropy');

entropyStats.posthoc.CONvsECC_75 = paired_compare(Y(:,1), Y(:,2), 'Entropy CON75 vs ECC75');
entropyStats.posthoc.CONvsECC_90 = paired_compare(Y(:,3), Y(:,4), 'Entropy CON90 vs ECC90');

pH = holm([entropyStats.posthoc.CONvsECC_75.p, entropyStats.posthoc.CONvsECC_90.p]);
entropyStats.posthoc.CONvsECC_75.p_holm = pH(1);
entropyStats.posthoc.CONvsECC_90.p_holm = pH(2);

fprintf('\n===== SCALAR ENTROPY, %d-%d%% CV =====\n', midWin(1), midWin(2));
print_rm_result(entropyStats.scalar.rm, 'Shannon spatial entropy');

fprintf('\nPaired comparisons (Holm corrected):\n');
fprintf('%-26s %8s %8s %9s %8s %8s\n','comparison','CON','ECC','diff','p','p_holm');
resRows = {};
for k = ["75","90"]
    H = entropyStats.posthoc.("CONvsECC_"+k);
    fprintf('%-26s %8.2f %8.2f %+9.3f %8.4f %8.4f  dz = %+.2f\n', ...
        H.label, H.mean1, H.mean2, H.meanDiff, H.p, H.p_holm, H.dz);
    resRows(end+1,:) = {H.label, H.mean1, H.sd1, H.mean2, H.sd2, H.meanDiff, ...
        H.stats.tstat, H.stats.df, H.p, H.p_holm, H.dz, H.ci(1), H.ci(2)}; %#ok<SAGROW>
end

writetable(cell2table(resRows, 'VariableNames', ...
    {'comparison','CON_mean','CON_sd','ECC_mean','ECC_sd','meanDiff', ...
     't','df','p','p_holm','dz','ci_lo','ci_hi'}), ...
    fullfile(ROOT,'stage6C4_results.csv'));

%% ---------------------------
%  (6) SPM time-series inference
%% ---------------------------
doSPM = ~isempty(which('spm1d.stats.anova2rm'));
if ~doSPM
    warning('spm1d not found. The scalar analysis is complete; SPM skipped.');
else
    entropyStats.spm = run_spm(entPct_aligned(:,:,midMask), xAxis_mid, ...
        alphaLevel, interpClusters);

    fprintf('\n===== SPM, ENTROPY =====\n');
    for e = {'contraction','intensity','interaction'}
        print_spm_summary(entropyStats.spm.effects.(e{1}));
    end

    rows = {};
    for e = {'contraction','intensity','interaction'}
        E = entropyStats.spm.effects.(e{1});
        if isempty(E.clusters)
            rows(end+1,:) = {string(E.effectLabel), 0, NaN, NaN, NaN, NaN}; %#ok<SAGROW>
        else
            for k = 1:numel(E.clusters)
                rows(end+1,:) = {string(E.effectLabel), k, ...
                    E.clusters(k).startPct, E.clusters(k).endPct, ...
                    E.clusters(k).extentPct, E.clusters(k).P}; %#ok<SAGROW>
            end
        end
    end
    writetable(cell2table(rows, 'VariableNames', ...
        {'effect','cluster','startPct','endPct','extentPct','p'}), ...
        fullfile(ROOT,'stage6C4_clusters.csv'));
end

%% ---------------------------
%  (7) Figures, one per effect
%% ---------------------------
% Every effect is drawn and saved. Which of them reach the supplement is decided
% in section (8); keeping the full set here means a result can be produced later
% without rerunning the analysis.

if CFG.doPlots
    figs = struct();

    figs.scalar = figure('Color','w','Name','Entropy | scalar');
    plot_scalar_panel(axes(figs.scalar), Y, ...
        ttl(CFG.showTitles, sprintf('Shannon spatial entropy | %d-%d%% CV', ...
        midWin(1), midWin(2))), 'Normalized spatial entropy [%]', CFG, STY);

    figs.timeseries = figure('Color','w','Name','Entropy | time series');
    tl = tiledlayout(figs.timeseries,1,2,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl);
    plot_pair(ax1, xAxis, squeeze(entPct_aligned(:,1,:)), squeeze(entPct_aligned(:,2,:)), ...
        'CON 75','ECC 75','Spatial entropy [%]', ...
        ttl(CFG.showTitles,'75% MVC'), midWin, STY);

    ax2 = nexttile(tl);
    plot_pair(ax2, xAxis, squeeze(entPct_aligned(:,3,:)), squeeze(entPct_aligned(:,4,:)), ...
        'CON 90','ECC 90','', ttl(CFG.showTitles,'90% MVC'), midWin, STY);

    link_ylim([ax1 ax2]);

    if doSPM
        for e = {'contraction','intensity','interaction'}
            figs.(['spm_' e{1}]) = plot_spm(entropyStats.spm.effects.(e{1}), ...
                entropyStats.spm.descriptives, CFG.showTitles, STY);
        end
    end

    if CFG.saveFigures
        names = fieldnames(figs);
        for i = 1:numel(names)
            exportgraphics(figs.(names{i}), ...
                fullfile(figDir, sprintf('Stage6C4_%s.png', names{i})), 'Resolution', 300);
        end
        fprintf('\nFigures saved to %s\n', figDir);
    end
end

%% ---------------------------
%  (8) Composite supplementary figure
%% ---------------------------
%  A  main effect of contraction mode
%  B  main effect of intensity

if CFG.doPlots && doSPM
    figE = figure('Color','w','Name','Supplementary: entropy SPM', ...
                  'Position',[100 100 900 950]);
    outerE = tiledlayout(figE,2,1,'TileSpacing','compact','Padding','compact');

    panelsE = { entropyStats.spm.effects.contraction, 'A', false, false;
                entropyStats.spm.effects.intensity,   'B', true,  true };

    for k = 1:size(panelsE,1)
        innerE = tiledlayout(outerE,2,1,'TileSpacing','compact','Padding','tight');
        innerE.Layout.Tile = k;

        optsE = struct('panelLetter',  panelsE{k,2}, ...
                       'yLim',         CFG.comp_entYLim, ...
                       'yTickStep',    CFG.comp_entStep, ...
                       'spmYLim',      CFG.comp_spmYLim, ...
                       'spmYTickStep', CFG.comp_spmStep, ...
                       'showXLabel',   panelsE{k,3}, ...
                       'showXTicks',   panelsE{k,4});

        draw_entropy_pair(innerE, panelsE{k,1}, entropyStats.spm.descriptives, ...
            false, STY, optsE);
    end

    if CFG.saveFigures
        exportgraphics(figE, fullfile(figDir,'Stage6C4_SFig_entropy_composite.png'), ...
            'Resolution', 300);
        exportgraphics(figE, fullfile(figDir,'Stage6C4_SFig_entropy_composite.pdf'), ...
            'ContentType','vector');
        fprintf('Composite entropy figure saved to %s\n', figDir);
    end
end

%% ---------------------------
%  (9) Save
%% ---------------------------
outName = fullfile(ROOT,'stage6C4_entropy_mid20to80.mat');
save(outName, 'entropyStats', '-v7.3');
fprintf('Saved results : %s\n', outName);
fprintf('Saved scalars : %s\n', fullfile(ROOT,'stage6C4_results.csv'));
if doSPM
    fprintf('Saved clusters: %s\n', fullfile(ROOT,'stage6C4_clusters.csv'));
end

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

       % A folder holding only the processed group file is a valid root for the
    % Stage 6 analyses, which read nothing else.
    hasGroupFile = isfile(fullfile(ROOT,'groupData_stage6A.mat'));

    assert(hasParticipant || hasGroupFile, ...
        ['The selected folder does not look like the expected data root.\n' ...
         'Expected either the raw layout\n' ...
         '  <ROOT>/<participant>/<CON_75|CON_90|ECC_75|ECC_90>/<Set_N>/\n' ...
         'or a folder containing groupData_stage6A.mat\n' ...
         'Selected: %s'], ROOT);
end

function t = ttl(show, txt)
    if show, t = txt; else, t = ''; end
end

function link_ylim(axes_)
    lims = cell2mat(get(axes_,'YLim'));
    set(axes_,'YLim',[min(lims(:,1)) max(lims(:,2))]);
end

function pAdj = holm(p)
    [ps, idx] = sort(p(:));
    m = numel(ps);
    adj = min(cummax(ps .* (m - (1:m)' + 1)), 1);
    pAdj = nan(size(ps));
    pAdj(idx) = adj;
    pAdj = reshape(pAdj, size(p));
end

function Yout = within_subject(Ycell)
% Removes each participant's own mean across the groups supplied, so the
% resulting SD reflects within-subject variability. The grand mean is added back
% to keep the values on their original scale.
    nG = numel(Ycell);
    ref = 0;
    for g = 1:nG, ref = ref + mean(Ycell{g}, 2, 'omitnan'); end
    ref = ref / nG;
    gm = mean(ref, 'omitnan');

    Yout = cell(1,nG);
    for g = 1:nG
        Yout{g} = Ycell{g} - ref + gm;
    end
end

function S = run_rm2x2(Y, participants, varLabel)
    T = table(string(participants(:)), Y(:,1), Y(:,2), Y(:,3), Y(:,4), ...
        'VariableNames', {'Subject','CON_75','ECC_75','CON_90','ECC_90'});

    within = table(categorical({'CON';'ECC';'CON';'ECC'}), ...
                   categorical({'75';'75';'90';'90'}), ...
                   'VariableNames', {'Contraction','Intensity'});

    rm  = fitrm(T, 'CON_75-ECC_90 ~ 1', 'WithinDesign', within);
    tbl = ranova(rm, 'WithinModel', 'Contraction*Intensity');

    S = struct('label', varLabel, 'data', Y, 'ranova', tbl, ...
        'mean', mean(Y,1,'omitnan'), 'sd', std(Y,0,1,'omitnan'), ...
        'df_effect', 1, 'df_error', size(Y,1)-1);

    S.p_contraction = getval(tbl, "(Intercept):Contraction", 'pValue');
    S.p_intensity   = getval(tbl, "(Intercept):Intensity",   'pValue');
    S.p_interaction = getval(tbl, "(Intercept):Contraction:Intensity", 'pValue');

    S.F_contraction = getval(tbl, "(Intercept):Contraction", 'F');
    S.F_intensity   = getval(tbl, "(Intercept):Intensity",   'F');
    S.F_interaction = getval(tbl, "(Intercept):Contraction:Intensity", 'F');

    S.eta2p_contraction = eta2p(tbl, "(Intercept):Contraction", "Error(Contraction)");
    S.eta2p_intensity   = eta2p(tbl, "(Intercept):Intensity",   "Error(Intensity)");
    S.eta2p_interaction = eta2p(tbl, "(Intercept):Contraction:Intensity", ...
                                     "Error(Contraction:Intensity)");
end

function v = getval(tbl, rowName, colName)
    v = NaN;
    i = find(string(tbl.Properties.RowNames) == rowName, 1);
    if ~isempty(i) && any(strcmp(tbl.Properties.VariableNames, colName))
        v = tbl.(colName)(i);
    end
end

function v = eta2p(tbl, effRow, errRow)
% Partial eta squared from the effect and its own error term.
    v = NaN;
    rn = string(tbl.Properties.RowNames);
    ie = find(rn == effRow, 1);
    ir = find(rn == errRow, 1);
    if ~isempty(ie) && ~isempty(ir)
        v = tbl.SumSq(ie) / (tbl.SumSq(ie) + tbl.SumSq(ir));
    end
end

function print_rm_result(S, labelText)
    fprintf('\n--- %s ---\n', labelText);
    fprintf('  means [CON75 ECC75 CON90 ECC90]: %s\n', num2str(S.mean,'%8.2f'));
    fprintf('  SDs                            : %s\n', num2str(S.sd,  '%8.2f'));
    fprintf('  Contraction : F(%d,%d) = %7.2f, p = %.5f, eta2p = %.3f\n', ...
        S.df_effect, S.df_error, S.F_contraction, S.p_contraction, S.eta2p_contraction);
    fprintf('  Intensity   : F(%d,%d) = %7.2f, p = %.5f, eta2p = %.3f\n', ...
        S.df_effect, S.df_error, S.F_intensity, S.p_intensity, S.eta2p_intensity);
    fprintf('  Interaction : F(%d,%d) = %7.2f, p = %.5f, eta2p = %.3f\n', ...
        S.df_effect, S.df_error, S.F_interaction, S.p_interaction, S.eta2p_interaction);
end

function H = paired_compare(y1, y2, labelText)
    [~, pv, ci, st] = ttest(y1, y2);
    d = y1 - y2;
    H = struct('label', labelText, ...
        'mean1', mean(y1,'omitnan'), 'mean2', mean(y2,'omitnan'), ...
        'sd1', std(y1,0,'omitnan'), 'sd2', std(y2,0,'omitnan'), ...
        'meanDiff', mean(d,'omitnan'), 'sdDiff', std(d,0,'omitnan'), ...
        'dz', mean(d,'omitnan')/std(d,0,'omitnan'), ...
        'ci', ci(:).', 'p', pv, 'p_holm', NaN, 'stats', st);
end

function R = run_spm(T, xAxis, alphaLevel, interpClusters)
% Reshapes [nSub x 4 x nT] into the long form spm1d expects, runs the two-way
% repeated measures ANOVA, and packages each effect with its descriptives.
    [nSub, nCond, nT] = size(T);
    assert(nCond == 4);

    Y = zeros(nSub*4, nT);
    A = zeros(nSub*4,1); B = zeros(nSub*4,1); S = zeros(nSub*4,1);

    row = 0;
    for s = 1:nSub
        for c = 1:4
            row = row + 1;
            Y(row,:) = squeeze(T(s,c,:)).';
            switch c
                case 1, A(row)=1; B(row)=1;   % CON 75
                case 2, A(row)=2; B(row)=1;   % ECC 75
                case 3, A(row)=1; B(row)=2;   % CON 90
                case 4, A(row)=2; B(row)=2;   % ECC 90
            end
            S(row) = s;
        end
    end

    sl = spm1d.stats.anova2rm(Y, A, B, S);

    D = struct('xAxis', xAxis, 'varLabel', 'Spatial entropy', 'units', '%', ...
        'condNames', {{'CON 75','ECC 75','CON 90','ECC 90'}}, 'Y', {cell(1,4)});
    for c = 1:4, D.Y{c} = squeeze(T(:,c,:)); end

    R = struct('descriptives', D);
    R.effects.contraction = build(infer(sl('A'), alphaLevel, interpClusters), ...
        xAxis, 'Main effect: contraction type', {[1 3],[2 4]}, {'CON','ECC'});
    R.effects.intensity   = build(infer(sl('B'), alphaLevel, interpClusters), ...
        xAxis, 'Main effect: intensity', {[1 2],[3 4]}, {'75% MVC','90% MVC'});
    R.effects.interaction = build(infer(sl('AB'), alphaLevel, interpClusters), ...
        xAxis, 'Interaction: contraction x intensity', {1,2,3,4}, D.condNames);
end

function spmi = infer(spmObj, alphaLevel, interpClusters)
% Older spm1d versions do not accept the interp option, so the call falls back
% to the plain form.
    try
        spmi = spmObj.inference(alphaLevel, 'interp', interpClusters);
    catch
        spmi = spmObj.inference(alphaLevel);
    end
end

function E = build(spmi, xAxis, effectLabel, groups, groupNames)
    E = struct('varLabel', 'Spatial entropy', 'units', '%', ...
        'effectLabel', effectLabel, 'xAxis', xAxis, 'spmi', spmi, ...
        'z', getf(spmi,'z'), 'zstar', getf(spmi,'zstar'), ...
        'h0reject', getf(spmi,'h0reject'), 'df', getf(spmi,'df'), ...
        'clusters', clusters_of(spmi, xAxis), 'statLabel','SPM{F}', ...
        'groups', {groups}, 'groupNames', {groupNames});
end

function val = getf(S, fieldName)
% spm1d exposes some quantities as properties and some as fields, depending on
% the version, so both are tried.
    val = [];
    try
        if isprop(S, fieldName) || isfield(S, fieldName), val = S.(fieldName); end
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

            C(k).startPct = s1; C(k).endPct = s2; C(k).extentPct = ex; C(k).P = pv;
        end
    catch
    end
end

function print_spm_summary(R)
    fprintf('\n--- %s ---\n', R.effectLabel);
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

function plot_scalar_panel(ax, Y, ttlText, ylab, CFG, STY)
% Individual participants connected within each intensity, with the group mean
% and between-participant SD drawn over them. Between-subject SD is the right
% descriptor here, because the panel plots the individuals themselves.
    hold(ax,'on');
    xPos = [1 2 4 5];
    nSub = size(Y,1);
    cols = distinguishable_colors(nSub, {'w','k'});

    for i = 1:nSub
        plot(ax, xPos(1:2), Y(i,1:2), '-', 'Color',[0.75 0.75 0.75], 'LineWidth',0.6);
        plot(ax, xPos(3:4), Y(i,3:4), '-', 'Color',[0.75 0.75 0.75], 'LineWidth',0.6);
        for c = 1:4
            scatter(ax, xPos(c), Y(i,c), 60, 'MarkerFaceColor', cols(i,:), ...
                'MarkerEdgeColor','k', 'LineWidth',0.4);
        end
    end

    m = mean(Y,1,'omitnan');
    s = std(Y,0,1,'omitnan');
    capW = 0.25;
    for c = 1:4
        plot(ax, [xPos(c) xPos(c)], [m(c)-s(c) m(c)+s(c)], '-', ...
            'Color',[0.35 0.35 0.35], 'LineWidth',1.0);
        plot(ax, xPos(c)+[-capW capW], [m(c)-s(c) m(c)-s(c)], '-', ...
            'Color',[0.35 0.35 0.35], 'LineWidth',1.5);
        plot(ax, xPos(c)+[-capW capW], [m(c)+s(c) m(c)+s(c)], '-', ...
            'Color',[0.35 0.35 0.35], 'LineWidth',1.5);
    end
    plot(ax, xPos(1:2), m(1:2), '--k', 'LineWidth',1.2);
    plot(ax, xPos(3:4), m(3:4), '--k', 'LineWidth',1.2);
    scatter(ax, xPos, m, 42, 'MarkerFaceColor','k','MarkerEdgeColor','k');

    xlim(ax,[0.4 5.6]);
    ylim(ax, CFG.entYLim);
    yticks(ax, CFG.entYLim(1):CFG.entStep:CFG.entYLim(2));

    set(ax,'XTick',xPos,'XTickLabel',{'CON75','ECC75','CON90','ECC90'}, ...
        'FontSize',STY.FontSize,'FontWeight','bold', ...
        'LineWidth',STY.AxLineWidth,'TickDir','out', ...
        'TickLength',[STY.TickLen STY.TickLen]);
    ylabel(ax, ylab, 'FontWeight','bold','FontSize',STY.FontSize);
    if ~isempty(ttlText)
        title(ax, ttlText, 'FontWeight','bold','Interpreter','none', ...
            'FontSize',STY.FontSize);
    end
    box(ax,'off');
end

function plot_pair(ax, xAxis, Y1, Y2, name1, name2, ylab, ttlText, winCV, STY)
% CON against ECC across the whole constant velocity phase, with the analysis
% window marked.
    hold(ax,'on');

    Yp = within_subject({Y1, Y2});

    m1 = mean(Yp{1},1,'omitnan'); s1 = std(Yp{1},0,1,'omitnan');
    m2 = mean(Yp{2},1,'omitnan'); s2 = std(Yp{2},0,1,'omitnan');

    col1 = [0.00 0.4470 0.7410];   % CON
    col2 = [0.4660 0.6740 0.1880]; % ECC

    % Entropy is bounded at 100% of its theoretical maximum, so the upper band is
    % clipped there. Mean plus SD can exceed the ceiling even when no individual
    % observation does.
    up1 = min(m1+s1, 100);
    up2 = min(m2+s2, 100);

    fill(ax,[xAxis fliplr(xAxis)],[m1-s1 fliplr(up1)], col1, ...
        'FaceAlpha',STY.FaceAlpha,'EdgeColor','none','HandleVisibility','off');
    fill(ax,[xAxis fliplr(xAxis)],[m2-s2 fliplr(up2)], col2, ...
        'FaceAlpha',STY.FaceAlpha,'EdgeColor','none','HandleVisibility','off');
    h1 = plot(ax, xAxis, m1, 'Color', col1, 'LineWidth', STY.TraceWidth);
    h2 = plot(ax, xAxis, m2, 'Color', col2, 'LineWidth', STY.TraceWidth);

    xline(ax, winCV(1), '--k','LineWidth',2,'HandleVisibility','off');
    xline(ax, winCV(2), '--k','LineWidth',2,'HandleVisibility','off');

    xlabel(ax,'Constant velocity phase [%]','FontWeight','bold','FontSize',STY.FontSize);
    if ~isempty(ylab)
        ylabel(ax, ylab, 'FontWeight','bold','FontSize',STY.FontSize);
    end
    if ~isempty(ttlText)
        title(ax, ttlText, 'FontWeight','bold','FontSize',STY.FontSize);
    end

    legend(ax,[h1 h2],{name1,name2},'Box','off','Location','best', ...
        'FontSize',STY.FontSize);
    box(ax,'off');
    set(ax,'TickDir','out','LineWidth',STY.AxLineWidth, ...
        'TickLength',[STY.TickLen STY.TickLen], ...
        'FontWeight','bold','FontSize',STY.FontSize);
    xlim(ax,[xAxis(1) xAxis(end)]);
end

function fig = plot_spm(R, D, showTitles, STY)
% Standalone figure: one trace panel above one SPM panel.
    fig = figure('Color','w','Name',sprintf('Entropy | %s', R.effectLabel));
    tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
    draw_entropy_pair(tl, R, D, showTitles, STY, struct());
end

function draw_entropy_pair(tl, R, D, showTitles, STY, opts)
% Draws the trace panel and the SPM panel into an existing 2x1 tiledlayout.
% opts fields, all optional:
%   panelLetter  : 'A', 'B', ... drawn inside the trace panel, top left
%   yLim         : [lo hi] forced limits for the trace panel
%   yTickStep    : tick spacing for the trace panel when yLim is given
%   spmYLim      : [lo hi] forced limits for the SPM panel
%   spmYTickStep : tick spacing for the SPM panel, ticks anchored at zero
%   showXLabel   : logical, default true
%   showXTicks   : logical, default true

    if ~isfield(opts,'panelLetter'),  opts.panelLetter  = ''; end
    if ~isfield(opts,'yLim'),         opts.yLim         = []; end
    if ~isfield(opts,'yTickStep'),    opts.yTickStep    = []; end
    if ~isfield(opts,'spmYLim'),      opts.spmYLim      = []; end
    if ~isfield(opts,'spmYTickStep'), opts.spmYTickStep = []; end
    if ~isfield(opts,'showXLabel'),   opts.showXLabel   = true; end
    if ~isfield(opts,'showXTicks'),   opts.showXTicks   = true; end

    tickLen = [STY.TickLen STY.TickLen];

    nG = numel(R.groups);
    Yg = cell(1,nG);
    for g = 1:nG
        acc = 0;
        for c = R.groups{g}, acc = acc + D.Y{c}; end
        Yg{g} = acc / numel(R.groups{g});
    end
    Yg = within_subject(Yg);

    M = cell(1,nG); S = cell(1,nG);
    for g = 1:nG
        M{g} = mean(Yg{g},1,'omitnan');
        S{g} = std(Yg{g},0,1,'omitnan');
    end

    % CON blue and ECC green; 75% purple and 90% orange, matching the other
    % stages so a factor always looks the same across the figure set.
    cols = [0.00 0.4470 0.7410; 0.4660 0.6740 0.1880;
            0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560];
    if nG == 2 && contains(lower(R.effectLabel),'main effect: intensity')
        cols = cols([4 3],:);
    end

    % ---------- upper panel ----------
    ax1 = nexttile(tl,1); hold(ax1,'on');
    for g = 1:nG
        upG = min(M{g}+S{g}, 100);
        fill(ax1,[D.xAxis fliplr(D.xAxis)],[M{g}-S{g} fliplr(upG)], ...
            cols(g,:),'FaceAlpha',STY.FaceAlpha,'EdgeColor','none', ...
            'HandleVisibility','off');
    end
    h = gobjects(1,nG);
    for g = 1:nG
        h(g) = plot(ax1, D.xAxis, M{g}, 'Color', cols(g,:), 'LineWidth', STY.TraceWidth);
    end

    ylabel(ax1, sprintf('%s [%s]', D.varLabel, D.units), ...
        'FontWeight','bold','FontSize',STY.FontSize);
    if showTitles
        title(ax1, R.effectLabel, 'FontWeight','bold','FontSize',STY.FontSize);
    end
    legend(ax1, h, R.groupNames, 'Box','off','Location','southeast', ...
        'FontSize',STY.FontSize);
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
    end

    if ~isempty(opts.panelLetter)
        text(ax1, 0.01, 0.97, opts.panelLetter, 'Units','normalized', ...
            'FontSize',STY.FontSize+4, 'FontWeight','bold', ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
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
            'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
            'Clipping','off');
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
            yticks(ax2, 0:opts.spmYTickStep:opts.spmYLim(2));
        end
   else
        % Autoscaled panels still need a small negative floor, or an F trace
        % bounded at zero runs flat along the axis line.
        yl = ylim(ax2);
        ylim(ax2, [-0.03*yl(2), yl(2)]);
    end

    if ~opts.showXTicks
        set(ax2,'XTick',[]);
    end
end
