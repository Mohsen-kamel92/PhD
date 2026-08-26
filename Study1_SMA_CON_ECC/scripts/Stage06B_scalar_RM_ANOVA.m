%% ===== Stage 06B: Scalar 2x2 RM-ANOVA over the analysis window =====
%
% PURPOSE
%   Test the four scalar outcomes, averaged over the analysis window, in a two
%   by two repeated measures design, and draw the main text figures.
%
% SUPPORTS
%   The torque matching, EMG amplitude, and centroid in
%   the Results section of the manuscript.
%
% INPUT
%   <ROOT>/groupData_stage6A.mat                       (Stage 06A)
%
% OUTPUT
%   <ROOT>/stage6B_stats_mid20to80.mat
%   <ROOT>/Figures_Stage6B/*.png
%   <ROOT>/stage6B_results.csv
%
% DESIGN
%   2 x 2 repeated measures: contraction type (CON, ECC) by intensity (75, 90).
%   Condition order is fixed: 1 = CON_75, 2 = ECC_75, 3 = CON_90, 4 = ECC_90.
%
% ON SPHERICITY
%   Each factor has two levels, so only one difference score exists per factor
%   and sphericity is satisfied by definition. Greenhouse-Geisser and
%   Huynh-Feldt corrections are not applicable, and uncorrected p values are
%   reported.
%
% ON MULTIPLE COMPARISONS
%   Eight paired follow-up tests are run. Holm-corrected p values are reported
%   alongside the uncorrected ones.
%
% ON THE REVERSED ECC TRACES
%   ECC series are reversed for display so that both modes are shown at matched
%   joint angles. The scalar means are unaffected by the reversal.
%
% ON TORQUE
%   Torque is a manipulation check rather than an outcome. It is tested in the
%   same model so the degree of matching can be judged, not so that a
%   difference in it can be interpreted.
%
% ON THE SHADED BANDS IN THE TIME-SERIES FIGURES
%   Centroid position differs by up to 11 mm between participants, because the
%   grid was positioned relative to each participant's own muscle-tendon
%   junction.That between-subject variance is irrelevant to a within-subject
%   comparison but would dominate a conventional SD band and hide the effect being tested.
%   For the centroid panels the bands are therefore computed after removing each participant's own mean across the
%   two conditions shown, so they reflect within-subject variability. Torque and
%   EMG panels keep conventional between-subject SDs, since their absolute
%   values are directly interpretable.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped.
%
% DEPENDENCIES
%   Statistics and Machine Learning Toolbox   fitrm, ranova, ttest
%   swtest.m                  Shapiro-Wilk         (MATLAB File Exchange)
%   distinguishable_colors.m  participant colours  (MATLAB File Exchange)

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'mainBin',         'all');
CFG = set_default(CFG, 'midWin',          [20 80]);   % analysis window, % of CV
CFG = set_default(CFG, 'alphaLevel',      0.05);
CFG = set_default(CFG, 'doWorkAnalysis',  true);      % section (6)
CFG = set_default(CFG, 'doPlots',         true);
CFG = set_default(CFG, 'saveFigures',     true);
CFG = set_default(CFG, 'flipECCforPlots', true);      % display only
CFG = set_default(CFG, 'showTitles',      true);      % off for the manuscript figures

mainBin    = CFG.mainBin;
midWin     = CFG.midWin;
alphaLevel = CFG.alphaLevel;

%% ---------------------------
%  (2) Data root and group data
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 6B: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

for dep = {'swtest','distinguishable_colors','fitrm'}
    assert(~isempty(which(dep{1})), ...
        'Stage 6B: %s not found. See the dependencies in the header.', dep{1});
end

srcFile = fullfile(ROOT, 'groupData_stage6A.mat');
assert(isfile(srcFile), ...
    'Stage 6B: %s not found. Run Stage 6A first.', srcFile);

L = load(srcFile);
assert(isfield(L,'groupData'), 'File does not contain "groupData".');
groupData = L.groupData;

condNames     = string(groupData.conditions(:));
expectedOrder = ["CON_75"; "ECC_75"; "CON_90"; "ECC_90"];
assert(all(condNames == expectedOrder), ...
    'Condition order mismatch. Expected [CON_75, ECC_75, CON_90, ECC_90].');

nSub  = groupData.nSub;
nCond = groupData.nCond;
nT    = groupData.nTime;
xAxis = groupData.xAxis(:).';

figDir = fullfile(ROOT, 'Figures_Stage6B');
if CFG.saveFigures && ~exist(figDir,'dir'), mkdir(figDir); end

fprintf('\n===== Stage 6B: scalar statistics =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Participants: %d | bin: %s | window: %d-%d %%CV\n', ...
    nSub, mainBin, midWin(1), midWin(2));

%% ---------------------------
%  (3) Stats structure
%% ---------------------------
stats6B = struct();
stats6B.meta = struct( ...
    'createdOn', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'script', mfilename, 'sourceFile', srcFile, ...
    'mainBin', mainBin, 'nSub', nSub, ...
    'conditions', condNames, 'alpha', alphaLevel, ...
    'midWin_percentCV', midWin, ...
    'sphericity_note', ["two levels per factor, so sphericity holds by " ...
                        "definition and no correction is applicable"]);

stats6B.design = struct('factor1_name','contractionType', ...
    'factor1_levels',{{'CON','ECC'}}, 'factor2_name','intensity', ...
    'factor2_levels',{{'75','90'}});

%% ---------------------------
%  (4) Time series and window scalars
%% ---------------------------
TS = groupData.timeseries.(mainBin);
for k = {'torque_101','emgMeanNorm_101','centroidXc_101','centroidYc_101','angle_101'}
    assert(isfield(TS,k{1}), 'groupData.timeseries.%s.%s missing.', mainBin, k{1});
end

T_torque = TS.torque_101;
T_emg    = TS.emgMeanNorm_101;
T_cx     = TS.centroidXc_101;
T_cy     = TS.centroidYc_101;
T_angle  = TS.angle_101;

for V = {T_torque,T_emg,T_cx,T_cy,T_angle}
    assert(isequal(size(V{1}),[nSub nCond nT]), 'Unexpected time-series size.');
end

midMask = xAxis >= midWin(1) & xAxis <= midWin(2);
assert(any(midMask), 'No samples inside the requested window.');

Y_torque = mean(T_torque(:,:,midMask), 3, 'omitnan');
Y_emg    = mean(T_emg(:,:,midMask),    3, 'omitnan');
Y_cx     = mean(T_cx(:,:,midMask),     3, 'omitnan');
Y_cy     = mean(T_cy(:,:,midMask),     3, 'omitnan');

stats6B.mid20to80 = struct('mask', midMask, 'xAxis', xAxis(midMask), ...
    'Y_torque', Y_torque, 'Y_emg', Y_emg, 'Y_cx', Y_cx, 'Y_cy', Y_cy);

% Angle range per contraction mode. Averaging across modes is meaningless here,
% because CON runs from plantarflexion to dorsiflexion and ECC runs the other
% way, so the two cancel to a flat mid-range value.
angCON = squeeze(mean(mean(T_angle(:,[1 3],:),1,'omitnan'),2,'omitnan')).';
angECC = squeeze(mean(mean(T_angle(:,[2 4],:),1,'omitnan'),2,'omitnan')).';
stats6B.angleAxis = struct('CON', angCON, 'ECC', angECC);
fprintf('Angle range: CON %.1f to %.1f deg | ECC %.1f to %.1f deg\n', ...
    angCON(1), angCON(end), angECC(1), angECC(end));

%% ---------------------------
%  (5) Torque matching check
%% ---------------------------
fprintf('\n===== TORQUE MATCHING (manipulation check) =====\n');

d75 = Y_torque(:,2) - Y_torque(:,1);          % ECC minus CON
d90 = Y_torque(:,4) - Y_torque(:,3);
rel75 = (d75 ./ Y_torque(:,1)) * 100;
rel90 = (d90 ./ Y_torque(:,3)) * 100;

mm = struct();
for I = ["75","90"]
    if I == "75", d = d75; rel = rel75; else, d = d90; rel = rel90; end
    [~, pv, ci, st] = ttest(d);
    dz = mean(d,'omitnan') / std(d,0,'omitnan');

    fprintf('%s%%: ECC-CON = %+.3f +/- %.3f Nm (%+.2f +/- %.2f %%), ', ...
        I, mean(d,'omitnan'), std(d,0,'omitnan'), ...
        mean(rel,'omitnan'), std(rel,0,'omitnan'));
    fprintf('t(%d) = %.2f, p = %.4f, dz = %.2f, 95%% CI [%.3f %.3f]\n', ...
        st.df, st.tstat, pv, dz, ci(1), ci(2));

    mm.("d"+I)   = d;    mm.("rel"+I) = rel;
    mm.("p"+I)   = pv;   mm.("t"+I)   = st.tstat;
    mm.("dz"+I)  = dz;   mm.("ci"+I)  = ci(:).';
end
stats6B.mid20to80.mismatch = mm;

%% ---------------------------
%  (6) Mechanical work
%% ---------------------------
% Work is torque integrated over joint angle. CON and ECC traverse the range of
% motion in opposite directions, so the integral carries opposite signs and the
% magnitudes are compared. Angle is converted to radians so that work is in
% joules.
%
% Work is computed over the whole constant velocity phase and over the analysis
% window, and the angular excursion is reported alongside, because a difference
% in excursion between modes would contribute to any difference in work.

if CFG.doWorkAnalysis
    fprintf('\n===== MECHANICAL WORK (torque integrated over joint angle) =====\n');

    W_full = nan(nSub, nCond);
    W_mid  = nan(nSub, nCond);
    ROMdeg = nan(nSub, nCond);

    for s = 1:nSub
        for c = 1:nCond
            th = squeeze(T_angle(s,c,:)) * pi/180;
            tq = squeeze(T_torque(s,c,:));
            okv = isfinite(th) & isfinite(tq);

            if nnz(okv) >= 5
                W_full(s,c)  = abs(trapz(th(okv), tq(okv)));
                ROMdeg(s,c)  = abs(th(find(okv,1,'last')) - th(find(okv,1,'first'))) * 180/pi;
            end

            mv = okv(:) & midMask(:);
            if nnz(mv) >= 5
                W_mid(s,c) = abs(trapz(th(mv), tq(mv)));
            end
        end
    end

    stats6B.work = struct('W_full_J', W_full, 'W_mid_J', W_mid, 'ROM_deg', ROMdeg, ...
        'note', ["magnitude of torque integrated over joint angle; CON and ECC " ...
                 "integrate in opposite directions"]);

    wk = struct();
    for W = ["full","mid"]
        if W == "full", Wm = W_full; else, Wm = W_mid; end

        fprintf('\n--- %s CVP ---\n', upper(W));
        fprintf('%-8s %14s %14s %10s %10s %10s\n', ...
            'level','CON (J)','ECC (J)','diff (J)','diff (%)','p');

        for I = ["75","90"]
            if I == "75", a = 1; b = 2; else, a = 3; b = 4; end
            d   = Wm(:,b) - Wm(:,a);                 % ECC minus CON
            rel = 100 * d ./ Wm(:,a);
            [~, pv, ci, st] = ttest(d);
            dz = mean(d,'omitnan') / std(d,0,'omitnan');

            fprintf('%-8s %6.2f +/-%5.2f %6.2f +/-%5.2f %+10.3f %+10.2f %10.4f\n', ...
                I+"%", mean(Wm(:,a),'omitnan'), std(Wm(:,a),0,'omitnan'), ...
                mean(Wm(:,b),'omitnan'), std(Wm(:,b),0,'omitnan'), ...
                mean(d,'omitnan'), mean(rel,'omitnan'), pv);

            wk.(W).("i"+I) = struct( ...
                'meanCON', mean(Wm(:,a),'omitnan'), 'sdCON', std(Wm(:,a),0,'omitnan'), ...
                'meanECC', mean(Wm(:,b),'omitnan'), 'sdECC', std(Wm(:,b),0,'omitnan'), ...
                'diff', mean(d,'omitnan'), 'sdDiff', std(d,0,'omitnan'), ...
                'relPct', mean(rel,'omitnan'), 'sdRelPct', std(rel,0,'omitnan'), ...
                't', st.tstat, 'df', st.df, 'p', pv, 'dz', dz, 'ci', ci(:).');
        end
    end
    stats6B.work.tests = wk;

    % How much of the mismatch the window restriction removes
    fprintf('\n--- mismatch reduction from restricting to %d-%d%% CVP ---\n', ...
        midWin(1), midWin(2));
    for I = ["75","90"]
        if I == "75", a = 1; b = 2; else, a = 3; b = 4; end
        rFull = mean(abs(100*(W_full(:,b)-W_full(:,a))./W_full(:,a)), 'omitnan');
        rMid  = mean(abs(100*(W_mid(:,b) -W_mid(:,a)) ./W_mid(:,a)),  'omitnan');
        fprintf('%s%% MVC: %.2f%% over the full CVP, %.2f%% over the window\n', ...
            I, rFull, rMid);
    end

    % A difference in excursion would contribute to the work difference
    fprintf('\nROM: CON %.1f +/- %.1f deg | ECC %.1f +/- %.1f deg\n', ...
        mean(ROMdeg(:,[1 3]),'all','omitnan'), std(ROMdeg(:,[1 3]),0,'all','omitnan'), ...
        mean(ROMdeg(:,[2 4]),'all','omitnan'), std(ROMdeg(:,[2 4]),0,'all','omitnan'));
end

%% ---------------------------
%  (7) Normality of the paired differences
%% ---------------------------
fprintf('\n===== NORMALITY (Shapiro-Wilk on paired differences) =====\n');

outcomes = struct('Torque',Y_torque,'EMG',Y_emg,'CentroidX',Y_cx,'CentroidY',Y_cy);
fn = fieldnames(outcomes);
stats6B.normality = struct();

for i = 1:numel(fn)
    Y = outcomes.(fn{i});
    [~, p75, W75] = swtest(Y(:,1)-Y(:,2), alphaLevel);
    [~, p90, W90] = swtest(Y(:,3)-Y(:,4), alphaLevel);

    fprintf('%-10s 75%%: W = %.3f, p = %.4f%s\n', fn{i}, W75, p75, flag(p75, alphaLevel));
    fprintf('%-10s 90%%: W = %.3f, p = %.4f%s\n', fn{i}, W90, p90, flag(p90, alphaLevel));

    stats6B.normality.(fn{i}) = struct('W75',W75,'p75',p75,'W90',W90,'p90',p90);
end

%% ---------------------------
%  (8) RM-ANOVA
%% ---------------------------
stats6B.mid20to80.scalar.meanTorque    = run_rm2x2(Y_torque, groupData.participants, 'meanTorque');
stats6B.mid20to80.scalar.meanEMG       = run_rm2x2(Y_emg,    groupData.participants, 'meanEMG');
stats6B.mid20to80.scalar.meanCentroidX = run_rm2x2(Y_cx,     groupData.participants, 'meanCentroidX');
stats6B.mid20to80.scalar.meanCentroidY = run_rm2x2(Y_cy,     groupData.participants, 'meanCentroidY');

fprintf('\n===== RM-ANOVA (2 x 2, %d-%d%% CV window) =====\n', midWin(1), midWin(2));
print_rm_result(stats6B.mid20to80.scalar.meanTorque,    'Mean torque');
print_rm_result(stats6B.mid20to80.scalar.meanEMG,       'Mean EMG');
print_rm_result(stats6B.mid20to80.scalar.meanCentroidX, 'Mean centroid X');
print_rm_result(stats6B.mid20to80.scalar.meanCentroidY, 'Mean centroid Y');

%% ---------------------------
%  (9) Paired follow-up tests, Holm corrected
%% ---------------------------
ph = struct();
ph.meanTorque.CONvsECC_75    = paired_compare(Y_torque(:,1), Y_torque(:,2), 'Torque CON75 vs ECC75');
ph.meanTorque.CONvsECC_90    = paired_compare(Y_torque(:,3), Y_torque(:,4), 'Torque CON90 vs ECC90');
ph.meanEMG.CONvsECC_75       = paired_compare(Y_emg(:,1),    Y_emg(:,2),    'EMG CON75 vs ECC75');
ph.meanEMG.CONvsECC_90       = paired_compare(Y_emg(:,3),    Y_emg(:,4),    'EMG CON90 vs ECC90');
ph.meanCentroidX.CONvsECC_75 = paired_compare(Y_cx(:,1),     Y_cx(:,2),     'Centroid X CON75 vs ECC75');
ph.meanCentroidX.CONvsECC_90 = paired_compare(Y_cx(:,3),     Y_cx(:,4),     'Centroid X CON90 vs ECC90');
ph.meanCentroidY.CONvsECC_75 = paired_compare(Y_cy(:,1),     Y_cy(:,2),     'Centroid Y CON75 vs ECC75');
ph.meanCentroidY.CONvsECC_90 = paired_compare(Y_cy(:,3),     Y_cy(:,4),     'Centroid Y CON90 vs ECC90');

% Holm correction is applied across all eight tests together, so the family is
% the whole set of follow-ups rather than each outcome separately.
outer = fieldnames(ph);
praw = []; ref = {};
for a = 1:numel(outer)
    inner = fieldnames(ph.(outer{a}));
    for b = 1:numel(inner)
        praw(end+1) = ph.(outer{a}).(inner{b}).p;  %#ok<SAGROW>
        ref{end+1}  = {outer{a}, inner{b}};        %#ok<SAGROW>
    end
end

pHolm = holm(praw);
for i = 1:numel(ref)
    ph.(ref{i}{1}).(ref{i}{2}).p_holm = pHolm(i);
end
stats6B.mid20to80.posthoc = ph;

fprintf('\n===== PAIRED COMPARISONS (Holm corrected across %d tests) =====\n', numel(praw));
fprintf('%-28s %9s %9s %8s %8s %20s\n','comparison','mean diff','t','p','p_holm','95% CI');
for i = 1:numel(ref)
    H = ph.(ref{i}{1}).(ref{i}{2});
    fprintf('%-28s %9.3f %9.2f %8.4f %8.4f  [%7.3f %7.3f]  dz = %.2f\n', ...
        H.label, H.meanDiff, H.stats.tstat, H.p, H.p_holm, H.ci(1), H.ci(2), H.dz);
end

%% ---------------------------
%  (10) Results table
%% ---------------------------
res = table('Size',[numel(ref) 9], ...
    'VariableTypes',{'string','double','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'comparison','mean1','mean2','meanDiff','t','df','p','p_holm','dz'});
for i = 1:numel(ref)
    H = ph.(ref{i}{1}).(ref{i}{2});
    res(i,:) = {string(H.label), H.mean1, H.mean2, H.meanDiff, ...
                H.stats.tstat, H.stats.df, H.p, H.p_holm, H.dz};
end
writetable(res, fullfile(ROOT,'stage6B_results.csv'));

%% ---------------------------
%  (11) Figures
%% ---------------------------
if CFG.doPlots
    figs = struct();
    showTitles = CFG.showTitles;
    doFlip     = CFG.flipECCforPlots;

    % --- scalar panels, one per outcome ---
    figs.scalar_torque = scalar_figure(Y_torque, 'Torque [Nm]',                  'Mean torque', showTitles);
    figs.scalar_emg    = scalar_figure(Y_emg,    'RMS EMG [% MVA]',              'Mean EMG',    showTitles);
    figs.scalar_cx     = scalar_figure(Y_cx,     'Centroid X displacement [mm]', 'Centroid X',  showTitles);
    figs.scalar_cy     = scalar_figure(Y_cy,     'Centroid Y displacement [mm]', 'Centroid Y',  showTitles);

    % Centroid panels share a fixed range so the two axes can be compared
    for fh = {figs.scalar_cx, figs.scalar_cy}
        axh = findobj(fh{1}, 'Type', 'axes');
        ylim(axh, [-7 8]);
        yticks(axh, -6:2:6);
    end

    % --- torque and EMG, both intensities in one figure ---
    figs.torque_emg = figure('Color','w','Name','Stage6B | torque and EMG');
    tl = tiledlayout(figs.torque_emg,2,2,'TileSpacing','compact','Padding','compact');

    ax = gobjects(2,2);
    ax(1,1) = nexttile(tl,1);
    [hC, hE] = plot_condition_pair(ax(1,1), xAxis, squeeze(T_torque(:,1,:)), ...
        maybe_flip(squeeze(T_torque(:,2,:)), doFlip), ...
        'Torque [Nm]', ttl(showTitles,'75% MVC'), false, false, midWin, false);

    ax(1,2) = nexttile(tl,2);
    plot_condition_pair(ax(1,2), xAxis, squeeze(T_torque(:,3,:)), ...
        maybe_flip(squeeze(T_torque(:,4,:)), doFlip), ...
        '', ttl(showTitles,'90% MVC'), false, false, midWin, false);

    ax(2,1) = nexttile(tl,3);
    plot_condition_pair(ax(2,1), xAxis, squeeze(T_emg(:,1,:)), ...
        maybe_flip(squeeze(T_emg(:,2,:)), doFlip), ...
        'RMS EMG [% MVA]','', false, true, midWin, false);

    ax(2,2) = nexttile(tl,4);
    plot_condition_pair(ax(2,2), xAxis, squeeze(T_emg(:,3,:)), ...
        maybe_flip(squeeze(T_emg(:,4,:)), doFlip), ...
        '','', false, true, midWin, false);

    link_ylim(ax(1,:));      % torque row
    link_ylim(ax(2,:));      % EMG row

    ylim(ax(1,1), [0 40]); ylim(ax(1,2), [0 40]);
    yticks(ax(1,1), 0:10:40); yticks(ax(1,2), 0:10:40);
    set(ax(1,2), 'YTick', []);
    set(ax(2,2), 'YTick', []);

    legend(ax(1,2), [hC hE], {'CON','ECC'}, 'Box','off', ...
        'Location','northeast', 'FontSize',12);

    add_panel_letters([ax(1,1) ax(1,2) ax(2,1) ax(2,2)]);

    % --- centroid, both coordinates and both intensities ---
    figs.centroid = figure('Color','w','Name','Stage6B | centroid time series');
    tl = tiledlayout(figs.centroid,2,2,'TileSpacing','compact','Padding','compact');

    axc = gobjects(2,2);
    axc(1,1) = nexttile(tl,1);
    [hCx, hEx] = plot_condition_pair(axc(1,1), xAxis, squeeze(T_cx(:,1,:)), ...
        maybe_flip(squeeze(T_cx(:,2,:)), doFlip), ...
        'Centroid X [mm]', ttl(showTitles,'75% MVC'), true, false, midWin, true);

    axc(1,2) = nexttile(tl,2);
    plot_condition_pair(axc(1,2), xAxis, squeeze(T_cx(:,3,:)), ...
        maybe_flip(squeeze(T_cx(:,4,:)), doFlip), ...
        '', ttl(showTitles,'90% MVC'), true, false, midWin, true);

    axc(2,1) = nexttile(tl,3);
    plot_condition_pair(axc(2,1), xAxis, squeeze(T_cy(:,1,:)), ...
        maybe_flip(squeeze(T_cy(:,2,:)), doFlip), ...
        'Centroid Y [mm]','', true, true, midWin, true);

    axc(2,2) = nexttile(tl,4);
    plot_condition_pair(axc(2,2), xAxis, squeeze(T_cy(:,3,:)), ...
        maybe_flip(squeeze(T_cy(:,4,:)), doFlip), ...
        '','', true, true, midWin, true);

    legend(axc(1,1), [hCx hEx], {'CON','ECC'}, 'Box','off', ...
        'Location','northeast', 'FontSize',12);

    for a = axc(:)'
        ylim(a, [-6 6]); yticks(a, -6:2:6);
    end
    set(axc(1,2), 'YTick', []);
    set(axc(2,2), 'YTick', []);

    add_panel_letters([axc(1,1) axc(1,2) axc(2,1) axc(2,2)]);

    if CFG.saveFigures
        names = fieldnames(figs);
        for i = 1:numel(names)
            exportgraphics(figs.(names{i}), ...
                fullfile(figDir, sprintf('Stage6B_%s.png', names{i})), 'Resolution', 600);
        end
        fprintf('\nFigures saved to %s\n', figDir);
    end
end

%% ---------------------------
%  (12) Save
%% ---------------------------
outName = fullfile(ROOT, 'stage6B_stats_mid20to80.mat');
save(outName, 'stats6B', '-v7.3');
fprintf('Saved stats  : %s\n', outName);
fprintf('Saved results: %s\n', fullfile(ROOT,'stage6B_results.csv'));

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
% selection fails immediately with a useful message.
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

function s = flag(pv, alphaLevel)
    if pv < alphaLevel, s = '   *** normality rejected'; else, s = ''; end
end

function t = ttl(show, txt)
    if show, t = txt; else, t = ''; end
end

function link_ylim(axes_)
    lims = cell2mat(get(axes_,'YLim'));
    set(axes_,'YLim',[min(lims(:,1)) max(lims(:,2))]);
end

function add_panel_letters(axList)
% Letters sit inside each panel, so they render correctly regardless of where
% the panel falls in the layout.
    letters = {'A','B','C','D','E','F'};
    for i = 1:min(numel(axList), numel(letters))
        text(axList(i), 0.01, 0.97, letters{i}, 'Units','normalized', ...
            'FontWeight','bold', 'FontSize',16, ...
            'HorizontalAlignment','left', 'VerticalAlignment','top');
    end
end

function pAdj = holm(p)
    [ps, idx] = sort(p(:));
    m = numel(ps);
    adj = min(cummax(ps .* (m - (1:m)' + 1)), 1);
    pAdj = nan(size(ps));
    pAdj(idx) = adj;
    pAdj = reshape(pAdj, size(p));
end

function S = run_rm2x2(Y, participants, varLabel)
    T = table(string(participants(:)), Y(:,1), Y(:,2), Y(:,3), Y(:,4), ...
        'VariableNames', {'Subject','CON_75','ECC_75','CON_90','ECC_90'});

    within = table(categorical({'CON';'ECC';'CON';'ECC'}), ...
                   categorical({'75';'75';'90';'90'}), ...
                   'VariableNames', {'Contraction','Intensity'});

    rm  = fitrm(T, 'CON_75-ECC_90 ~ 1', 'WithinDesign', within);
    tbl = ranova(rm, 'WithinModel', 'Contraction*Intensity');

    S = struct('label', varLabel, 'data', Y, 'subjects', string(participants(:)), ...
        'table', T, 'within', within, 'rm', rm, 'ranova', tbl);

    S.p_contraction = get_p(tbl, "(Intercept):Contraction");
    S.p_intensity   = get_p(tbl, "(Intercept):Intensity");
    S.p_interaction = get_p(tbl, "(Intercept):Contraction:Intensity");

    S.F_contraction = get_f(tbl, "(Intercept):Contraction");
    S.F_intensity   = get_f(tbl, "(Intercept):Intensity");
    S.F_interaction = get_f(tbl, "(Intercept):Contraction:Intensity");

    S.eta2p_contraction = eta2p(tbl, "(Intercept):Contraction",  "Error(Contraction)");
    S.eta2p_intensity   = eta2p(tbl, "(Intercept):Intensity",    "Error(Intensity)");
    S.eta2p_interaction = eta2p(tbl, "(Intercept):Contraction:Intensity", ...
                                     "Error(Contraction:Intensity)");

    S.df_effect = 1;
    S.df_error  = size(Y,1) - 1;
    S.mean = mean(Y,1,'omitnan');
    S.sd   = std(Y,0,1,'omitnan');
end

function v = get_p(tbl, rowName)
    v = NaN;
    i = find(string(tbl.Properties.RowNames) == rowName, 1);
    if isempty(i), return; end
    if any(strcmp(tbl.Properties.VariableNames,'pValue')) && isfinite(tbl.pValue(i))
        v = tbl.pValue(i);
    end
end

function v = get_f(tbl, rowName)
    v = NaN;
    i = find(string(tbl.Properties.RowNames) == rowName, 1);
    if ~isempty(i) && any(strcmp(tbl.Properties.VariableNames,'F'))
        v = tbl.F(i);
    end
end

function v = eta2p(tbl, effRow, errRow)
% Partial eta squared, from the sums of squares of the effect and its error.
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
    fprintf('  means [CON75 ECC75 CON90 ECC90]: %s\n', num2str(S.mean, '%8.2f'));
    fprintf('  SDs                            : %s\n', num2str(S.sd,   '%8.2f'));
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
        'sd1', std(y1,0,'omitnan'),  'sd2', std(y2,0,'omitnan'), ...
        'meanDiff', mean(d,'omitnan'), 'sdDiff', std(d,0,'omitnan'), ...
        'dz', mean(d,'omitnan')/std(d,0,'omitnan'), ...
        'ci', ci(:).', 'p', pv, 'p_holm', NaN, 'stats', st);
end

function fh = scalar_figure(Y, ylab, titleText, showTitles)
    fh = figure('Color','w','Name',titleText);
    ax = axes(fh);
    plot_scalar_points_panel(ax, Y, ttl(showTitles, titleText), ylab);
end

function Yout = maybe_flip(Yin, doFlip)
% ECC is reversed for display so both modes are shown at matched joint angles.
    Yout = Yin;
    if doFlip, Yout = fliplr(Yin); end
end

function plot_scalar_points_panel(ax, Y, ttlText, ylab)
% Every participant shown individually, connected within each intensity, with
% the group mean and SD overlaid.
    hold(ax,'on');
    xPos = [1 2 4 5];

    yMin0 = min(Y(:),[],'omitnan');
    yMax0 = max(Y(:),[],'omitnan');
    yRange0 = yMax0 - yMin0;
    if yRange0 == 0, yRange0 = 1; end
    pad = 0.12 * yRange0;
    yLo = yMin0 - pad;
    yHi = yMax0 + pad;

    nSub = size(Y,1);
    subjColors = distinguishable_colors(nSub, {'w','k'});
    subjColors = subjColors(1:nSub,:);

    lineColor    = [0.75 0.75 0.75];
    lineWidthInd = 0.5;
    markerSize   = 100;

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
                    'MarkerFaceColor', thisColor, 'MarkerEdgeColor','k', 'LineWidth',0.4);
            end
        end
    end

    m = mean(Y,1,'omitnan');
    s = std(Y,0,1,'omitnan');

    capW = 0.25;
    for c = 1:4
        if ~isnan(m(c)) && ~isnan(s(c))
            plot(ax, [xPos(c) xPos(c)], [m(c)-s(c), m(c)+s(c)], ...
                '-', 'Color',[0.35 0.35 0.35], 'LineWidth',2.0);
            plot(ax, [xPos(c)-capW xPos(c)+capW], [m(c)-s(c) m(c)-s(c)], ...
                '-', 'Color',[0.35 0.35 0.35], 'LineWidth',2.0);
            plot(ax, [xPos(c)-capW xPos(c)+capW], [m(c)+s(c) m(c)+s(c)], ...
                '-', 'Color',[0.35 0.35 0.35], 'LineWidth',2.0);
        end
    end

    plot(ax, xPos(1:2), m(1:2), '--', 'Color','k','LineWidth',1.2);
    plot(ax, xPos(3:4), m(3:4), '--', 'Color','k','LineWidth',1.2);
    scatter(ax, xPos, m, 42, 'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',0.5);

    xlim(ax, [0.5 5.5]);
    ylim(ax, [yLo yHi]);

    set(ax, 'XTick', xPos, 'XTickLabel', {'CON75','ECC75','CON90','ECC90'}, ...
        'FontSize',12, 'FontWeight','bold', 'LineWidth',2.5, 'TickDir','out');

    ylabel(ax, ylab, 'FontWeight','bold','FontSize',12);
    if ~isempty(ttlText)
        title(ax, ttlText, 'Interpreter','none','FontWeight','bold','FontSize',12);
    end
    box(ax,'off');

    if yLo < 0 && yHi > 0
        yline(ax, 0, ':k', 'LineWidth',2);
    end
end

function [hCON, hECC] = plot_condition_pair(ax, xAxis, Y1, Y2, ylab, ttlText, ...
                             addZeroLine, showXlabel, winCV, normWithin)
% CON against ECC for one outcome at one intensity, with the analysis window
% marked. See the header for why the centroid bands are within-subject.
    if nargin < 7  || isempty(addZeroLine), addZeroLine = false; end
    if nargin < 8  || isempty(showXlabel),  showXlabel  = true;  end
    if nargin < 9,  winCV = []; end
    if nargin < 10 || isempty(normWithin), normWithin = false; end

    hold(ax,'on');

    if normWithin
        % Remove each participant's own offset across the two conditions, so the
        % band reflects within-subject variability rather than differences in
        % where the grid sat on each leg. The grand mean is added back so the
        % plotted values stay on the original scale.
        G  = (mean(Y1,2,'omitnan') + mean(Y2,2,'omitnan')) / 2;
        gm = mean(G,'omitnan');
        Y1p = Y1 - G + gm;
        Y2p = Y2 - G + gm;
    else
        Y1p = Y1;  Y2p = Y2;
    end

    m1 = mean(Y1p,1,'omitnan'); s1 = std(Y1p,0,1,'omitnan');
    m2 = mean(Y2p,1,'omitnan'); s2 = std(Y2p,0,1,'omitnan');

    colCON = [0.00 0.4470 0.7410];
    colECC = [0.4660 0.6740 0.1880];

    fill(ax, [xAxis fliplr(xAxis)], [m1-s1 fliplr(m1+s1)], colCON, ...
        'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
    fill(ax, [xAxis fliplr(xAxis)], [m2-s2 fliplr(m2+s2)], colECC, ...
        'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');

    hCON = plot(ax, xAxis, m1, 'Color', colCON, 'LineWidth', 2.2);
    hECC = plot(ax, xAxis, m2, 'Color', colECC, 'LineWidth', 2.2);

    if ~isempty(winCV) && numel(winCV) == 2
        xline(ax, winCV(1), '--k', 'LineWidth',2.0, 'HandleVisibility','off');
        xline(ax, winCV(2), '--k', 'LineWidth',2.0, 'HandleVisibility','off');
    end

    if addZeroLine, yline(ax, 0, ':k', 'LineWidth',2); end

    if ~isempty(ylab), ylabel(ax, ylab, 'FontWeight','bold','FontSize',12); end
    if ~isempty(ttlText)
        title(ax, ttlText, 'Interpreter','none', 'FontWeight','bold','FontSize',12);
    end

    box(ax,'off');
    set(ax, 'TickDir','out', 'LineWidth',2.5, 'FontWeight','bold', 'FontSize',12);
    xlim(ax,[0 100]); xticks(ax,0:20:100);

    if showXlabel
        xlabel(ax,'Constant velocity phase [%]','FontWeight','bold','FontSize',12);
    else
        set(ax,'XTick',[]);
    end
end
