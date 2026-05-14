%% ===== Stage 6C1: MID 20-80% Time-Series Inference with SPM1D (RM ANOVA) =====
% PURPOSE:
%   Perform 2-way repeated-measures SPM ANOVA on the MID 20-80% time series for:
%       - torque
%       - mean EMG
%       - centroid X
%       - centroid Y
%
% FACTORS:
%   A) Contraction type: CON vs ECC
%   B) Intensity: 75% vs 90%
%
% INPUT:
%   groupData_stage6A.mat
%
% OUTPUT:
%   stage6C1_spmResults_mid20to80.mat
%
% IMPORTANT:
%   - Uses mainBin = 'all'
%   - Restricts analysis to 20-80% CV
%   - If flipECCforAnalysis = true, ECC is reversed left-right before both
%     plotting and inference.

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
alphaLevel = 0.05;
interpClusters = true;
flipECCforAnalysis = true;
doPlots = true;
saveFigures = false;

midWin = [20 80];

condNames = string(groupData.conditions(:));
expectedOrder = ["CON_75"; "ECC_75"; "CON_90"; "ECC_90"];
assert(all(condNames == expectedOrder), ...
    'Condition order mismatch. Expected [CON_75, ECC_75, CON_90, ECC_90].');

nSub = groupData.nSub;
xAxis_full = groupData.xAxis(:).';

midMask = xAxis_full >= midWin(1) & xAxis_full <= midWin(2);
assert(any(midMask), 'No samples found inside the requested mid window.');

xAxis = xAxis_full(midMask);
nT = numel(xAxis);

fprintf('\n===== Stage 6C1 MID loaded =====\n');
fprintf('Participants: %d\n', nSub);
fprintf('Main bin: %s\n', mainBin);
fprintf('Alpha: %.3f\n', alphaLevel);
fprintf('Flip ECC for analysis: %d\n', flipECCforAnalysis);
fprintf('Analysis window: %d-%d %%CV\n', midWin(1), midWin(2));

%% ---------------------------
%  (3) Check spm1d availability
%% ---------------------------
spmPath = which('spm1d.stats.anova2rm');

assert(~isempty(spmPath), ...
    ['spm1d toolbox not found on MATLAB path. ', ...
     'Add spm1d first, e.g. addpath(genpath(''path_to_spm1d''));']);

fprintf('SPM1D found at:\n%s\n', spmPath);

%% ---------------------------
%  (4) Extract MID 20-80% time series
%% ---------------------------
T_torque = groupData.timeseries.(mainBin).torque_101(:,:,midMask);
T_emg    = groupData.timeseries.(mainBin).emgMeanNorm_101(:,:,midMask);
T_cx     = groupData.timeseries.(mainBin).centroidXc_101(:,:,midMask);
T_cy     = groupData.timeseries.(mainBin).centroidYc_101(:,:,midMask);

assert(all(size(T_torque) == [nSub 4 nT]));
assert(all(size(T_emg)    == [nSub 4 nT]));
assert(all(size(T_cx)     == [nSub 4 nT]));
assert(all(size(T_cy)     == [nSub 4 nT]));

%% ---------------------------
%  (5) Optionally flip ECC conditions
%% ---------------------------
if flipECCforAnalysis
    T_torque(:,2,:) = flip(T_torque(:,2,:), 3);
    T_torque(:,4,:) = flip(T_torque(:,4,:), 3);

    T_emg(:,2,:) = flip(T_emg(:,2,:), 3);
    T_emg(:,4,:) = flip(T_emg(:,4,:), 3);

    T_cx(:,2,:) = flip(T_cx(:,2,:), 3);
    T_cx(:,4,:) = flip(T_cx(:,4,:), 3);

    T_cy(:,2,:) = flip(T_cy(:,2,:), 3);
    T_cy(:,4,:) = flip(T_cy(:,4,:), 3);
end

%% ---------------------------
%  (6) Build output struct
%% ---------------------------
spmResults = struct();
spmResults.meta = struct();
spmResults.meta.createdOn = datestr(now);
spmResults.meta.script = mfilename;
spmResults.meta.sourceFile = fullfile(p,f);
spmResults.meta.mainBin = mainBin;
spmResults.meta.nSub = nSub;
spmResults.meta.nTime = nT;
spmResults.meta.xAxis = xAxis;
spmResults.meta.alpha = alphaLevel;
spmResults.meta.interpClusters = interpClusters;
spmResults.meta.flipECCforAnalysis = flipECCforAnalysis;
spmResults.meta.analysisWindow = midWin;
spmResults.meta.design = '2x2 repeated-measures ANOVA';
spmResults.meta.factorA = 'Contraction';
spmResults.meta.factorALevels = {'CON','ECC'};
spmResults.meta.factorB = 'Intensity';
spmResults.meta.factorBLevels = {'75','90'};

%% ---------------------------
%  (7) Run SPM RM ANOVA
%% ---------------------------
spmResults.torque = run_spm_anova2rm( ...
    T_torque, xAxis, alphaLevel, interpClusters, ...
    'Torque', 'Nm');

spmResults.emg = run_spm_anova2rm( ...
    T_emg, xAxis, alphaLevel, interpClusters, ...
    'Mean EMG', '%MVC');

spmResults.cx = run_spm_anova2rm( ...
    T_cx, xAxis, alphaLevel, interpClusters, ...
    'Centroid X', 'mm');

spmResults.cy = run_spm_anova2rm( ...
    T_cy, xAxis, alphaLevel, interpClusters, ...
    'Centroid Y', 'mm');

%% ---------------------------
%  (8) Print summary
%% ---------------------------
fprintf('\n===== Stage 6C1 MID 20-80%% SPM RM-ANOVA SUMMARY =====\n');

print_spm_summary(spmResults.torque.effects.contraction);
print_spm_summary(spmResults.torque.effects.intensity);
print_spm_summary(spmResults.torque.effects.interaction);

print_spm_summary(spmResults.emg.effects.contraction);
print_spm_summary(spmResults.emg.effects.intensity);
print_spm_summary(spmResults.emg.effects.interaction);

print_spm_summary(spmResults.cx.effects.contraction);
print_spm_summary(spmResults.cx.effects.intensity);
print_spm_summary(spmResults.cx.effects.interaction);

print_spm_summary(spmResults.cy.effects.contraction);
print_spm_summary(spmResults.cy.effects.intensity);
print_spm_summary(spmResults.cy.effects.interaction);

%% ---------------------------
%  (9) Figures
%% ---------------------------
if doPlots
    figs = struct();

    figs.torque_contraction = plot_spm_result(spmResults.torque.effects.contraction, spmResults.torque.descriptives);
    figs.torque_intensity   = plot_spm_result(spmResults.torque.effects.intensity,   spmResults.torque.descriptives);
    figs.torque_interaction = plot_spm_result(spmResults.torque.effects.interaction, spmResults.torque.descriptives);

    figs.emg_contraction = plot_spm_result(spmResults.emg.effects.contraction, spmResults.emg.descriptives);
    figs.emg_intensity   = plot_spm_result(spmResults.emg.effects.intensity,   spmResults.emg.descriptives);
    figs.emg_interaction = plot_spm_result(spmResults.emg.effects.interaction, spmResults.emg.descriptives);

    figs.cx_contraction = plot_spm_result(spmResults.cx.effects.contraction, spmResults.cx.descriptives);
    figs.cx_intensity   = plot_spm_result(spmResults.cx.effects.intensity,   spmResults.cx.descriptives);
    figs.cx_interaction = plot_spm_result(spmResults.cx.effects.interaction, spmResults.cx.descriptives);

    figs.cy_contraction = plot_spm_result(spmResults.cy.effects.contraction, spmResults.cy.descriptives);
    figs.cy_intensity   = plot_spm_result(spmResults.cy.effects.intensity,   spmResults.cy.descriptives);
    figs.cy_interaction = plot_spm_result(spmResults.cy.effects.interaction, spmResults.cy.descriptives);

    if saveFigures
        outFigDir = uigetdir(pwd, 'Select folder to save Stage6C1 figures');
        if ~isequal(outFigDir,0)
            exportgraphics(figs.torque_contraction, fullfile(outFigDir, 'Stage6C1_MID_SPM_torque_mainContraction.png'), 'Resolution', 300);
            exportgraphics(figs.torque_intensity,   fullfile(outFigDir, 'Stage6C1_MID_SPM_torque_mainIntensity.png'),   'Resolution', 300);
            exportgraphics(figs.torque_interaction, fullfile(outFigDir, 'Stage6C1_MID_SPM_torque_interaction.png'),     'Resolution', 300);

            exportgraphics(figs.emg_contraction, fullfile(outFigDir, 'Stage6C1_MID_SPM_emg_mainContraction.png'), 'Resolution', 300);
            exportgraphics(figs.emg_intensity,   fullfile(outFigDir, 'Stage6C1_MID_SPM_emg_mainIntensity.png'),   'Resolution', 300);
            exportgraphics(figs.emg_interaction, fullfile(outFigDir, 'Stage6C1_MID_SPM_emg_interaction.png'),     'Resolution', 300);

            exportgraphics(figs.cx_contraction, fullfile(outFigDir, 'Stage6C1_MID_SPM_centroidX_mainContraction.png'), 'Resolution', 300);
            exportgraphics(figs.cx_intensity,   fullfile(outFigDir, 'Stage6C1_MID_SPM_centroidX_mainIntensity.png'),   'Resolution', 300);
            exportgraphics(figs.cx_interaction, fullfile(outFigDir, 'Stage6C1_MID_SPM_centroidX_interaction.png'),     'Resolution', 300);

            exportgraphics(figs.cy_contraction, fullfile(outFigDir, 'Stage6C1_MID_SPM_centroidY_mainContraction.png'), 'Resolution', 300);
            exportgraphics(figs.cy_intensity,   fullfile(outFigDir, 'Stage6C1_MID_SPM_centroidY_mainIntensity.png'),   'Resolution', 300);
            exportgraphics(figs.cy_interaction, fullfile(outFigDir, 'Stage6C1_MID_SPM_centroidY_interaction.png'),     'Resolution', 300);
        end
    end
end

%% ---------------------------
%  (10) Save results
%% ---------------------------
[outFile, outPath] = uiputfile('stage6C1_spmResults_mid20to80.mat', 'Save Stage6C1 SPM results as');
if isequal(outFile,0)
    outPath = pwd;
    outFile = 'stage6C1_spmResults_mid20to80.mat';
end

save(fullfile(outPath, outFile), 'spmResults', '-v7.3');
fprintf('\n✅ Saved Stage 6C1 MID SPM results:\n%s\n', fullfile(outPath, outFile));

%% =========================================================
%  Local functions
%% =========================================================
function R = run_spm_anova2rm(T, xAxis, alphaLevel, interpClusters, varLabel, units)
% T is [nSub x 4 x nT]

    [nSub, nCond, nT] = size(T);
    assert(nCond == 4, 'T must have 4 conditions in dimension 2.');

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
                case 1
                    A(row) = 1; B(row) = 1;
                case 2
                    A(row) = 2; B(row) = 1;
                case 3
                    A(row) = 1; B(row) = 2;
                case 4
                    A(row) = 2; B(row) = 2;
            end
            SUBJ(row) = s;
        end
    end

    spmList = spm1d.stats.anova2rm(Y, A, B, SUBJ);

    spmA  = spmList('A');
    spmB  = spmList('B');
    spmAB = spmList('AB');

    spmiA  = local_infer_spm(spmA,  alphaLevel, interpClusters);
    spmiB  = local_infer_spm(spmB,  alphaLevel, interpClusters);
    spmiAB = local_infer_spm(spmAB, alphaLevel, interpClusters);

    D = struct();
    D.xAxis = xAxis;
    D.varLabel = varLabel;
    D.units = units;
    D.condNames = {'CON_75','ECC_75','CON_90','ECC_90'};
    D.Y = cell(1,4);
    D.mean = cell(1,4);
    D.sd = cell(1,4);

    for c = 1:4
        Yc = squeeze(T(:,c,:));
        D.Y{c} = Yc;
        D.mean{c} = mean(Yc,1,'omitnan');
        D.sd{c} = std(Yc,0,1,'omitnan');
    end

    R = struct();
    R.varLabel = varLabel;
    R.units = units;
    R.xAxis = xAxis;
    R.descriptives = D;

    R.effects = struct();
    R.effects.contraction = build_effect_struct( ...
        spmA, spmiA, xAxis, varLabel, units, ...
        'Main effect: Contraction', 'Contraction');

    R.effects.intensity = build_effect_struct( ...
        spmB, spmiB, xAxis, varLabel, units, ...
        'Main effect: Intensity', 'Intensity');

    R.effects.interaction = build_effect_struct( ...
        spmAB, spmiAB, xAxis, varLabel, units, ...
        'Interaction: Contraction × Intensity', 'Interaction');
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
        nodeAxis = 0:(nT-1);   % SPM node coordinates for this cropped domain

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

            % --- Preferred: endpoints from cluster object
            try
                xy = cl.endpoints;

                % Convert from node space -> xAxis space if needed
                if all(isfinite(xy))
                    startPct = interp1(nodeAxis, xAxis, xy(1), 'linear', 'extrap');
                    endPct   = interp1(nodeAxis, xAxis, xy(2), 'linear', 'extrap');
                    extentPct = endPct - startPct;
                end
            catch
                % fallback
                try
                    inds = cl.indices;

                    % indices are node indices; convert to xAxis
                    inds = inds(:);
                    inds = max(min(round(inds), nT-1), 0);   % keep within 0..nT-1
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

            C(k).startPct  = startPct;
            C(k).endPct    = endPct;
            C(k).extentPct = extentPct;
            C(k).P         = pval;
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

function fig = plot_spm_result(R, D)
    fig = figure('Color','w', 'Name', sprintf('SPM | %s | %s', R.varLabel, R.effectLabel));
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl,1); hold(ax1,'on');

    cols = [...
        0.00 0.4470 0.7410;
        0.4660 0.6740 0.1880;
        0.8500 0.3250 0.0980;
        0.4940 0.1840 0.5560];

    h = gobjects(1,4);
    for c = 1:4
        fill(ax1, [D.xAxis fliplr(D.xAxis)], ...
            [D.mean{c}-D.sd{c} fliplr(D.mean{c}+D.sd{c})], ...
            cols(c,:), 'FaceAlpha',0.10, 'EdgeColor','none', ...
            'HandleVisibility','off');

        h(c) = plot(ax1, D.xAxis, D.mean{c}, 'Color', cols(c,:), 'LineWidth',2.0);
    end

    if contains(lower(R.varLabel), 'centroid')
        yline(ax1, 0, '--k', 'LineWidth',0.9, 'HandleVisibility','off');
    end

    ylabel(ax1, sprintf('%s (%s)', R.varLabel, R.units), 'FontWeight','bold');
    title(ax1, sprintf('%s | %s', R.varLabel, R.effectLabel), ...
        'Interpreter','none', 'FontWeight','bold');

    legend(ax1, h, strrep(D.condNames,'_',' '), 'Box','off', 'Location','best');
    box(ax1,'off');
    set(ax1,'TickDir','out', 'LineWidth',1.2, 'FontWeight','bold');
    xlim(ax1, [D.xAxis(1) D.xAxis(end)]);

    ax2 = nexttile(tl,2); hold(ax2,'on');

    if ~isempty(R.z)
        hSPM = plot(ax2, R.xAxis, R.z, 'k', 'LineWidth',1.8, ...
            'DisplayName','SPM{F}');
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