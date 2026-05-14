%% === Stage 1: Matching QC per Rep (raw torque vs raw ramp trace), saves repQC ===

clear; clc; close all;

%% --- Load RAMP (.txt) ---
[txtFile, txtPath] = uigetfile('*.txt', 'Select the Ramp .txt file');
if isequal(txtFile, 0)
    error('No Ramp file selected.');
end

ramp_data = readmatrix(fullfile(txtPath, txtFile), 'Delimiter', ',');
t_ramp = ramp_data(:,1);
y_ramp = ramp_data(:,2);

% Remove duplicates & zero-base
[t_ramp, idx_unique] = unique(t_ramp, 'stable');
y_ramp = y_ramp(idx_unique);
t_ramp = t_ramp - t_ramp(1);

ramp_duration = t_ramp(end);

%% --- Load synced struct (.mat) ---
[matFile, matPath] = uigetfile('*.mat', 'Select *_sync_rawVolts.mat');
if isequal(matFile, 0)
    return;
end

S = load(fullfile(matPath, matFile));
if isfield(S, 'sync')
    sync = S.sync;
else
    error('This file does not contain a "sync" struct. Re-save sync stage.');
end

time   = sync.time(:);
torque = sync.torque_raw(:);
angle  = sync.angle(:);

Fs = 2000;
condition = 'CON';   % <-- CHANGE to used condition e.g CON or ECC

%% --- Smooth for window detection only (QC) ---
fc_angle  = 6;
fc_torque = 20;

angle_smooth  = Wfilt(angle,  fc_angle,  'low', Fs);
torque_smooth = Wfilt(torque, fc_torque, 'low', Fs);

%% --- Onset detection (choose trigger onsets saved by sync) ---
if isfield(sync, 'dac_on') && ~isempty(sync.dac_on)
    onset_times = sync.dac_on(:);
elseif isfield(sync, 'aux_on') && ~isempty(sync.aux_on)
    onset_times = sync.aux_on(:);
else
    error('No onset markers found in sync (dac_on / aux_on).');
end

onset_idx = knnsearch(time, onset_times);
nReps = numel(onset_idx);

%% --- Active/CV window detection via angular velocity ---
% omega > 0 is ECC, omega < 0 is CON
angle_velocity = [0; diff(angle_smooth) ./ diff(time)];

velocity_thresh = 5;    % deg/s
min_duration    = 0.10; % sec

repIdx = nan(nReps, 2);   % [startIdx endIdx]
CVwin  = nan(nReps, 2);   % [tStart tEnd] in rep-relative time

for r = 1:nReps
    idx_start = onset_idx(r);

    % End rep using ramp duration from same onset
    t_end_target = time(idx_start) + ramp_duration;
    idx_end_by_ramp = find(time >= t_end_target, 1, 'first');

    if isempty(idx_end_by_ramp)
        idx_end_by_ramp = numel(time);
    end

    % Also prevent overlap with next rep
    if r < nReps
        idx_end = min(idx_end_by_ramp, onset_idx(r+1) - 1);
    else
        idx_end = idx_end_by_ramp;
    end

    repIdx(r,:) = [idx_start idx_end];

    t_rep = time(idx_start:idx_end);
    t_rel = t_rep - t_rep(1);
    v_rep = angle_velocity(idx_start:idx_end);

    switch upper(condition)
        case 'ECC'
            above = v_rep > velocity_thresh;
        case 'CON'
            above = v_rep < -velocity_thresh;
        otherwise
            error('Unknown condition. Use ''ECC'' or ''CON''.');
    end

    % Detect contiguous regions above threshold
    d = diff([0; above; 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    t_active_start = NaN;
    t_active_end   = NaN;

    for j = 1:numel(starts)
        t0 = t_rel(starts(j));
        t1 = t_rel(ends(j));

        if (t1 - t0) >= min_duration
            t_active_start = t0;
            t_active_end   = t1;
            break
        end
    end

    CVwin(r,:) = [t_active_start t_active_end];
end

%% --- Matching error per rep (raw vs raw within CV window) ---
matching_errors = nan(nReps,1);
n_valid_points  = nan(nReps,1);
is_valid_rep    = false(nReps,1);

validity_thresh_pct = 100;
eps_den = 1e-6;

repPlot = struct();
repPlot.t_rep    = cell(nReps,1);
repPlot.tq_rep   = cell(nReps,1);
repPlot.ramp_rep = cell(nReps,1);
repPlot.CVg      = nan(nReps,2);

for r = 1:nReps
    t_start = CVwin(r,1);
    t_end   = CVwin(r,2);

    if ~isfinite(t_start) || ~isfinite(t_end)
        fprintf('Rep %02d | skipped: invalid CV window\n', r);
        continue
    end

    idx_start = repIdx(r,1);
    idx_end   = repIdx(r,2);

    if ~isfinite(idx_start) || ~isfinite(idx_end) || idx_end <= idx_start
        fprintf('Rep %02d | skipped: invalid rep indices\n', r);
        continue
    end

    t_rep = time(idx_start:idx_end);
    t_rel = t_rep - t_rep(1);

    % Store traces for plotting
    repPlot.t_rep{r}    = t_rep;
    repPlot.tq_rep{r}   = torque_smooth(idx_start:idx_end);
    repPlot.ramp_rep{r} = interp1(t_ramp, y_ramp, t_rel, 'linear', NaN);

    % Convert CV window from rep-relative to global time
    repPlot.CVg(r,:) = t_rep(1) + CVwin(r,:);

    idx_rel = find(t_rel >= t_start & t_rel <= t_end);
    if numel(idx_rel) < 5
        fprintf('Rep %02d | skipped: too few samples in CV window\n', r);
        continue
    end

    idx = idx_start + idx_rel - 1;

    tq_active   = torque_smooth(idx);
    ramp_active = interp1(t_ramp, y_ramp, t_rel(idx_rel), 'linear', NaN);

    valid = isfinite(tq_active) & isfinite(ramp_active);
    if nnz(valid) < 5
        fprintf('Rep %02d | skipped: too few valid comparison samples\n', r);
        continue
    end

    tqv = tq_active(valid);
    rav = ramp_active(valid);

    den = max(abs(tqv) + abs(rav), eps_den);
    err_vec = abs(tqv - rav) ./ den * 100;

    matching_errors(r) = median(err_vec);
    n_valid_points(r)  = nnz(valid);
    is_valid_rep(r)    = matching_errors(r) <= validity_thresh_pct;

    fprintf('Rep %02d | err = %.2f%% | nPts = %d | valid = %d\n', ...
        r, matching_errors(r), n_valid_points(r), is_valid_rep(r));
end

%% --- Summary Output (ALL + VALID ONLY) ---
fprintf('\n=== Torque-to-Ramp Matching Error Summary (Median-based) ===\n');
fprintf('Condition: %s\n', upper(condition));
fprintf('Validity threshold: ≤ %.2f%%\n', validity_thresh_pct);

% All reps
fprintf('\n[All reps]\n');
fprintf('Mean of medians: %.2f%%\n', mean(matching_errors, 'omitnan'));
fprintf('Median of medians: %.2f%%\n', median(matching_errors, 'omitnan'));
fprintf('Min–Max: %.2f%% to %.2f%%\n', ...
    min(matching_errors, [], 'omitnan'), max(matching_errors, [], 'omitnan'));
fprintf('Std (of medians): %.2f%%\n', std(matching_errors, 'omitnan'));

% Valid reps only
validMask = is_valid_rep & isfinite(matching_errors);
nValid = nnz(validMask);

fprintf('\n[Valid reps only]\n');
fprintf('Valid reps: %d / %d (%.1f%%)\n', nValid, nReps, 100*nValid/nReps);

if nValid > 0
    fprintf('Mean of medians (valid): %.2f%%\n', mean(matching_errors(validMask), 'omitnan'));
    fprintf('Median of medians (valid): %.2f%%\n', median(matching_errors(validMask), 'omitnan'));
    fprintf('Min–Max (valid): %.2f%% to %.2f%%\n', ...
        min(matching_errors(validMask), [], 'omitnan'), ...
        max(matching_errors(validMask), [], 'omitnan'));
    fprintf('Std (valid): %.2f%%\n', std(matching_errors(validMask), 'omitnan'));
else
    fprintf('No valid reps under threshold.\n');
end

%% --- Plotting (Angle + Torque vs Ramp per rep, highlight invalid reps) ---
figure('Color','w', 'Name', 'Matching QC - Torque vs Ramp (per rep)');
tl = tiledlayout(2,1);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

% --- Upper: Angle
ax1 = nexttile(tl, 1);
hold(ax1, 'on');
plot(ax1, time, angle_smooth, 'g-', 'LineWidth', 1.3);
ylabel(ax1, 'Angle (°)', 'FontWeight', 'bold');
set(ax1, 'Box', 'off', 'XColor', 'none', 'TickDir', 'out', ...
    'LineWidth', 1.2, 'FontSize', 11, 'FontWeight', 'bold');

% --- Lower: Torque + Ramp
ax2 = nexttile(tl, 2);
hold(ax2, 'on');

% Lock y-limits using torque trace
plot(ax2, time, torque_smooth, 'k', 'LineWidth', 0.1, 'HandleVisibility', 'off');
y_lim_fixed = ylim(ax2);

for r = 1:nReps
    t_rep    = repPlot.t_rep{r};
    tq_rep   = repPlot.tq_rep{r};
    ramp_rep = repPlot.ramp_rep{r};

    if isempty(t_rep) || isempty(tq_rep) || isempty(ramp_rep)
        continue
    end

    valid = isfinite(tq_rep) & isfinite(ramp_rep);
    if nnz(valid) < 2
        continue
    end

    t_valid    = t_rep(valid);
    tq_valid   = tq_rep(valid);
    ramp_valid = ramp_rep(valid);

    if is_valid_rep(r)
        tqColor = [0.5 0 0.5];
        lw = 1.3;
    else
        tqColor = [0.65 0.65 0.65];
        lw = 1.0;
    end

    plot(ax2, t_valid, ramp_valid, 'b-', 'LineWidth', 1.3);
    plot(ax2, t_valid, tq_valid, '-', 'Color', tqColor, 'LineWidth', lw);

    % CV window in global time
    t_start = repPlot.CVg(r,1);
    t_end   = repPlot.CVg(r,2);

    if isfinite(t_start) && isfinite(t_end)
        mask_active = (t_valid >= t_start) & (t_valid <= t_end);

        if nnz(mask_active) > 2
            x_fill = [t_valid(mask_active); flipud(t_valid(mask_active))];
            y_fill = [ramp_valid(mask_active); flipud(tq_valid(mask_active))];
            fill(ax2, x_fill, y_fill, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        end

        plot(ax2, [t_start t_start], y_lim_fixed, 'k--', 'HandleVisibility', 'off');
        plot(ax2, [t_end   t_end],   y_lim_fixed, 'k--', 'HandleVisibility', 'off');
    end

    if isfinite(matching_errors(r))
        txt = sprintf('R%02d: %.1f%%', r, matching_errors(r));
        text(ax2, t_valid(1), y_lim_fixed(2) * 0.95, txt, ...
            'FontSize', 7, 'FontWeight', 'bold', 'Color', tqColor, ...
            'VerticalAlignment', 'top');
    end
end

xlabel(ax2, 'Time (s)', 'FontWeight', 'bold');
ylabel(ax2, 'Torque (Nm)', 'FontWeight', 'bold');
set(ax2, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 11, ...
    'FontWeight', 'bold', 'Box', 'off');

mean_err_all = mean(matching_errors, 'omitnan');
title(ax1, sprintf('Torque-Ramp Matching Error: %.2f%% SPDs | Threshold ≤ %.0f%%', ...
    mean_err_all, validity_thresh_pct), ...
    'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');

% Legend (proxy handles)
hAng     = plot(ax2, NaN, NaN, 'g-', 'LineWidth', 1.3);
hRamp    = plot(ax2, NaN, NaN, 'b-', 'LineWidth', 1.3);
hTQvalid = plot(ax2, NaN, NaN, '-', 'Color', [0.5 0 0.5], 'LineWidth', 1.3);
hTQinv   = plot(ax2, NaN, NaN, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0);
hFill    = fill(ax2, NaN, NaN, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
hWin     = plot(ax2, NaN, NaN, 'k--', 'LineWidth', 1.2);

legend(ax2, [hAng hRamp hTQvalid hTQinv hFill hWin], ...
    {'Angle', 'Ramp', 'Torque (valid)', 'Torque (invalid)', 'Error area (CV)', 'CV window'}, ...
    'Location', 'northeast', 'FontSize', 7, 'Box', 'off');

%% --- Package outputs into a repQC struct and save ---
repQC = struct();
repQC.onset_times      = onset_times;
repQC.onset_idx        = onset_idx;
repQC.repIdx           = repIdx;
repQC.CVwin_time       = CVwin;
repQC.CVwin_globalTime = repPlot.CVg;
repQC.matchErr_pct_med = matching_errors;
repQC.n_valid_points   = n_valid_points;
repQC.is_valid_rep     = is_valid_rep;
repQC.valid_rep_idx    = find(is_valid_rep);

repQC.params = struct();
repQC.params.condition            = condition;
repQC.params.velocity_thresh      = velocity_thresh;
repQC.params.min_duration         = min_duration;
repQC.params.validity_thresh_pct  = validity_thresh_pct;
repQC.params.eps_den              = eps_den;
repQC.params.fc_angle             = fc_angle;
repQC.params.fc_torque            = fc_torque;
repQC.params.ramp_duration        = ramp_duration;

% Save alongside sync
[pth, name, ~] = fileparts(fullfile(matPath, matFile));
outname = fullfile(pth, sprintf('%s_stage1_matchQC.mat', name));
save(outname, 'sync', 'repQC', '-v7.3');

fprintf('\n✅ Saved Stage-1 match QC: %s\n', outname);
fprintf('Valid reps: %d / %d (threshold %.2f%%)\n', nnz(is_valid_rep), nReps, validity_thresh_pct);