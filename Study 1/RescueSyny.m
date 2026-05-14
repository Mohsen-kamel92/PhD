%% ===== Stage 0 RESCUE (FIXED): Interpolate SPIKE Torque/Angle to EMG Timebase + TRUE common-timeline save =====
% Summary:
% 1) Load SPIKE (.mat) containing DAC1/Torque/Angle, crop to DAC1>0.85 window, re-zero SPIKE time.
% 2) Load MAECS (HD-EMG), crop by AUX>0 window, re-zero EMG time.
% 3) Estimate constant onset offset (DAC1 vs AUX), shift SPIKE time by -dt.
% 4) Interpolate SPIKE torque/angle/DAC1 onto EMG timebase (no extrapolation).
% 5) Crop all signals to common valid overlap (remove NaN edges).
% 6) Optionally estimate constant EMG lag relative to mechanics.
% 7) If rescue lag exists, ACTUALLY realign EMG and AUX onto the common timeline.
% 8) Build DAC1-based onset markers on the final saved timeline.
% 9) Plot QC figure using only the final common timeline.
% 10) Save synchronized dataset in the same practical format expected by the original downstream pipeline.

clear; clc; close all;

%% --- Optional rescue for problematic participant only ---
apply_emg_lag_correction = true;   % false for normal participants
manual_extra_emg_shift_s = 0.200;   % extra left shift of EMG/AUX in seconds
                                      % positive value = EMG is moved earlier
%% --- Load SPIKE .mat (expects fields: DAC1, Torque, Angle with .times/.values) ---
[SpikeFile, pathFile] = uigetfile('*.mat','Pick SPIKE mat file','MultiSelect','off');
if isequal(SpikeFile,0)
    return
end
S = load(fullfile(pathFile, SpikeFile));

mustHave = {'DAC1','Torque','Angle'};
for k = 1:numel(mustHave)
    if ~isfield(S, mustHave{k})
        error('Missing "%s" in SPIKE file.', mustHave{k});
    end
end

%% --- Crop SPIKE by DAC1 > 0.85 (first to last), re-zero time ---
th_dac = 0.85;
dac_values = S.DAC1.values(:);

i1 = find(dac_values > th_dac, 1, 'first');
i2 = find(dac_values > th_dac, 1, 'last');

if isempty(i1) || isempty(i2) || i2 <= i1
    error('Could not find a valid high segment in DAC1 using threshold %.2f.', th_dac);
end

SPIKE_vars = fieldnames(S);
for i = 1:numel(SPIKE_vars)
    fn = SPIKE_vars{i};
    if isstruct(S.(fn)) && isfield(S.(fn),'times') && isfield(S.(fn),'values')
        S.(fn).times  = S.(fn).times(i1:i2);
        S.(fn).values = S.(fn).values(i1:i2);
        S.(fn).times  = S.(fn).times - S.(fn).times(1);
    end
end

Fs_spike = 2000;   % nominal SPIKE sampling frequency
Fs_emg   = 2048;   % MAECS sampling frequency

%% --- Load MAECS (HD-EMG) and crop by AUX > 0, re-zero time ---
HDEMG = MAECS_read();

gNames = fieldnames(HDEMG);
assert(isscalar(gNames), 'Expected exactly one grid, found %d.', numel(gNames));
gridName = gNames{1};

if ~isfield(HDEMG.(gridName), 'AUX')
    error('AUX not found in HDEMG.%s', gridName);
end

if ~isfield(HDEMG.(gridName), 'Times')
    error('No Time vector found in HDEMG.%s.', gridName);
end

timeEMG = HDEMG.(gridName).Times(:);
aux     = HDEMG.(gridName).AUX(:);
emg     = HDEMG.(gridName).EMG;

assert(numel(aux) == numel(timeEMG), 'AUX length must match time length.');
assert(size(emg,1) == numel(timeEMG), 'EMG rows must match time length.');

aux(~isfinite(aux)) = 0;

j1 = find(aux > 0, 1, 'first');
j2 = find(aux > 0, 1, 'last');

if isempty(j1) || isempty(j2) || j2 <= j1
    error('Could not find a valid high segment in AUX (need AUX>0 span).');
end

emg_vals = emg(j1:j2, :);
emg_time = timeEMG(j1:j2);
emg_time = emg_time - emg_time(1);
AUX_cut  = aux(j1:j2);

%% --- Onset alignment (DAC1 vs AUX), interpolate SPIKE -> EMG timebase ---
tAux0_pre = emg_time(find(AUX_cut > 0, 1, 'first'));
tDAC0_pre = S.DAC1.times(find(S.DAC1.values > th_dac, 1, 'first'));

dt_align = tDAC0_pre - tAux0_pre;   % (+) DAC lags AUX

tTorque = S.Torque.times - dt_align;
tAngle  = S.Angle.times  - dt_align;
tDAC    = S.DAC1.times   - dt_align;

torque_i = interp1(tTorque, S.Torque.values, emg_time, 'linear');
angle_i  = interp1(tAngle,  S.Angle.values,  emg_time, 'linear');
dac1_i   = interp1(tDAC,    S.DAC1.values,   emg_time, 'linear');

%% --- Crop to common valid overlap (remove NaN edges) ---
valid = isfinite(torque_i) & isfinite(angle_i) & isfinite(dac1_i);

firstIdx = find(valid, 1, 'first');
lastIdx  = find(valid, 1, 'last');

if isempty(firstIdx) || isempty(lastIdx) || lastIdx <= firstIdx
    error('No valid overlap found after interpolation. Check dt_align and trigger windows.');
end

keep = firstIdx:lastIdx;

time0     = emg_time(keep);     % provisional common timeline
emg0      = emg_vals(keep,:);
AUX0      = AUX_cut(keep);
torque0   = torque_i(keep);
angle0    = angle_i(keep);
dac10     = dac1_i(keep);

%% ============================================================
% OPTIONAL EMG LAG CORRECTION (mechanics-based)
% IMPORTANT:
% We do NOT save separate EMG time.
% We estimate a constant lag, then REALIGN EMG/AUX onto the same final time axis.
% ============================================================
lag_const = NaN;
lag_total = NaN;
raw_aux_on_pre = [time0(1); time0(find(diff(AUX0 > 0) == 1) + 1)];
if isempty(raw_aux_on_pre)
    raw_aux_on_pre = time0(1);
end

if apply_emg_lag_correction

    fprintf('\n===== EMG LAG ESTIMATION =====\n');

    % Smooth angle for velocity estimation
    angle_s = Wfilt(angle0, 4, 'low', Fs_emg);

    % Angular velocity
    vel = [0; diff(angle_s)] * Fs_emg;

    % Smooth velocity
    vel = Wfilt(vel, 4, 'low', Fs_emg);

    % Threshold for movement detection
    vel_thr = 5;                 % deg/s
    minDur  = round(0.150*Fs_emg);

    vel_bin = vel > vel_thr;
    vel_bin = keepRuns(vel_bin, minDur);

    cv_idx = find(vel_bin);

    if isempty(cv_idx)
        warning('Could not detect CV onset. No EMG lag correction applied.');
    else
        t_cv_on = time0(cv_idx(1));

        % EMG envelope from provisional aligned/cropped EMG
        emg_env_tmp = sqrt(movmean(emg0.^2, round(0.050*Fs_emg), 1, 'Endpoints','shrink'));
        emg_env_mean_tmp = mean(emg_env_tmp, 2, 'omitnan');

        % EMG baseline
        base_idx = time0 < (t_cv_on - 0.5);
        if ~any(base_idx)
            base_idx = time0 < min(time0(1) + 0.5, time0(end));
        end
        baseline = emg_env_mean_tmp(base_idx);

        thr = mean(baseline, 'omitnan') + 3*std(baseline, 'omitnan');

        % EMG onset
        emg_bin = emg_env_mean_tmp > thr;
        emg_bin = keepRuns(emg_bin, round(0.050*Fs_emg));

        emg_idx = find(emg_bin);

        if isempty(emg_idx)
            warning('Could not detect EMG onset. No EMG lag correction applied.');
        else
            t_emg_on = time0(emg_idx(1));
            lag_const = t_emg_on - t_cv_on;
lag_total = lag_const + manual_extra_emg_shift_s;

fprintf('CV onset         : %.3f s\n', t_cv_on);
fprintf('EMG onset        : %.3f s\n', t_emg_on);
fprintf('Estimated lag    : %.3f s\n', lag_const);
fprintf('Manual extra lag : %.3f s\n', manual_extra_emg_shift_s);
fprintf('Total realign    : %.3f s\n', lag_total);
fprintf('EMG/AUX will be realigned onto common timeline by %.3f s\n', lag_total);
        end
    end

    fprintf('================================\n');
end

%% ============================================================
% REALIGN EMG and AUX ONTO THE FINAL COMMON TIMELINE
% ============================================================
% Original downstream scripts assume:
%   sync.time, sync.emg, sync.AUX, sync.torque_raw, sync.angle, sync.dac1
% are all aligned sample-by-sample on ONE timeline.
%
% Therefore, if lag_const exists, we resample EMG/AUX onto the same time axis
% rather than storing a second EMG time vector.

time     = time0;
torque_i = torque0;
angle_i  = angle0;
dac1_i   = dac10;

if apply_emg_lag_correction && isfinite(lag_total) && abs(lag_total) > 0
    emg_src_time = time0 - lag_total;
    aux_src_time = time0 - lag_total;
    % Interpolate EMG onto the final mechanics/DAC timeline
    emg_aligned = interp1(emg_src_time, emg0, time, 'linear', NaN);

    % Interpolate AUX onto the same timeline (nearest preserves trigger-like shape better)
    AUX_aligned = interp1(aux_src_time, double(AUX0), time, 'nearest', NaN);

    % Keep only overlap where all saved signals are valid
    valid2 = all(isfinite(emg_aligned), 2) & isfinite(AUX_aligned) & ...
             isfinite(torque_i) & isfinite(angle_i) & isfinite(dac1_i);

    firstIdx2 = find(valid2, 1, 'first');
    lastIdx2  = find(valid2, 1, 'last');

    if isempty(firstIdx2) || isempty(lastIdx2) || lastIdx2 <= firstIdx2
        error('No valid overlap remained after rescue EMG/AUX realignment.');
    end

    keep2 = firstIdx2:lastIdx2;

    time     = time(keep2);
    emg_vals = emg_aligned(keep2,:);
    AUX_cut  = AUX_aligned(keep2);
    torque_i = torque_i(keep2);
    angle_i  = angle_i(keep2);
    dac1_i   = dac1_i(keep2);

else
    emg_vals = emg0;
    AUX_cut  = AUX0;
    firstIdx2 = NaN;
    lastIdx2  = NaN;
end

% Clean AUX after interpolation
AUX_cut(~isfinite(AUX_cut)) = 0;
AUX_cut = double(AUX_cut > 0.5);   % restore clean digital-ish AUX

%% --- Compute onset markers on FINAL saved timeline ---
thr_dac_plot = th_dac;

% DAC-based markers = final reference markers for plotting + downstream
dac_on = [time(1); time(find(diff(dac1_i > thr_dac_plot) == 1) + 1)];
if isempty(dac_on)
    dac_on = time(1);
end

% Raw/final AUX detections on saved timeline (for metadata only)
raw_aux_on = [time(1); time(find(diff(AUX_cut > 0) == 1) + 1)];
if isempty(raw_aux_on)
    raw_aux_on = time(1);
end

% For downstream consistency with scripts expecting aux_on:
aux_on = dac_on;

% Offset between saved AUX and DAC on FINAL timeline (for reporting only)
tAux0_post = time(find(AUX_cut > 0,            1, 'first'));
tDAC0_post = time(find(dac1_i  > thr_dac_plot, 1, 'first'));

if isempty(tAux0_post), tAux0_post = time(1); end
if isempty(tDAC0_post), tDAC0_post = time(1); end

dt_ms = (tDAC0_post - tAux0_post) * 1e3;

fprintf('dt_align = %.3f ms; kept %.3f s overlap.\n', dt_align*1e3, time(end)-time(1));
if isfinite(lag_const)
    fprintf('Final saved AUX-DAC onset offset = %.3f ms\n', dt_ms);
end

%% --- QC plot (smoothed torque/angle for display only; EMG RMS mean) ---
% NOTE: Saved torque/angle remain unfiltered.
fc_angle  = 6;
fc_torque = 20;

angle_smooth  = Wfilt(angle_i,  fc_angle,  'low', Fs_emg);
torque_smooth = Wfilt(torque_i, fc_torque, 'low', Fs_emg);

win_sec = 0.050;
win = max(1, ceil(win_sec*Fs_emg));
emg_env = sqrt(movmean(emg_vals.^2, win, 1, 'Endpoints','shrink'));
emg_env_mean = mean(emg_env, 2, 'omitnan');

%% --- Plot tiled layout ---
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ================= Top panel =================
ax1 = nexttile(1);
plot(time, angle_smooth, 'Color',[0 0.55 0], 'LineWidth',1.8);
hold(ax1,'on');
ylabel('Angle (°)','FontWeight','bold');
set(ax1,'XTickLabel',[],'XColor','none');

xline(ax1, dac_on(1), '--k', 'LineWidth',1.2, 'HandleVisibility','on');
for i = 2:numel(dac_on)
    xline(ax1, dac_on(i), '--k', 'LineWidth',1.2, 'HandleVisibility','off');
end

% ================= Bottom panel =================
ax2 = nexttile(2);

yyaxis(ax2,'left')
plot(time, emg_env_mean, 'r', 'LineWidth',1.8);
hold(ax2,'on');
ylabel('EMG RMS mean (V)','FontWeight','bold');
ax2.YAxis(1).Color = [1 0 0];

yyaxis(ax2,'right')
plot(time, torque_smooth, 'Color',[0.5 0 0.5], 'LineWidth',1.8);
ylabel('Torque (Nm)','FontWeight','bold');
ax2.YAxis(2).Color = [0.5 0 0.5];
xlabel('Time (s)','FontWeight','bold');

yyaxis(ax2,'left')
for i = 1:numel(dac_on)
    xline(ax2, dac_on(i), '--k', 'LineWidth',1.2, 'HandleVisibility','off');
end

linkaxes([ax1,ax2],'x');
set(gcf,'Color','w');
set([ax1,ax2], 'Color','w', 'XGrid','off', 'YGrid','off', 'Box','off', ...
    'TickDir','out', 'FontWeight','bold', 'FontSize',12, 'LineWidth',1.2);

if apply_emg_lag_correction && isfinite(lag_total)
    sgtitle(sprintf('Trigger Alignment (saved AUX - DAC1 = %.2f ms) | rescue EMG/AUX realigned by %.3f s', dt_ms, lag_total), ...
        'FontWeight','bold');
else
    sgtitle(sprintf('Trigger Alignment (AUX - DAC1 = %.2f ms)', dt_ms), ...
        'FontWeight','bold');
end

%% --- Legend (DAC only) ---
hAngProxy = plot(ax2, NaN, NaN, '-',  'Color',[0 0.55 0],   'LineWidth',1.8);
hTQProxy  = plot(ax2, NaN, NaN, '-',  'Color',[0.5 0 0.5],  'LineWidth',1.8);
hEMGProxy = plot(ax2, NaN, NaN, '-',  'Color',[1 0 0],      'LineWidth',1.8);
hDACproxy = plot(ax2, NaN, NaN, '--', 'Color',[0 0 0],      'LineWidth',1.2);

lgd = legend(ax2, ...
    [hAngProxy hTQProxy hEMGProxy hDACproxy], ...
    {'Angle','Torque','EMG RMS-mean','DAC1'}, ...
    'Location','best', 'Box','off', 'FontSize',8);

lgd.ItemTokenSize = [10 6];

%% --- Folder structure assumed: .../Participant/CON_75/Set_1/ ---
[setPath, setFolder] = fileparts(pathFile(1:end-1));
[condPath, conditionFolder] = fileparts(setPath);
[~, participantFolder] = fileparts(condPath);

participantName = participantFolder;
conditionName   = conditionFolder;
setName         = setFolder;

%% --- Define method string ---
if apply_emg_lag_correction && isfinite(lag_total)
    method = 'onset-shift + interp-to-EMG + crop-overlap + rescue EMG/AUX realignment onto common timeline (estimated lag + manual extra shift)';
else
    method = 'onset-shift + interp-to-EMG + crop-overlap (no extrap)';
end

%% --- Save synchronized dataset as a single struct ---
% IMPORTANT:
% Save only the standard single-timeline contract expected by downstream scripts.
sync = struct();

sync.time     = time;
sync.Fs_emg   = Fs_emg;
sync.Fs_spike = Fs_spike;

sync.emg = emg_vals;     % NOW truly aligned to sync.time
sync.AUX = AUX_cut;      % NOW truly aligned to sync.time

sync.torque_raw = torque_i;
sync.angle      = angle_i;
sync.dac1       = dac1_i;

sync.aux_on = aux_on;    % intentionally DAC-based for downstream consistency
sync.dac_on = dac_on;

%% --- Metadata ---
sync.meta = struct();
sync.meta.method      = method;
sync.meta.thr_dac     = th_dac;
sync.meta.dt_align_s  = dt_align;
sync.meta.dt_align_ms = dt_align * 1e3;

sync.meta.participantName = participantName;
sync.meta.conditionName   = conditionName;
sync.meta.setName         = setName;

sync.meta.spike_file = SpikeFile;
sync.meta.spike_path = pathFile;

sync.meta.crop = struct();
sync.meta.crop.spike_i1         = i1;
sync.meta.crop.spike_i2         = i2;
sync.meta.crop.emg_j1           = j1;
sync.meta.crop.emg_j2           = j2;
sync.meta.crop.overlap_firstIdx = firstIdx;
sync.meta.crop.overlap_lastIdx  = lastIdx;
sync.meta.crop.realign_firstIdx = firstIdx2;
sync.meta.crop.realign_lastIdx  = lastIdx2;

% Rescue-specific tracking (metadata only)
sync.meta.reference_marker_source    = 'DAC1';
sync.meta.raw_aux_on                 = raw_aux_on;
sync.meta.saved_aux_on_equals_dac_on = true;

sync.meta.emg_lag_correction_applied = apply_emg_lag_correction && isfinite(lag_total);
sync.meta.emg_lag_estimated_s        = lag_const;
sync.meta.emg_lag_manual_extra_s     = manual_extra_emg_shift_s;
sync.meta.emg_lag_total_s            = lag_total;

outname = fullfile(pathFile, sprintf('%s_%s_sync_rawVolts.mat', participantName, conditionName));

try
    save(outname, 'sync', '-v7.3');
catch
    save(outname, 'sync');
end

fprintf('✅ Saved synchronized struct to: %s\n', outname);

%% ===== Helper =====
function bin_out = keepRuns(bin_in, minLen)
bin_in = logical(bin_in(:));
d = diff([0; bin_in; 0]);
runStarts = find(d == 1);
runEnds   = find(d == -1) - 1;

bin_out = false(size(bin_in));

for k = 1:numel(runStarts)
    runLen = runEnds(k) - runStarts(k) + 1;
    if runLen >= minLen
        bin_out(runStarts(k):runEnds(k)) = true;
    end
end
end