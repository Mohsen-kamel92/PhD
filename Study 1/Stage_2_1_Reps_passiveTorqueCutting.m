%% ===== Stage 2.1: TRIAL PASSIVE TORQUE CORRECTION + ACTIVE WINDOW DETECTION =====
% INPUT : *_stage1_matchQC.mat   (contains sync + repQC)
% OUTPUT: *_ActiveTorqueCorrected.mat
%
% PURPOSE:
%   - Load Stage-1 output (sync + repQC)
%   - Filter torque and angle
%   - Estimate passive torque from full trial
%   - Subtract passive torque to obtain active torque
%   - Detect active window inside each rep using angular velocity
%   - Save repQC.CVwin_idx as numeric [nReps x 2] for Stage 3
%
% LOGIC:
%   - Same style as MVC Stage 2.2
%   - No fractions
%   - No complex helper functions
%   - No longest-segment logic
%   - First valid velocity-based segment inside each rep is used

clear; clc; close all;

%% ---------------- User options ----------------
CutOff_torque = 20;
CutOff_angle  = 6;
doQCplot      = true;

% ---- Active window detection ----
velocity_thresh = 5;      % deg/s
min_duration    = 0.10;   % s

%% ---------------- Load Stage-1 file ----------------
[file, path] = uigetfile('*_stage1_matchQC.mat', ...
    'Select Stage-1 matchQC file (*_stage1_matchQC.mat)');
if isequal(file,0), return; end

L = load(fullfile(path, file));

if ~isfield(L,'sync')
    error('Selected file does not contain "sync" struct.');
end
if ~isfield(L,'repQC')
    error('Selected file does not contain "repQC" struct.');
end

sync  = L.sync;
repQC = L.repQC;

%% ---------------- Core signals ----------------
time       = sync.time(:);
torque_raw = sync.torque_raw(:);
angle_raw  = sync.angle(:);

if isfield(sync,'emg')
    emg_raw = sync.emg;
else
    error('sync.emg is missing. Stage 3 requires emg_raw.');
end

if isfield(sync,'AUX')
    AUX = sync.AUX(:);
else
    AUX = [];
end

if isfield(sync,'dac1')
    dac1 = sync.dac1(:);
elseif isfield(sync,'DAC1')
    dac1 = sync.DAC1(:);
else
    dac1 = [];
end

Fs = 1 / median(diff(time));

%% ---------------- Metadata ----------------
if isfield(sync,'meta') && isfield(sync.meta,'participantName')
    participantName = string(sync.meta.participantName);
else
    participantName = "Unknown";
end

if isfield(sync,'meta') && isfield(sync.meta,'conditionName')
    conditionName = string(sync.meta.conditionName);
else
    conditionName = "Trial";
end

%% ---------------- Validate repQC ----------------
assert(isfield(repQC,'repIdx') && ~isempty(repQC.repIdx), ...
    'repQC.repIdx missing/empty.');

repIdx = repQC.repIdx;
nReps  = size(repIdx,1);

if size(repIdx,2) ~= 2
    error('repQC.repIdx must be [nReps x 2].');
end

if isfield(repQC,'is_valid_rep') && ~isempty(repQC.is_valid_rep)
    is_valid_rep = logical(repQC.is_valid_rep(:));
else
    is_valid_rep = true(nReps,1);
    repQC.is_valid_rep = is_valid_rep;
end

if isfield(repQC,'valid_rep_idx') && ~isempty(repQC.valid_rep_idx)
    valid_rep_idx = repQC.valid_rep_idx(:);
else
    valid_rep_idx = find(is_valid_rep);
    repQC.valid_rep_idx = valid_rep_idx;
end

%% ---------------- Filtering ----------------
torque_filt = real(Wfilt(torque_raw, CutOff_torque, 'low', Fs));
angle_filt  = real(Wfilt(angle_raw,  CutOff_angle,  'low', Fs));

%% ---------------- Passive subtraction ----------------
passive_p = Auto_RestTorque(torque_filt, angle_filt, Fs);
[tau_pass_hat, tau_active] = Normalize_Torque(torque_filt, angle_filt, passive_p);

% Consistent naming
torque_active = tau_active;
angle         = angle_filt;

%% ---------------- Angular velocity ----------------
angle_velocity = [0; diff(angle) ./ diff(time)];

%% ---------------- Detect active window inside each rep ----------------
activeWin_idx        = nan(nReps,2);   % global indices
activeWin_globalTime = nan(nReps,2);   % global time
activeWin_time       = nan(nReps,2);   % rep-relative time

% For Stage 3, use same active window directly as CV window
CVwin_idx            = nan(nReps,2);
CVwin_globalTime     = nan(nReps,2);
CVwin_time           = nan(nReps,2);

for r = 1:nReps
    idx0_rep = repIdx(r,1);
    idx1_rep = repIdx(r,2);

    if ~isfinite(idx0_rep) || ~isfinite(idx1_rep) || idx1_rep <= idx0_rep || ...
            idx0_rep < 1 || idx1_rep > numel(time)
        warning('Rep %d has invalid repIdx. Window left as NaN.', r);
        continue
    end

    idx_rep = idx0_rep:idx1_rep;
    t_rep   = time(idx_rep);
    v_rep   = angle_velocity(idx_rep);

    above = abs(v_rep) > velocity_thresh;
    d = diff([0; above; 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    rep_found = false;

    for j = 1:numel(starts)
        t0 = t_rep(starts(j));
        t1 = t_rep(ends(j));

        if (t1 - t0) >= min_duration
            ig0 = idx_rep(starts(j));
            ig1 = idx_rep(ends(j));

            activeWin_idx(r,:)        = [ig0 ig1];
            activeWin_globalTime(r,:) = [time(ig0) time(ig1)];
            activeWin_time(r,:)       = [time(ig0)-time(idx0_rep), time(ig1)-time(idx0_rep)];

            % same window used as CV window
            CVwin_idx(r,:)            = [ig0 ig1];
            CVwin_globalTime(r,:)     = [time(ig0) time(ig1)];
            CVwin_time(r,:)           = [time(ig0)-time(idx0_rep), time(ig1)-time(idx0_rep)];

            rep_found = true;
            break
        end
    end

    if ~rep_found
        warning('Rep %d: no valid velocity-based active window detected.', r);
    end
end

%% ---------------- Save back into repQC ----------------
repQC.activeWin_idx        = activeWin_idx;
repQC.activeWin_globalTime = activeWin_globalTime;
repQC.activeWin_time       = activeWin_time;

repQC.CVwin_idx            = CVwin_idx;         % Stage 3 expects this
repQC.CVwin_globalTime     = CVwin_globalTime;
repQC.CVwin_time           = CVwin_time;

repQC.detector.method           = 'velocity-based active window inside repIdx';
repQC.detector.velocity_thresh  = velocity_thresh;
repQC.detector.min_duration     = min_duration;

%% ---------------- Optional matched rep-wise storage ----------------
time_matched         = cell(nReps,1);
angle_matched        = cell(nReps,1);
torque_matched       = cell(nReps,1);
emg_raw_matched      = cell(nReps,1);

torque_raw_matched   = cell(nReps,1);
torque_filt_matched  = cell(nReps,1);
angle_raw_matched    = cell(nReps,1);
angle_filt_matched   = cell(nReps,1);

for r = 1:nReps
    idx0_rep = repIdx(r,1);
    idx1_rep = repIdx(r,2);

    if ~isfinite(idx0_rep) || ~isfinite(idx1_rep) || idx1_rep <= idx0_rep || ...
            idx0_rep < 1 || idx1_rep > numel(time)
        continue
    end

    idx = idx0_rep:idx1_rep;

    time_matched{r}        = time(idx);
    angle_matched{r}       = angle(idx);
    torque_matched{r}      = torque_active(idx);

    torque_raw_matched{r}  = torque_raw(idx);
    torque_filt_matched{r} = torque_filt(idx);
    angle_raw_matched{r}   = angle_raw(idx);
    angle_filt_matched{r}  = angle(idx);

    emg_raw_matched{r}     = emg_raw(idx,:);
end

%% ---------------- Flags ----------------
angle_is_filtered  = true;
torque_is_filtered = true;
angle_filter_Hz    = CutOff_angle;
torque_filter_Hz   = CutOff_torque;
filter_method      = "Wfilt low-pass";

method = ['Trial passive fit (Auto_RestTorque + Normalize_Torque) ' ...
          '+ active window detected by angular velocity inside each rep'];

sourceFile = fullfile(path, file);

%% ---------------- QC plot ----------------
if doQCplot
    figure('Color','w','Name','QC Trial Passive Subtraction + Active Window');
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl,1); hold(ax1,'on');
    hAng = plot(ax1, time, angle, 'Color',[0 0.55 0], 'LineWidth',1.5);
    ylabel(ax1,'Angle (deg)','FontWeight','bold');
    set(ax1,'XTickLabel',[],'Box','off','TickDir','out','FontWeight','bold');

    ax2 = nexttile(tl,2); hold(ax2,'on');
    hVel = plot(ax2, time, angle_velocity, 'Color',[0.25 0.25 0.25], 'LineWidth',1.2);
    yline(ax2, velocity_thresh,  'm--', 'LineWidth',1.0, 'HandleVisibility','off');
    yline(ax2, -velocity_thresh, 'm--', 'LineWidth',1.0, 'HandleVisibility','off');
    ylabel(ax2,'Vel (deg/s)','FontWeight','bold');
    set(ax2,'XTickLabel',[],'Box','off','TickDir','out','FontWeight','bold');

    ax3 = nexttile(tl,3); hold(ax3,'on');
    hRaw = plot(ax3, time, torque_filt,   'LineWidth',1.2);
    hNet = plot(ax3, time, torque_active, 'LineWidth',1.6);

    ylabel(ax3,'Torque (Nm)','FontWeight','bold');
    xlabel(ax3,'Time (s)','FontWeight','bold');
    set(ax3,'Box','off','TickDir','out','FontWeight','bold');

    yl2 = ylim(ax2);
    yl3 = ylim(ax3);

    for r = 1:nReps
        idx0_rep = repIdx(r,1);
        idx1_rep = repIdx(r,2);

        if isfinite(idx0_rep) && isfinite(idx1_rep) && idx1_rep > idx0_rep && ...
                idx0_rep >= 1 && idx1_rep <= numel(time)

            if is_valid_rep(r)
                repColor = 'k';
            else
                repColor = [0.65 0.65 0.65];
            end

            % original rep window
            plot(ax2, [time(idx0_rep) time(idx0_rep)], yl2, '--', ...
                'Color', repColor, 'LineWidth',0.8, 'HandleVisibility','off');
            plot(ax2, [time(idx1_rep) time(idx1_rep)], yl2, '--', ...
                'Color', repColor, 'LineWidth',0.8, 'HandleVisibility','off');

            plot(ax3, [time(idx0_rep) time(idx0_rep)], yl3, '--', ...
                'Color', repColor, 'LineWidth',0.8, 'HandleVisibility','off');
            plot(ax3, [time(idx1_rep) time(idx1_rep)], yl3, '--', ...
                'Color', repColor, 'LineWidth',0.8, 'HandleVisibility','off');
        end

        % detected active/CV window
        if all(isfinite(activeWin_idx(r,:))) && activeWin_idx(r,2) > activeWin_idx(r,1)
            t0 = time(activeWin_idx(r,1));
            t1 = time(activeWin_idx(r,2));

            plot(ax2, [t0 t0], yl2, 'r--', 'LineWidth',1.2, 'HandleVisibility','off');
            plot(ax2, [t1 t1], yl2, 'r--', 'LineWidth',1.2, 'HandleVisibility','off');

            plot(ax3, [t0 t0], yl3, 'm--', 'LineWidth',1.2, 'HandleVisibility','off');
            plot(ax3, [t1 t1], yl3, 'm--', 'LineWidth',1.2, 'HandleVisibility','off');
        end
    end

    linkaxes([ax1 ax2 ax3],'x');

    sgtitle(sprintf('%s | %s | Trial passive subtraction QC', participantName, conditionName), ...
        'FontWeight','bold','Interpreter','none');

    hRepValid   = plot(ax3, NaN, NaN, 'k--', 'LineWidth', 1.0);
    hRepInvalid = plot(ax3, NaN, NaN, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0);
    hWinProxy   = plot(ax3, NaN, NaN, 'm--', 'LineWidth', 1.2);

    legend(ax3, [hAng hVel hRaw hNet hRepValid hRepInvalid hWinProxy], ...
        {'Angle','Angular velocity','Raw torque (filt)','Active torque', ...
         'Rep window (valid)','Rep window (invalid)','Detected active/CV window'}, ...
        'Location','best','Box','off');
end

%% ---------------- Save ----------------
saveName = sprintf('%s_%s_ActiveTorqueCorrected.mat', participantName, conditionName);
savePath = fullfile(path, saveName);

save(savePath, ...
    'participantName','conditionName','method','sourceFile', ...
    'time','Fs', ...
    'angle_raw','angle_filt','angle','angle_velocity', ...
    'torque_raw','torque_filt','torque_active', ...
    'passive_p','tau_pass_hat', ...
    'emg_raw','AUX','dac1', ...
    'time_matched','angle_matched','torque_matched','emg_raw_matched', ...
    'torque_raw_matched','torque_filt_matched','angle_raw_matched','angle_filt_matched', ...
    'repQC','is_valid_rep','valid_rep_idx', ...
    'angle_is_filtered','torque_is_filtered','angle_filter_Hz','torque_filter_Hz','filter_method', ...
    'velocity_thresh','min_duration', ...
    '-v7.3');

fprintf('✅ Saved trial corrected file: %s\n', saveName);
fprintf('Valid reps carried forward: %d / %d\n', nnz(is_valid_rep), nReps);
fprintf('Detected windows: %d / %d\n', nnz(all(isfinite(CVwin_idx),2)), nReps);










%% ===== DEBUG: Mean ACTIVE TORQUE across reps (CV-normalized) =====
% Goal:
%   - See torque behavior across the contraction cycle
%   - Identify where mismatch occurs (early vs late)

nPts = 101;  % match Stage 3 convention
torque_101 = nan(nReps, nPts);

for r = 1:nReps

    if ~is_valid_rep(r), continue; end

    idx0 = CVwin_idx(r,1);
    idx1 = CVwin_idx(r,2);

    if ~isfinite(idx0) || ~isfinite(idx1) || idx1 <= idx0
        continue
    end

    idx = idx0:idx1;

    t_seg = time(idx);
    tau_seg = torque_active(idx);

    % Normalize to 0–100% CV
    t_norm = linspace(t_seg(1), t_seg(end), nPts);
    tau_norm = interp1(t_seg, tau_seg, t_norm, 'linear');

    torque_101(r,:) = tau_norm;
end

% Keep only valid rows
valid_rows = all(isfinite(torque_101),2);
torque_101 = torque_101(valid_rows,:);

% Mean + SD
mean_tau = mean(torque_101,1);
sd_tau   = std(torque_101,0,1);

cv = linspace(0,100,nPts);

%% ===== Plot =====
figure('Color','w'); hold on;

% Shaded SD
fill([cv fliplr(cv)], ...
     [mean_tau+sd_tau fliplr(mean_tau-sd_tau)], ...
     [0.7 0.7 0.7], 'FaceAlpha',0.3, 'EdgeColor','none');

% Mean
plot(cv, mean_tau, 'k', 'LineWidth',2);

xlabel('Contraction Cycle (%)','FontWeight','bold');
ylabel('Active Torque (Nm)','FontWeight','bold');
title(sprintf('%s | %s | Mean ACTIVE Torque (CV-normalized)', ...
    participantName, conditionName), ...
    'Interpreter','none','FontWeight','bold');

set(gca,'Box','off','TickDir','out','FontWeight','bold');

%% ===== SCALAR EXTRACTION FROM MEAN ACTIVE TORQUE =====

% --- Define windows ---
win_full  = cv >= 0  & cv <= 100;
win_mid   = cv >= 10 & cv <= 80;   % recommended main window
win_early = cv >= 0  & cv <= 20;
win_late  = cv >= 80 & cv <= 100;

% --- Mean values ---
mean_full  = mean(mean_tau(win_full));
mean_mid   = mean(mean_tau(win_mid));
mean_early = mean(mean_tau(win_early));
mean_late  = mean(mean_tau(win_late));

% --- SD across CV (optional descriptive) ---
sd_full  = std(mean_tau(win_full));
sd_mid   = std(mean_tau(win_mid));
sd_early = std(mean_tau(win_early));
sd_late  = std(mean_tau(win_late));

%% ===== PRINT RESULTS =====
fprintf('\n===== ACTIVE TORQUE SCALARS (%s | %s) =====\n', participantName, conditionName);

fprintf('Full cycle (0–100%%):  %.2f Nm (SD=%.2f)\n', mean_full, sd_full);
fprintf('Mid phase  (10–80%%):  %.2f Nm (SD=%.2f)\n', mean_mid,  sd_mid);
fprintf('Early      (0–20%%):   %.2f Nm (SD=%.2f)\n', mean_early, sd_early);
fprintf('Late       (80–100%%): %.2f Nm (SD=%.2f)\n', mean_late,  sd_late);