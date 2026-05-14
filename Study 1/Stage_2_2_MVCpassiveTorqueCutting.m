%% ===== Stage 2.2: MVC PASSIVE TORQUE CORRECTION (MVC Stage 2) — CONSISTENT PIPELINE =====
% INPUT : *_sync_rawVolts.mat  (MVC recording; contains sync struct)
% OUTPUT: *_MVC_ActiveTorqueCorrected.mat
%
% - No ramp reference, no repQC.
% - Passive torque subtraction performed on full trial.
% - MVC window detected via angular velocity
% - QC plot: angle (green), torque_filt + torque_active (default colors), dashed window bounds

clear; clc;

%% --- User options (define explicitly) ---
useVelocityWindow = true;
useFixedWindow    = false;   % set true if you want fixed
search_start      = -Inf;     % seconds (global). set e.g. 0
search_end        = +Inf;     % seconds (global). set e.g. 20
fixed_t0          = 0;        % seconds
fixed_t1          = 5;        % seconds

velocity_thresh = 5;          % deg/s
min_duration    = 0.10;       % s

CutOff_torque = 20;
CutOff_angle  = 6;

doQCplot = true;

%% --- Load MVC synced file (sync struct) ---
[file, path] = uigetfile('*_sync_rawVolts.mat', 'Select MVC synchronized file (*_sync_rawVolts.mat)');
if isequal(file,0), return; end
L = load(fullfile(path, file));

if ~isfield(L,'sync')
    error('Selected file does not contain "sync" struct.');
end
sync = L.sync;

%% --- Core signals (same naming as trials) ---
time   = sync.time(:);
torque_raw = sync.torque_raw(:);
angle_raw  = sync.angle(:);

if isfield(sync,'emg')
    emg_raw = sync.emg;
else
    emg_raw = [];
end

Fs = 1 / median(diff(time));  % consistent with your other stages

%% --- Filtering (same logic as trials) ---
torque_filt = real(Wfilt(torque_raw, CutOff_torque, 'low', Fs));
angle_filt  = real(Wfilt(angle_raw,  CutOff_angle,  'low', Fs));

%% --- Passive subtraction (unchanged) ---
passive_p = Auto_RestTorque(torque_filt, angle_filt, Fs);
[tau_pass_hat, tau_active] = Normalize_Torque(torque_filt, angle_filt, passive_p);

% Keep consistent naming
torque_active = tau_active;   % passive removed torque (net)
angle         = angle_filt;   % use filtered angle downstream

%% --- Define MVC window (time-based, then convert to indices) ---
if useVelocityWindow
    % angular velocity from filtered angle
    angle_velocity = [0; diff(angle) ./ diff(time)];

    idx_search = find(time >= search_start & time <= search_end);
    if isempty(idx_search)
        warning('Search window empty; falling back to fixed window.');
        useVelocityWindow = false; useFixedWindow = true;
    else
        v = angle_velocity(idx_search);
        t = time(idx_search);

        above = abs(v) > velocity_thresh;
        d = diff([0; above; 0]);
        starts = find(d==1);
        ends   = find(d==-1)-1;

        mvc_t0 = NaN; mvc_t1 = NaN;

        for j = 1:numel(starts)
            t0 = t(starts(j));
            t1 = t(ends(j));
            if (t1 - t0) >= min_duration
                mvc_t0 = t0;
                mvc_t1 = t1;
                break
            end
        end

        if ~isfinite(mvc_t0) || ~isfinite(mvc_t1)
            warning('⚠️ No valid velocity-based MVC window detected. Falling back to fixed window.');
            useVelocityWindow = false; useFixedWindow = true;
        end
    end
end

if useFixedWindow
    mvc_t0 = fixed_t0;
    mvc_t1 = fixed_t1;
end

if ~(isfinite(mvc_t0) && isfinite(mvc_t1) && mvc_t1 > mvc_t0)
    error('MVC window invalid (mvc_t0/mvc_t1). Check your settings.');
end

% Convert time bounds to indices (global)
idx0 = find(time >= mvc_t0, 1, 'first');
idx1 = find(time <= mvc_t1, 1, 'last');

if isempty(idx0) || isempty(idx1) || idx1 <= idx0
    error('MVC window invalid after conversion to indices. Check mvc_t0/mvc_t1.');
end

%% --- Extract MVC segment (same naming as trials, but vectors not cells) ---
idx = idx0:idx1;

time_matched   = time(idx);
angle_matched  = angle(idx);
torque_matched = torque_active(idx);

if ~isempty(emg_raw)
    emg_raw_matched = emg_raw(idx,:);
else
    emg_raw_matched = [];
end

%% --- Flags consistent with pipeline ---
angle_is_filtered  = true;
torque_is_filtered = true;
angle_filter_Hz    = CutOff_angle;
torque_filter_Hz   = CutOff_torque;
filter_method      = "Wfilt low-pass";

method = 'MVC passive fit (Auto_RestTorque + Normalize_Torque) + MVC window (velocity/fixed)';

%% --- QC Plot: Raw vs Active Torque + Angle + MVC window ---
if doQCplot
    figure('Color','w','Name','QC MVC Passive Subtraction (Raw vs Active) + MVC window');
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl,1); hold(ax1,'on');
    hAng = plot(ax1, time, angle, 'Color',[0 0.55 0], 'LineWidth',1.5);
    ylabel(ax1,'Angle (deg)','FontWeight','bold');
    set(ax1,'XTickLabel',[],'Box','off','TickDir','out','FontWeight','bold');

    ax2 = nexttile(tl,2); hold(ax2,'on');
    hRaw = plot(ax2, time, torque_filt,   'LineWidth', 1.2);   % filtered raw
    hNet = plot(ax2, time, torque_active, 'LineWidth', 1.6);   % net

    ylabel(ax2,'Torque (Nm)','FontWeight','bold');
    xlabel(ax2,'Time (s)','FontWeight','bold');
    set(ax2,'Box','off','TickDir','out','FontWeight','bold');

    y_lim_fixed = ylim(ax2);
    t0 = time(idx0); t1 = time(idx1);
    plot(ax2, [t0 t0], y_lim_fixed, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
    plot(ax2, [t1 t1], y_lim_fixed, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');

    uistack([hRaw hNet],'top');
    linkaxes([ax1 ax2],'x');

    if isfield(sync,'meta') && isfield(sync.meta,'participantName')
        participantName = string(sync.meta.participantName);
    else
        participantName = "Unknown";
    end
    if isfield(sync,'meta') && isfield(sync.meta,'conditionName')
        conditionName = string(sync.meta.conditionName);
    else
        conditionName = "MVC";
    end

    sgtitle(sprintf('%s | %s | MVC Passive subtraction QC', participantName, conditionName), ...
        'FontWeight','bold','Interpreter','none');

    hWinProxy = plot(ax2, NaN, NaN, 'k--', 'LineWidth', 1.2);
    legend(ax2, [hAng hRaw hNet hWinProxy], ...
        {'Angle','Raw torque (filt)','Active torque','MVC window'}, ...
        'Location','best', 'Box','off');
end

%% --- Save (MVC-specific filename) ---
if isfield(sync,'meta') && isfield(sync.meta,'participantName')
    participantName = string(sync.meta.participantName);
else
    participantName = "Unknown";
end
if isfield(sync,'meta') && isfield(sync.meta,'conditionName')
    conditionName = string(sync.meta.conditionName);
else
    conditionName = "MVC";
end

saveName = sprintf('%s_%s_MVC_ActiveTorqueCorrected.mat', participantName, conditionName);
savePath = fullfile(path, saveName);

sourceFile = fullfile(path, file);  % provenance

save(savePath, ...
    'participantName','conditionName','method','sourceFile', ...
    'time','Fs', ...
    'angle_raw','angle_filt','angle', ...
    'torque_raw','torque_filt','torque_active', ...
    'passive_p','tau_pass_hat', ...
    'emg_raw', ...
    'time_matched','angle_matched','torque_matched','emg_raw_matched', ...
    'angle_is_filtered','torque_is_filtered','angle_filter_Hz','torque_filter_Hz','filter_method', ...
    'mvc_t0','mvc_t1','idx0','idx1', ...
    'velocity_thresh','min_duration','useVelocityWindow','useFixedWindow','search_start','search_end','fixed_t0','fixed_t1', ...
    '-v7.3');

fprintf('✅ Saved MVC corrected trial: %s\n', saveName);