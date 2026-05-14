%% --- Load synced file (EMG-timebase) ---
[file, path] = uigetfile('*_sync_rawVolts.mat', 'Select synchronized (EMG-timebase) file');
if isequal(file,0), return; end
S = load(fullfile(path, file));


%% Get passive torque
Frequency = 1 / mean(diff(S.time));
CutOff = 6;
Filt_order = 4;

%correct order and cut off
[cu,order] = Wcu(Filt_order,CutOff,'low',Frequency);
[b,a] = butter(order, cu/(0.5*Frequency),'low');

S.angle = filtfilt(b,a,S.angle);
S.torque = filtfilt(b,a,S.torque);
%now following filtering, get the automatic the polyfit of the resting
%torque
passive_p = Auto_RestTorque(S.torque,S.angle,Frequency);
%calculate passive torque and active torque using the previous caluclated
%polyfit values
[pasTor, Active_Tor] = Normalize_Torque(S.torque,S.angle,passive_p);


%% --- Save corrected trial with matched active torque
time_matched       = time_matched_all;
angle_matched      = angle_matched_all;
torque_matched     = torque_matched_all;
emg_raw_matched    = emg_matched_all;  % ⬅️ Synchronized raw EMG to 90% active torque

% Meta
torque_corrected   = tau_active_display;  % Full filtered active torque trace
angle              = angle_smooth;        % Full smoothed angle trace
time               = time;                % Full time vector
dt_align           = dt_ms / 1000;
method             = 'Passive fit + baseline correction + 90% central torque-EMG sync';
emg_raw            = emg;  % Full original EMG (unsliced)

if ~exist('participantName','var'), participantName = "Unknown"; end
if ~exist('conditionName','var'),  conditionName  = "Unknown";  end
if isfield(S,'aux_on'), aux_on = S.aux_on; else, aux_on = []; end
if isfield(S,'dac_on'), dac_on = S.dac_on; else, dac_on = []; end

saveName = sprintf('%s_%s_ActiveTorqueCorrected.mat', ...
                   string(participantName), string(conditionName));
savePath = fullfile(path, saveName);

save(savePath, ...
    'participantName','conditionName','Fs_emg','dt_align','method', ...
    'time','angle','torque_corrected','emg_raw', ...
    'time_matched','angle_matched','torque_matched','emg_raw_matched', ...
    'active_windows','aux_on','dac_on');

fprintf('✅ Saved corrected trial: %s\n', saveName);


