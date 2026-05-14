%% === Load synced file with EMG frequency and signal ===
[file, path] = uigetfile('*_sync_rawVolts.mat', '🔁 Select synced file with passive return');
if isequal(file,0)
    disp('❌ User canceled file selection.');
    return;
end
load(fullfile(path, file));  % Loads: time, angle, torque, Fs_emg, etc.

%% === Estimate Passive Torque from Opposite Contraction's Passive Return ===
delta = round(1.15 * Fs_emg);                % ~1.15s window in samples
TT_cut = [8.1, 18.1, 28.1, 38.1];            % Passive return start times (s)
indxes = dsearchn(time, TT_cut');            % Convert to sample indices

% --- Gather angle–torque data from each segment
ang_x = [];
tor_y = [];
for i = 1:length(indxes)
    idx_start = indxes(i);
    idx_end = min(idx_start + delta, length(time));  % prevent overflow
    ang_x = [ang_x; angle(idx_start:idx_end)];
    tor_y = [tor_y; torque(idx_start:idx_end)];
end

% --- Fit linear model to passive return
ang_x = ang_x(:);
tor_y = tor_y(:);
p = polyfit(ang_x, tor_y, 1);                 % Linear fit: tau ≈ m·angle + q

% --- Evaluate model on test range (optional visualization)
fake_ang = -3 : 1 : 37;
rest_tor = polyval(p, fake_ang);

% --- Print model equation
stringP = ["polyP CON: " + num2str(p)];
disp(stringP);
