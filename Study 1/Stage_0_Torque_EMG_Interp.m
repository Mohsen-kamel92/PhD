
%% ===== Stage 0: Interpolate SPIKE Torque/Angle to EMG Timebase + Align DAC1 vs AUX + Save Sync File =====
% Summary of what this script does:
% 1) Load SPIKE (.mat) containing DAC1/Torque/Angle, crop to DAC1>0.85 window, re-zero SPIKE time.
% 2) Load MAECS (HD-EMG), crop to AUX>0 window, re-zero EMG time.
% 3) Estimate constant onset offset (DAC1 vs AUX), shift SPIKE time by -dt.
% 4) Interpolate SPIKE torque/angle/DAC1 onto EMG timebase (no extrapolation).
% 5) Crop all signals to common valid overlap (remove NaN edges).
% 6) Compute DAC/AUX onset markers on the cropped timeline (for QC + later reps).
% 7) Plot QC figure (angle + EMG RMS mean + torque) with onset markers.
% 8) Save synchronized dataset: raw EMG + aligned torque/angle + triggers + metadata.
clear; close all;

%% --- Load SPIKE .mat (expects fields: DAC1, Torque, Angle with .times/.values) ---
[SpikeFile, pathFile] = uigetfile('*.mat','Pick SPIKE mat file','MultiSelect','off'); % off allows just one file selection
if SpikeFile == 0 % Stop if user cancels Picking Spike File
    return
end
S = load(fullfile(pathFile, SpikeFile)); %build a struct and store it is S.

mustHave = {'DAC1','Torque','Angle'};
for k = 1:numel(mustHave)
    if ~isfield(S, mustHave{k}) %if S does not have the mustHave's fielsd, return false=true, then excute the error function directly
        error('Missing "%s" in SPIKE file.', mustHave{k});
    end
end

%% --- Crop SPIKE by DAC1 > 0.85 (first to last), re-zero time ---
th_dac = 0.85;
dac_values = S.DAC1.values(:); % Access the field values in the DAC1, and return it as column array using (:) colon operator 

i1 = find(dac_values > th_dac, 1, 'first'); %Find the first 1 index where dac_values > th_dac is true. 1> means return only one match  
i2 = find(dac_values > th_dac, 1, 'last');
if isempty(i1) || isempty(i2) || i2 <= i1 % || is scalar logical OR (short-circuit) used inside if statement. | used inside arrays operates emlement by element
    error('Could not find a valid high segment in DAC1 using threshold %.2f.', th_dac);
end

SPIKE_vars = fieldnames(S); % store the S struct fields names in this <
for i = 1:numel(SPIKE_vars) % loop through this fields 
    fn = SPIKE_vars{i}; % now return spike vars we stored whos fields name is {i}. Also, it is a dynamic field referencing S.(fn)
    if isstruct(S.(fn)) && isfield(S.(fn),'times') && isfield(S.(fn),'values') % Here we check if fn just has the struct fields not anything else
        S.(fn).times  = S.(fn).times(i1:i2); %crop the times values to return just range from i1 to i2
        S.(fn).values = S.(fn).values(i1:i2); % same here
        S.(fn).times  = S.(fn).times - S.(fn).times(1); % rezero the field times to start  from zero
    end
end

% % SPIKE nominal frequency (if present)
% if isfield(S.DAC1,'interval')
%     Fs_spike = 1 / S.DAC1.interval;
% elseif isfield(S.Torque,'interval')
%     Fs_spike = 1 / S.Torque.interval;
% else
%     Fs_spike = NaN; %Not-a-Number, so Fs could not be determined so far
% end

Fs_spike = 2000; %~2k Hz
%% --- Load MAECS (HDEMG) and crop by AUX > 0, re-zero time ---
Fs_emg = 2048;
HDEMG  = MAECS_read();

gNames = fieldnames(HDEMG);
assert(isscalar(gNames), 'Expected exactly one grid, found %d.', numel(gNames)); % isscala expected gNames should be 1×1 scalar, otherwise, errpr.
gridName = gNames{1};

if ~isfield(HDEMG.(gridName), 'AUX') %check wether AUX is here
    error('AUX not found in HDEMG.%s', gridName);
end

if ~isfield(HDEMG.(gridName), 'Times') %if Times field is not here through an error
    error('No Time vector found in HDEMG.%s.', gridName);
end

if     isfield(HDEMG.(gridName), 'Times') % and if found, return it to timeEMG
    timeEMG = HDEMG.(gridName).Times(:);
end

aux = HDEMG.(gridName).AUX(:); % force AUX to be transposed (:) to column
emg = HDEMG.(gridName).EMG;   %  already has been stored as [samples x channels] from MAESC_read func

assert(numel(aux) == numel(timeEMG), 'AUX length must match time length.'); % numel(aux) → total number of elements in AUX
assert(size(emg,1) == numel(timeEMG), 'EMG rows must match time length.'); % same for emg
% Return the size of array A along dimension dim. size(A,dim) So:  
% size(emg,1) → number of rows> samples
% size(emg,2) → number of columns > Channels

aux(~isfinite(aux)) = 0; %replaces all NaN and Inf values in aux with zero using logical indexing.

j1 = find(aux > 0, 1, 'first'); % find just the 1st value of aux greater than 0 as we already sure that's now  0 1 logical-based
j2 = find(aux > 0, 1, 'last'); % same here as last value of aux
if isempty(j1) || isempty(j2) || j2 <= j1 % check wether the Js returned full, and they are inexing in the correct order
    error('Could not find a valid high segment in AUX (need AUX>0 span).');
end

emg_vals = emg(j1:j2, :); %Take rows from index j1 to j2 > Keep all channels : as they are dim 2, So this extracts a time segment of the EMG.
emg_time = timeEMG(j1:j2); % As time is NX1, cut it from j1 to j2 as well
emg_time = emg_time - emg_time(1); % just shift emg_time label to start at 0
AUX_cut  = aux(j1:j2); % now do the same for Aux, so emgVals, emgTime, and Aux are all aligned and ready for next step!
% Note! we now have the same time onset label for DAC1 and emg_time, which
% ease aligning their signals latter...

%% --- Onset alignment (DAC1 vs AUX), interpolate SPIKE -> EMG timebase ---
tAux0 = emg_time(find(AUX_cut > 0, 1, 'first')); % Time of first AUX activation (onset), so the first time point corresponding to first Aux value
tDAC0 = S.DAC1.times(find(S.DAC1.values > th_dac, 1, 'first')); % Time of first DAC1 threshold crossing, so same for DAC1, but using the fixed threshold

dt_align = tDAC0 - tAux0;    % given that (+ => DAC lags AUX) according to the difference in Fs BTW Aux and DAC1
tTorque = S.Torque.times - dt_align; % now subtract this dt from all DAC1 signals
tAngle  = S.Angle.times  - dt_align;  % same here 
tDAC    = S.DAC1.times   - dt_align;  % same here

torque_i = interp1(tTorque, S.Torque.values, emg_time, 'linear');  % no extrap
angle_i  = interp1(tAngle,  S.Angle.values,  emg_time, 'linear');
dac1_i   = interp1(tDAC,    S.DAC1.values,   emg_time, 'linear');

% Explaination of iterp1 function: interp1(x, y, xq, method), whereas, x =
% original time vec, y = original signal values, xq = the new time base
% to estimate the new yq values linearly. And linear = is the method used
% assuming that: BTW two points, we assume a line, so (x1, y1) & (x2, y2)
% So, formula is like: If x1 =< xq =< x2, so 
% yq = y1 + xq - x1 / x2 - x1  (y2 -y1)
%% --- Crop to common valid overlap (remove NaN edges) ---
valid = isfinite(torque_i) & isfinite(angle_i) & isfinite(dac1_i); % Only true if all signals are normal number,
% so a protection against NaN 
firstIdx = find(valid, 1, 'first');
lastIdx  = find(valid, 1, 'last');

if isempty(firstIdx) || isempty(lastIdx) || lastIdx <= firstIdx % wether 1st or 2end index is not here, stop and through the error
    error('No valid overlap found after interpolation. Check dt_align and trigger windows.');
end

keep = firstIdx:lastIdx; % return the first and last index to KEEP , so the clean region of the normal signal (nothing NaN)

time     = emg_time(keep); % Trim all signal according to KEEP Nx1
emg_vals = emg_vals(keep,:); % same here, but note that we use keep for rows i.e. samples X all channel using indexing colon operator (:), so preserving the 32 base
AUX_cut  = AUX_cut(keep); % same Nx1
torque_i = torque_i(keep); % same Nx1
angle_i  = angle_i(keep); % same Nx1
dac1_i   = dac1_i(keep); % same Nx1

%% --- Compute onset markers on CROPPED timeline ---
thr_dac_plot = th_dac;
dac_on = [time(1); time(find(diff(dac1_i > thr_dac_plot) == 1) + 1)];
% create logical vector where 
% DAC1 greater than th, then true is 1s, otherwise; false is 0s
% then compute the difference BTW each consecutive elements in this
% true-false logical array, then detect the rising moments using diff(DAC
% true) == 1, then do find(diff(trueDAC1) + 1 because diff shortens the vector
% by 1 sample, so adding +1 for reshifting the onset to the true position in the
% logical array. then, time now got that true signal's position > then
% time(1) returns just the 1st arised signal of dac1 to dac_on

aux_on = [time(1); time(find(diff(AUX_cut > 0) == 1) + 1)];
% we are doing the same for Aux, but note that Aux is compared to 0 sice
% its th is a digital not an analog signal as per dac1

if isempty(dac_on), dac_on = time(1); end  % If no onset was detected → force onset = first time sample.
if isempty(aux_on), aux_on = time(1); end % Same here to avoid errors of indexing latter on

% Offset between first onsets (BTW Aux and DAC1) on cropped time
tAux0  = time(find(AUX_cut > 0,            1, 'first'));
tDAC0  = time(find(dac1_i  > thr_dac_plot, 1, 'first'));

dt_ms = (tDAC0 - tAux0) * 1e3;  % we convert it to ms by * 1000

%% Report the time difference and the kept total segmented duration
fprintf('dt_align = %.3f ms; kept %.3f s overlap.\n', dt_align*1e3, time(end)-time(1)); % we are accessing 1st time value as (1)
% and last as (end) to report total time of condition set

%% --- QC plot (smoothed torque/angle for display only; EMG RMS mean) ---
% NOTE: Only for visualization. Saved torque/angle remain unfiltered (torque_i/angle_i).

% angle and torque smoothing (butterworth) processed by Wfilt + Wcu order
% and fc correction to keep zero-phase shift filtering
fc_angle  = 6;
fc_torque = 20;
angle_smooth  = Wfilt(angle_i,  fc_angle,  'low', Fs_emg);
torque_smooth = Wfilt(torque_i, fc_torque, 'low', Fs_emg);


% EMG signal smoothed over Nchannel then mean across channels
win_sec = 0.050;
win = max(1, ceil(win_sec*Fs_emg)); % convert window to N samples, so that the window length of samples should always be integer number. ceil round upward the nearst integer number, and max(1)
emg_env = sqrt(movmean(emg_vals.^2, win, 1, 'Endpoints','shrink')); % Prevents edge distortion by using smaller windows at the beginning and end 
emg_env_mean = mean(emg_env, 2, 'omitnan'); % mean across channels, D2



%% plot tiled layout

%Create a 2-row tiled figure (top: Angle, bottom: EMG+Torque with two y-axes)
%Add vertical event markers (DAC and AUX) on both panels
%Make the figure clean and add a legend (using “proxy” handles so the legend looks right)


% ==========Upper panel

tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact'); % 2 rows x 1 column, reduce the space BTW, and pading boarders
ax1 = nexttile(1); %activates tile #1 (top plot) and returns its axes handle
hAng = plot(time, angle_smooth, ...
    'Color',[0 0.55 0], ...  %'Color' system used, [R G B] at sacale of 0 to 1
    'LineWidth', 1.8); 
hold on;
ylabel('Angle (°)','FontWeight','bold');
set(ax1,'XTickLabel',[],'XColor','none'); %because the bottom plot will show the time axis; top plot stays clean
if ~isempty(dac_on)

    % First DAC line (visible in legend)
    xline(ax1, dac_on(1), '--k', ...
        'LineWidth', 1.2, ...
        'HandleVisibility','on');

    % Remaining DAC lines (hidden from legend)
    for i = 2:length(dac_on)
        xline(ax1, dac_on(i), '--k', ...
            'LineWidth', 1.2, ...
            'HandleVisibility','off');
    end

end
if ~isempty(aux_on)

    % First AUX line (visible in legend)
    xline(ax1, aux_on(1), '--', ...
        'Color',[1 0 1], ...
        'LineWidth', 1.2, ...
        'HandleVisibility','on');

    % Remaining AUX lines (hidden from legend)
    for i = 2:length(aux_on)
        xline(ax1, aux_on(i), '--', ...
            'Color',[1 0 1], ...
            'LineWidth', 1.2, ...
            'HandleVisibility','off');
    end

end



%========== Lower panel

%----EMG>
ax2 = nexttile(2);
yyaxis left
hEMG = plot(time, emg_env_mean, 'r', 'LineWidth', 1.8); hold on;
ylabel('EMG RMS mean (V)','FontWeight','bold');
ax2.YAxis(1).Color = [1 0 0];   % EMG axis in red
yyaxis right
%----Torque>
hTorque = plot(time, torque_smooth, 'Color',[0.5 0 0.5], 'LineWidth', 1.8);
ylabel('Torque (Nm)','FontWeight','bold');
ax2.YAxis(2).Color = [0.5 0 0.5];   % Torque axis in purple
xlabel('Time (s)','FontWeight','bold');


% ======== DAC & AUX markers on ax2 (bottom panel)
yyaxis(ax2,'left')

if ~isempty(dac_on)
    for i = 2:length(dac_on)
        xline(ax2, dac_on(i), '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
    end
end

if ~isempty(aux_on)
    for i = 2:length(aux_on)
        xline(ax2, aux_on(i), '--k', 'LineWidth', 1.2, 'Color',[1 0 1], 'HandleVisibility','off');
    end
end
   
%=====Global formatting
linkaxes([ax1,ax2],'x'); %Links the x-axes: zoom/pan one → both move together
set(gcf,'Color','w'); %“get current figure” weight color
set([ax1,ax2], 'Color','w','XGrid','off','YGrid','off','Box','off', ...
    'TickDir','out','FontWeight','bold','FontSize',12,'LineWidth',1.2);

%=====Title of the figure
sgtitle(sprintf('Trigger Alignment (DAC1 - AUX = %.2f ms)', dt_ms), 'FontWeight','bold'); %super group title” for the whole tiled layout
%% --- Legend (with proxy handles to keep box and colors correct) ---
hAngProxy = plot(ax2, NaN, NaN, '-',  'Color',[0 0.55 0], 'LineWidth',1.8);
hTQProxy  = plot(ax2, NaN, NaN, '-',  'Color',[0.5 0 0.5], 'LineWidth',1.8);
hEMGProxy = plot(ax2, NaN, NaN, '-',  'Color',[1 0 0],   'LineWidth',1.8);
hDACproxy = plot(ax2, NaN, NaN, '--', 'Color',[0 0 0],   'LineWidth',1.2);
hAUXproxy = plot(ax2, NaN, NaN, '--', 'Color',[1 0 1],   'LineWidth',1.2);

lgd = legend(ax2, ...
    [hAngProxy hTQProxy hEMGProxy hDACproxy hAUXproxy], ...
    {'Angle','Torque','EMG RMS-mean','DAC1','AUX'}, ...
    'Location','best', 'Box','off', 'FontSize',8);

lgd.ItemTokenSize = [10 6];


%% Folder structure assumed:
% .../Participant/CON_75/Set_1/

[setPath, setFolder] = fileparts(pathFile(1:end-1));  %parent folder path, folder name      % e.g Set_1, so removes the last character (/) 
[condPath, conditionFolder] = fileparts(setPath);          % > CON_75
[~, participantFolder] = fileparts(condPath);              % > Participant
% So, % setPath   = 'D:\Data\P01\CON_75'
      % setFolder = 'Set_1'
participantName = participantFolder;
conditionName   = conditionFolder;   % <-- CORRECT LEVEL
setName         = setFolder;         % <-- store set separately


% --- Define method string (so it exists) ---
method = 'onset-shift + interp-to-EMG + crop-overlap (no extrap)';

% --- Save synchronized dataset as a single struct ---
sync = struct(); %create a struct to store everything!

sync.time       = time;
sync.Fs_emg     = Fs_emg;
sync.Fs_spike   = Fs_spike;

sync.emg        = emg_vals;
sync.AUX        = AUX_cut;

sync.torque_raw     = torque_i;
sync.angle          = angle_i;
sync.dac1           = dac1_i;

sync.aux_on     = aux_on;
sync.dac_on     = dac_on;

sync.meta = struct();
sync.meta.method      = method;
sync.meta.thr_dac     = th_dac;
sync.meta.dt_align_s  = dt_align;
sync.meta.dt_align_ms = dt_align*1e3;

sync.meta.participantName = participantName;
sync.meta.conditionName   = conditionName;

sync.meta.spike_file = SpikeFile;
sync.meta.spike_path = pathFile;

sync.meta.crop = struct();
sync.meta.crop.spike_i1 = i1;
sync.meta.crop.spike_i2 = i2;
sync.meta.crop.emg_j1   = j1;
sync.meta.crop.emg_j2   = j2;
sync.meta.crop.overlap_firstIdx = firstIdx;
sync.meta.crop.overlap_lastIdx  = lastIdx;

outname = fullfile(pathFile, sprintf('%s_%s_sync_rawVolts.mat', participantName, conditionName));

try
    save(outname, 'sync', '-v7.3'); % -v7.3 Allows large variables.
catch
    save(outname, 'sync');
end

fprintf('✅ Saved synchronized struct to: %s\n', outname);
