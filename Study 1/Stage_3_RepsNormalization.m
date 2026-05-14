%% ===== Stage 3: Bipolar EMG RMS (50 ms) + MVC normalization + CV extraction =====
% Outputs are rep-level (NO collapsing). SPM1D-ready (101 points) and reviewer-traceable.
% - Bipolar-first (32 mono -> 28 bipolars) using validated ChMap
% - Band-pass on bipolars (20–450 Hz, 4th order, zero-phase)
% - RMS envelope (50 ms ~ 103 samples at 2048 Hz)
% - MVC reference per bipolar channel: BOTH ref95 and refMax saved
% - Trial normalization uses either ref95 or refMax (flag useMVCpercentile)
% - CV windows are reused from repQC (CVwin_idx or CVwin_time -> idx)

clc; clear; close all;
%% --- Load trial (Stage 2) ---
[fileT, pathT] = uigetfile('*_ActiveTorqueCorrected.mat', 'Select TRIAL ActiveTorqueCorrected (Stage 2)');
if isequal(fileT,0), return; end
T = load(fullfile(pathT, fileT));

assert(isfield(T,'repQC'), 'Trial file missing repQC.');
assert(isfield(T,'emg_raw'), 'Trial file missing emg_raw.');
assert(isfield(T,'time') && isfield(T,'angle') && isfield(T,'torque_active'), ...
    'Trial file missing time/angle/torque_active.');

repQC = T.repQC;
time  = T.time(:);
angle = T.angle(:);
torque_active = T.torque_active(:);   % passive removed (Stage 2)
emg_raw = T.emg_raw;                        % [N x 32] monopolar raw

Fs_trial = 1 / median(diff(time), 'omitnan');

%% --- Load MVC (MVC Stage 2) ---
[fileM, pathM] = uigetfile('*_MVC_ActiveTorqueCorrected.mat', 'Select MVC ActiveTorqueCorrected');
if isequal(fileM,0), return; end
M = load(fullfile(pathM, fileM));

assert(isfield(M,'emg_raw_matched'), 'MVC file missing emg_raw_matched (MVC window).');
emg_mvc_raw = M.emg_raw_matched;            % [Nmvc x 32] within MVC window

Fs_mvc = Fs_trial;
if isfield(M,'Fs_emg') && ~isempty(M.Fs_emg) && isfinite(M.Fs_emg)
    Fs_mvc = M.Fs_emg;
end

% Safety: enforce same Fs for MVC and trials (recommended)
if abs(Fs_mvc - Fs_trial) > 1e-6
    error('Fs mismatch: Trial Fs=%.6f Hz, MVC Fs=%.6f Hz. Resample or ensure same Fs.', Fs_trial, Fs_mvc);
end

Fs = Fs_trial;

%% --- PARAMETERS ---
% Bipolar + filtering
bp_low   = 20;    % Hz
bp_high  = 450;   % Hz
bp_order = 4;     % 4th order Butterworth. it's been corrected in Wcu function to guarantee zero-phase shift and doual filteration effect

% RMS envelope
rms_win_ms   = 50;
rms_win_samp = max(1, round((rms_win_ms/1000) * Fs));
if mod(rms_win_samp,2)==0, rms_win_samp = rms_win_samp + 1; end % if window is EVEN make it ODD to better align with centered movement

% MVC reference definition (used for normalization)
useMVCpercentile = true;   % true -> use ref95; false -> use refMax
mvc_prc = 95;              % percentile for ref95 (robust)

% Time normalization points (SPM1D-ready)
nTimeNorm = 101;
x_new = linspace(0,1,nTimeNorm);
%% --- Channel map (validated snake layout) ---
% PROXIMAL row 1 (top), DISTAL row 8 (bottom); Lateral left, Medial right.
ChMap = [1 16 17 32;
         2 15 18 31;
         3 14 19 30;
         4 13 20 29;
         5 12 21 28;
         6 11 22 27;
         7 10 23 26;
         8  9 24 25];

[nRows, nCols] = size(ChMap);
nBipolar = (nRows-1) * nCols; % (8-1)*4 = 28

% Documentation flags for later centroid interpretation
bipolar_direction = "distal_minus_proximal"; % bipolar = row(i) - row(i+1)

% Precompute bipolar pairs and labels (saved for traceability)
bipolar_pairs  = nan(nBipolar,2);  % [28 x 2] = derivations are between [proxCh distCh]
bipolar_labels = strings(nBipolar,1); % put a label for each ch, e.g. 'CH1'
bipolar_rowcol = nan(nBipolar,2);  % [row col] to know "location", e.g. between row and row+1

k = 1;
for col = 1:nCols
    for row = 1:(nRows-1)
        ch1 = ChMap(row, col);       % proximal electrode
        ch2 = ChMap(row+1, col);     % distal electrode
        bipolar_pairs(k,:)  = [ch1 ch2];
        bipolar_labels(k)   = sprintf('R%dC%d', row, col);
        bipolar_rowcol(k,:) = [row col];
        k = k + 1; %move to the next row in the storage arrays
    end
end

%% --- Helper: compute bipolar -> bandpass -> RMS ---
fnyq = Fs/2;
assert(bp_high < fnyq, 'bp_high must be < Nyquist.');

[b_bp, a_bp] = butter(bp_order, [bp_low bp_high]/fnyq, 'bandpass');
compute_bipolar_bp_rms = @(emg32) local_bipolar_bp_rms(emg32, bipolar_pairs, b_bp, a_bp, rms_win_samp);

%% --- MVC: bipolar -> BP -> RMS, then references per bipolar channel ---
[bip_mvc_rms, ~] = compute_bipolar_bp_rms(emg_mvc_raw);  % [Nmvc x 28]

MVC_ref95  = nan(1,nBipolar);
MVC_refMax = nan(1,nBipolar);

for ch = 1:nBipolar
    x = bip_mvc_rms(:,ch);
    MVC_ref95(ch)  = prctile(x, mvc_prc);
    MVC_refMax(ch) = max(x, [], 'omitnan');
end

% Choose denominator used in normalization
if useMVCpercentile
    MVC_ref_used = MVC_ref95;
else
    MVC_ref_used = MVC_refMax; %#ok<*UNRCH>
end

% Guard against invalid denominators
MVC_ref_used(MVC_ref_used<=0 | ~isfinite(MVC_ref_used)) = NaN;
MVC_ref95(MVC_ref95<=0 | ~isfinite(MVC_ref95)) = NaN;
MVC_refMax(MVC_refMax<=0 | ~isfinite(MVC_refMax)) = NaN;

% QC summary for MVC reference stability
ratio95 = MVC_ref95 ./ MVC_refMax;
fprintf('\n=== MVC reference QC (bipolar RMS) ===\n');
fprintf('ref95/refMax across channels: median=%.3f, IQR=[%.3f %.3f]\n', ...
    median(ratio95,'omitnan'), prctile(ratio95,25), prctile(ratio95,75)); % check Interquartile range, the range of midle 50% of channels
fprintf('Channels with ref95/refMax < 0.75: %d / %d\n', sum(ratio95 < 0.75, 'omitnan'), nBipolar); 
% check whether 95th percentile MVC reference (ref95) is close to the true maximum (refMax) for each bipolar channel

%% --- Ensure CV indices exist (time->idx conversion if needed) ---
assert(isfield(repQC,'repIdx') && ~isempty(repQC.repIdx), 'repQC.repIdx missing/empty.');
nReps = size(repQC.repIdx,1);

% Prefer CV indices (should already exist from Stage-2)
if isfield(repQC,'CVwin_idx') && ~isempty(repQC.CVwin_idx)
    CV_idx = repQC.CVwin_idx;

else
    error(['repQC.CVwin_idx is missing/empty. ', ...
           'Re-run Stage-2 (or add the CVwin_time->idx conversion there) and re-save repQC.']);
end

% validity mask (default all true if not present)
if isfield(repQC,'is_valid_rep') && ~isempty(repQC.is_valid_rep)
    is_valid_rep = repQC.is_valid_rep(:);
else
    is_valid_rep = true(nReps,1);
end

%% --- Trial: compute bipolar/BP/RMS for FULL trial once (efficient) ---
[bip_trial_rms, ~] = compute_bipolar_bp_rms(emg_raw);  % [N,time.S x 28]. we use raw_em to  filter once for whole trial then just index the segments

% Normalize per bipolar channel (%MVC) using chosen reference
% (implicit expansion requires R2016b+; otherwise use bsxfun)
bip_trial_rms_norm = (bip_trial_rms ./ MVC_ref_used) * 100;  % [N x 28] normalized



%% --- Stage 3 output struct ---
stage3 = struct();

% Metadata / provenance
stage3.meta = struct();
if isfield(T,'participantName'), stage3.meta.participantName = T.participantName; end
if isfield(T,'conditionName'),   stage3.meta.conditionName   = T.conditionName; end
stage3.meta.trialFile = fullfile(pathT, fileT);
stage3.meta.mvcFile   = fullfile(pathM, fileM);
stage3.meta.createdOn = datestr(now); %#ok<TNOW1,DATST> % returns the current date and time as a serial date number
% Added provenance info
stage3.meta.script = mfilename;   % name of the script that created this file
stage3.meta.matlab = version;     % MATLAB version used
% Parameters
stage3.params = struct( ...
    'Fs', Fs, ...
    'bp_low', bp_low, ...
    'bp_high', bp_high, ...
    'bp_order', bp_order, ...
    'rms_win_ms', rms_win_ms, ...
    'rms_win_samp', rms_win_samp, ...
    'mvc_prc', mvc_prc, ...
    'useMVCpercentile', useMVCpercentile, ...
    'mvc_ref_used_type', ternary(useMVCpercentile, sprintf('prctile%d',mvc_prc), 'max'), ...
    'timeNorm_points', nTimeNorm, ...
    'timeNorm_axis', x_new, ...      % <-- store the normalized CV axis
    'interp_method', "linear", ...
    'torque_units', "Nm (passive-removed)", ...
    'angle_units', "deg (filtered in Stage 2)", ...
    'torque_norm', "none (absolute)", ...
    'angle_norm', "none (absolute)");

% Mapping / geometry
stage3.map = struct( ...
    'ChMap', ChMap, ...
    'bipolar_pairs', bipolar_pairs, ...
    'bipolar_labels', bipolar_labels, ...
    'bipolar_rowcol', bipolar_rowcol, ...
    'bipolar_direction', bipolar_direction, ...
    'proximal_row', 1, ...
    'distal_row', nRows, ...
    'lateral_column', 1, ...
    'medial_column', nCols, ...
    'row_meaning', "1 = proximal, 8 = distal", ...
    'column_meaning', "1 = lateral, 4 = medial");


% MVC references (+ optional RMS time series for QC/repro)
stage3.mvc = struct( ...
    'MVC_ref_used', MVC_ref_used, ...
    'MVC_ref95', MVC_ref95, ...
    'MVC_refMax', MVC_refMax, ...
    'ratio_ref95_refMax', ratio95, ...
    'bip_mvc_rms', bip_mvc_rms);   % optional; comment out if file size is too large


% Rep segmentation indices (traceability)
stage3.repIdx = repQC.repIdx;

% Rep-level outputs
stage3.rep = repmat(struct( ...
    'is_valid', false, ...
    'CV_idx', [NaN NaN], ...
    'torque_mean_CV', NaN, ...
    'angle_mean_CV', NaN, ...
    'emgRMS_mean_CV', NaN, ...
    'emgRMS_mean_CV_normMVC', NaN, ...
    'torque_101', [], ...
    'angle_101', [], ...
    'emgRMSmean_101', [], ...
    'emgRMSmeanNorm_101', [], ...
    'emgRMS_bip_101', [], ...
    'emgRMS_bipNorm_101', [] ), nReps, 1);


%% --- Extract rep-level CV segments + time-normalize to 101 points ---
x_new = linspace(0,1,nTimeNorm); %creates a vector for a normalized time base 0 - 1 as 100% of old axis

for r = 1:nReps
    stage3.rep(r).is_valid = logical(is_valid_rep(r));
    stage3.rep(r).CV_idx   = CV_idx(r,:);

    if ~stage3.rep(r).is_valid
        continue
    end

    i0 = CV_idx(r,1); i1 = CV_idx(r,2);
    if ~isfinite(i0) || ~isfinite(i1) || i1 <= i0
        continue
    end

    idx = i0:i1; % this the global old time vector
    
    % Scalars (within CV)
    stage3.rep(r).torque_mean_CV = mean(torque_active(idx), 'omitnan'); % ascalr for the entire trial just for QC and states later
    stage3.rep(r).angle_mean_CV  = mean(angle(idx), 'omitnan');

    % EMG mean across bipolars (time series + scalar)
    emgMean  = mean(bip_trial_rms(idx,:),      2, 'omitnan');  % [len(S) x 1c], RMS bp
    emgMeanN = mean(bip_trial_rms_norm(idx,:), 2, 'omitnan');  % [len(S) x 1c], norm RMS bp

    stage3.rep(r).emgRMS_mean_CV         = mean(emgMean,  'omitnan'); % scalar for the emg over CV
    stage3.rep(r).emgRMS_mean_CV_normMVC = mean(emgMeanN, 'omitnan'); % scalar for the Normalized emg

    % Time normalization to 101 points (CV phase)
    x_old = (0:numel(idx)-1) ./ (numel(idx)-1); %x_old becomes: 0 → 1 normalized scale with the same number of points as the CV segment

    stage3.rep(r).torque_101 = interp1(x_old, torque_active(idx), x_new, 'linear').'; %now, bring the torque values to the new time 101 pts
    stage3.rep(r).angle_101  = interp1(x_old, angle(idx),           x_new, 'linear').'; % same for angle

    stage3.rep(r).emgRMSmean_101     = interp1(x_old, emgMean,  x_new, 'linear').'; % same for emgRMS_mean [N x meanRMS]
    stage3.rep(r).emgRMSmeanNorm_101 = interp1(x_old, emgMeanN, x_new, 'linear').';  % same for emgRMS_normalized

    % Per-bipolar time series (28 x 101)
    stage3.rep(r).emgRMS_bip_101     = interp1(x_old, bip_trial_rms(idx,:),      x_new, 'linear').'; % 28x101 bp
    stage3.rep(r).emgRMS_bipNorm_101 = interp1(x_old, bip_trial_rms_norm(idx,:), x_new, 'linear').'; % 28x101 norm bp based
end



%% --- Save Stage 3 output ---
[~, baseName, ~] = fileparts(fileT);
outName = fullfile(pathT, sprintf('%s_stage3_EMGnorm.mat', baseName));
save(outName, 'stage3', 'repQC', '-v7.3');
fprintf('\n✅ Saved Stage 3 EMG normalization output:\n%s\n', outName);

%% ===== Local functions for EMG signal filtering =====
function [bip_rms, bip_bp] = local_bipolar_bp_rms(emg32, bipolar_pairs, b_bp, a_bp, rms_win_samp)

    N  = size(emg32,1);
    nB = size(bipolar_pairs,1);

    % bipolar raw: distal - proximal
    bip_raw = nan(N,nB);
    for k = 1:nB
        ch1 = bipolar_pairs(k,1);
        ch2 = bipolar_pairs(k,2);
        bip_raw(:,k) = emg32(:,ch2) - emg32(:,ch1);
    end

    % bandpass filtering (zero phase)
    bip_bp = nan(size(bip_raw));
    for k = 1:nB
        x = bip_raw(:,k);

        nanMask = isnan(x);
        if any(nanMask)
            x = fillmissing(x,'linear','EndValues','nearest');
        end

        y = filtfilt(b_bp, a_bp, x);

        if any(nanMask)
            y(nanMask) = NaN;
        end

        bip_bp(:,k) = y;
    end

    % RMS envelope
    bip_rms = sqrt(movmean(bip_bp.^2, rms_win_samp, 1, 'omitnan'));

end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end







% QC checking

fprintf('\n===== Stage 3 STRUCT CHECK =====\n')

% number of reps
nReps = numel(stage3.rep);
fprintf('Number of reps: %d\n', nReps);

%% check MVC references
fprintf('\n--- MVC reference check ---\n')

disp(size(stage3.mvc.MVC_ref_used))
disp(size(stage3.mvc.MVC_ref95))
disp(size(stage3.mvc.MVC_refMax))

fprintf('NaN channels (ref95): %d\n', sum(isnan(stage3.mvc.MVC_ref95)));
fprintf('NaN channels (refMax): %d\n', sum(isnan(stage3.mvc.MVC_refMax)));

%% check rep fields

validCount = 0;

for r = 1:nReps

    if stage3.rep(r).is_valid

        validCount = validCount + 1;

        fprintf('\nRep %d\n',r)

        disp(size(stage3.rep(r).torque_101))
        disp(size(stage3.rep(r).angle_101))
        disp(size(stage3.rep(r).emgRMS_bip_101))
        disp(size(stage3.rep(r).emgRMS_bipNorm_101))

    end

end

fprintf('\nValid reps: %d / %d\n', validCount, nReps);

%% visualize EMG map for one rep

r = find([stage3.rep.is_valid],1);

figure
imagesc(stage3.rep(r).emgRMS_bipNorm_101)
colorbar
xlabel('Time (101)')
ylabel('Bipolar channel')
title('Normalized RMS EMG map (%MVC)')