%% ===== Stage 03: Bipolar EMG, RMS envelope, MVC normalization and CV extraction =====
%
% PURPOSE
%   Convert the monopolar grid to bipolar derivations, band-pass filter and
%   RMS-envelope them, express amplitude relative to the participant's MVA, and
%   extract each repetition over its constant velocity window, time-normalized
%   to a common length.
%
% SUPPORTS
%   The EMG processing and normalization described in the Data Analysis section
%   of the manuscript.
%
% INPUT
%   <ROOT>/<P>/<COND>/<Set_N>/<tag>_ActiveTorqueCorrected.mat      (Stage 02A)
%   <ROOT>/<P>/MVC/<tag>_MVC_ActiveTorqueCorrected.mat             (Stage 02B MVC)
%
% OUTPUT
%   <set folder>/<tag>_stage3_EMGnorm.mat
%   <ROOT>/QC_Stage3/<tag>.png
%   <ROOT>/stage3_summary.csv
%
% NORMALIZATION
%   Amplitudes are expressed relative to a single scalar reference per
%   participant, the 95th percentile of the grid-averaged MVC RMS.
%
%   This differs from normalizing each bipolar channel to its own MVC value.
%   The centroid and entropy calculations downstream both depend on the
%   relative amplitudes between channels. A centroid is invariant to a global
%   scale factor but not to a per-channel one, so dividing each channel by its
%   own reference re-weights the spatial map by the reciprocal of the MVC map
%   and displaces the centroid. Using one scalar leaves the spatial pattern
%   untouched while keeping the %MVA units, and matches the convention in the
%   HDsEMG spatial distribution literature, where the mean RMS across channels
%   is normalized to the mean RMS-EMG %MVA.
%
%   Per-channel references are still computed and saved for transparency, but
%   are not used in the normalization.
%
%   A percentile is used rather than the peak because the MVC was dynamic and a
%   single-sample maximum is easily set by a transient. The margin between the
%   two is reported per participant so the choice can be judged.
%
% METHOD
%   32 monopolar channels become 28 bipolar derivations along the proximal to
%   distal axis, using the validated channel map below. Signals are band-pass
%   filtered between 20 and 450 Hz with a fourth-order Butterworth applied in
%   both directions, giving a zero-phase response. The RMS envelope is computed
%   over a centred 50 ms window. Constant velocity windows are reused from
%   repQC, detected at Stage 1 and converted to indices at Stage 2. Each
%   repetition is then time-normalized to 101 points.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped. A file
%   that errors is logged and the batch continues.
%
% DEPENDENCIES
%   Signal Processing Toolbox   butter, filtfilt
%   Statistics and Machine Learning Toolbox   prctile

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'bp_low',      20);    % Hz
CFG = set_default(CFG, 'bp_high',     450);   % Hz
CFG = set_default(CFG, 'bp_order',    4);     % applied in both directions
CFG = set_default(CFG, 'rms_win_ms',  50);    % centred RMS window
CFG = set_default(CFG, 'mvc_prc',     95);    % percentile of the grid-averaged MVC RMS
CFG = set_default(CFG, 'nTimeNorm',   101);
CFG = set_default(CFG, 'saveFigures', true);

%% ---------------------------
%  (2) Channel map
%% ---------------------------
% Row 1 is proximal and row 8 distal; column 1 is lateral and column 4 medial.

ChMap = [1 16 17 32;
         2 15 18 31;
         3 14 19 30;
         4 13 20 29;
         5 12 21 28;
         6 11 22 27;
         7 10 23 26;
         8  9 24 25];

[nRows, nCols] = size(ChMap);
nBipolar = (nRows-1) * nCols;

% Derivation is distal minus proximal. The sign does not affect the RMS, but
% the convention is recorded so the map can be interpreted.
bipolar_direction = "distal_minus_proximal";

bipolar_pairs  = nan(nBipolar,2);
bipolar_labels = strings(nBipolar,1);
bipolar_rowcol = nan(nBipolar,2);

kk = 1;
for col = 1:nCols
    for row = 1:(nRows-1)
        bipolar_pairs(kk,:)  = [ChMap(row,col) ChMap(row+1,col)];
        bipolar_labels(kk)   = sprintf('R%dC%d', row, col);
        bipolar_rowcol(kk,:) = [row col];
        kk = kk + 1;
    end
end

%% ---------------------------
%  (3) Data root and file list
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 3: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

for dep = {'butter','filtfilt','prctile'}
    assert(~isempty(which(dep{1})), ...
        'Stage 3: %s not available. See the toolbox dependencies in the header.', dep{1});
end

qcDir = fullfile(ROOT,'QC_Stage3');
if CFG.saveFigures && ~exist(qcDir,'dir'), mkdir(qcDir); end

F = dir(fullfile(ROOT, '**', '*_ActiveTorqueCorrected.mat'));
F = F(~contains({F.name},'MVC_ActiveTorqueCorrected'));
inSet = arrayfun(@(x) ~isempty(regexp(x.folder,'[\\/]Set_\d+$','once')), F);
F = F(inSet);
assert(~isempty(F), 'Stage 3: no Stage 2 trial files found under %s', ROOT);

fprintf('\n===== Stage 3: bipolar EMG and MVC normalization =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Found %d trial files.\n\n', numel(F));

%% ---------------------------
%  (4) Batch loop
%% ---------------------------
nF  = numel(F);
res = table('Size',[nF 12], ...
    'VariableTypes',{'string','string','string','double','double','double', ...
                     'double','double','double','double','double','string'}, ...
    'VariableNames',{'participant','condition','set','nReps','nValid', ...
                     'MVC_ref_global','MVC_maxOverP95_pct','meanPctMVA', ...
                     'maxPctMVA','pctSamplesOver100','ratio95max_median','status'});

for k = 1:nF
    fp  = fullfile(F(k).folder, F(k).name);
    tag = "";
    partName = ''; condName = ''; setName = '';

    try
        %% --- load trial ---
        T = load(fp);
        for fld = {'repQC','emg_raw','time','angle','torque_active'}
            assert(isfield(T, fld{1}), 'Trial file missing %s.', fld{1});
        end

        repQC         = T.repQC;
        time          = T.time(:);
        angle         = T.angle(:);
        torque_active = T.torque_active(:);
        emg_raw       = T.emg_raw;

        [condPath, setName]  = fileparts(F(k).folder);
        [partPath, condName] = fileparts(condPath);
        [~,        partName] = fileparts(partPath);
        if isfield(T,'participantName'), partName = char(T.participantName); end
        if isfield(T,'conditionName'),   condName = char(T.conditionName);   end
        if isfield(T,'setName'),         setName  = char(T.setName);         end
        tag = string(sprintf('%s_%s_%s', partName, condName, setName));

        if isfield(T,'Fs') && isfinite(T.Fs)
            Fs = T.Fs;
        else
            Fs = 1 / median(diff(time),'omitnan');
        end

        %% --- locate and load this participant's MVC ---
        mvcDir = fullfile(partPath,'MVC');
        assert(isfolder(mvcDir), 'No MVC folder at %s', mvcDir);

        Mf = dir(fullfile(mvcDir,'*_MVC_ActiveTorqueCorrected.mat'));
        assert(isscalar(Mf), 'Expected 1 MVC file in %s, found %d.', mvcDir, numel(Mf));
        M = load(fullfile(Mf.folder, Mf.name));
        assert(isfield(M,'emg_raw_matched') && ~isempty(M.emg_raw_matched), ...
            'MVC file missing emg_raw_matched.');

        % Guards against pairing a trial with the wrong participant's reference
        if isfield(M,'participantName')
            assert(strcmp(char(M.participantName), partName), ...
                'MVC participant "%s" does not match trial participant "%s".', ...
                char(M.participantName), partName);
        end

        if isfield(M,'Fs') && isfinite(M.Fs)
            assert(abs(M.Fs - Fs) < 1e-6, ...
                'Fs mismatch: trial %.6f Hz, MVC %.6f Hz.', Fs, M.Fs);
        end

        emg_mvc_raw = M.emg_raw_matched;

        %% --- filter design ---
        fnyq = Fs/2;
        assert(CFG.bp_high < fnyq, 'bp_high must be below Nyquist.');
        [b_bp, a_bp] = butter(CFG.bp_order, [CFG.bp_low CFG.bp_high]/fnyq, 'bandpass');

        rms_win_samp = max(1, round((CFG.rms_win_ms/1000)*Fs));
        if mod(rms_win_samp,2)==0, rms_win_samp = rms_win_samp + 1; end

        %% --- MVC reference ---
        bip_mvc_rms = bipolar_bp_rms(emg_mvc_raw, bipolar_pairs, ...
                                     b_bp, a_bp, rms_win_samp);

        % Primary reference: one scalar, from the grid-averaged MVC RMS.
        mvc_grid_mean  = mean(bip_mvc_rms, 2, 'omitnan');
        MVC_ref_global = prctile(mvc_grid_mean, CFG.mvc_prc);
        MVC_max_global = max(mvc_grid_mean, [], 'omitnan');
        assert(isfinite(MVC_ref_global) && MVC_ref_global > 0, ...
            'Invalid global MVC reference.');

        % By how much a single-sample maximum would have exceeded the
        % percentile, which is the margin the percentile choice avoids.
        maxOverP95_pct = 100 * (MVC_max_global - MVC_ref_global) / MVC_ref_global;

        % Per-channel references, saved for transparency but not used.
        MVC_ref95  = prctile(bip_mvc_rms, CFG.mvc_prc, 1);
        MVC_refMax = max(bip_mvc_rms, [], 1, 'omitnan');
        MVC_ref95(MVC_ref95   <= 0 | ~isfinite(MVC_ref95))  = NaN;
        MVC_refMax(MVC_refMax <= 0 | ~isfinite(MVC_refMax)) = NaN;
        ratio95 = MVC_ref95 ./ MVC_refMax;

        %% --- CV windows ---
        assert(isfield(repQC,'repIdx') && ~isempty(repQC.repIdx), 'repQC.repIdx missing.');
        nReps = size(repQC.repIdx,1);
        assert(isfield(repQC,'CVwin_idx') && ~isempty(repQC.CVwin_idx), ...
            'repQC.CVwin_idx missing. Re-run Stage 2.');
        CV_idx = repQC.CVwin_idx;

        % Every repetition is carried forward. Sets were accepted or repeated
        % during testing, so there is no repetition-level exclusion rule here.
        is_used_rep = true(nReps,1);

        %% --- trial EMG, filtered once for the whole recording ---
        bip_trial_rms      = bipolar_bp_rms(emg_raw, bipolar_pairs, ...
                                            b_bp, a_bp, rms_win_samp);
        bip_trial_rms_norm = (bip_trial_rms ./ MVC_ref_global) * 100;

        %% --- output struct ---
        stage3 = struct();

        stage3.meta = struct( ...
            'participantName', string(partName), ...
            'conditionName',   string(condName), ...
            'setName',         string(setName), ...
            'trialFile',       string(fp), ...
            'mvcFile',         string(fullfile(Mf.folder, Mf.name)), ...
            'createdOn',       string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'script',          string(mfilename), ...
            'matlab',          string(version));

        stage3.params = struct( ...
            'Fs', Fs, 'bp_low', CFG.bp_low, 'bp_high', CFG.bp_high, ...
            'bp_order', CFG.bp_order, ...
            'filter_note', "butter applied with filtfilt, zero-phase response", ...
            'rms_win_ms', CFG.rms_win_ms, 'rms_win_samp', rms_win_samp, ...
            'mvc_prc', CFG.mvc_prc, ...
            'norm_scope', "global", ...
            'norm_note', ["%MVA computed with a single scalar reference (95th " ...
                          "percentile of the grid-averaged MVC RMS) so that centroid " ...
                          "and entropy are invariant to the normalization"], ...
            'timeNorm_points', CFG.nTimeNorm, ...
            'timeNorm_axis', linspace(0,100,CFG.nTimeNorm), ...
            'timeNorm_units', "percent of the constant velocity window", ...
            'interp_method', "linear", ...
            'torque_units', "Nm, passive component removed", ...
            'angle_units', "deg, filtered at Stage 2", ...
            'rep_inclusion_logic', "all repetitions used; matching QC is set-level");

        stage3.map = struct('ChMap', ChMap, ...
            'bipolar_pairs', bipolar_pairs, 'bipolar_labels', bipolar_labels, ...
            'bipolar_rowcol', bipolar_rowcol, 'bipolar_direction', bipolar_direction, ...
            'proximal_row', 1, 'distal_row', nRows, ...
            'lateral_column', 1, 'medial_column', nCols, ...
            'row_meaning', "1 = proximal, 8 = distal", ...
            'column_meaning', "1 = lateral, 4 = medial");

        stage3.mvc = struct( ...
            'MVC_ref_global', MVC_ref_global, ...
            'MVC_ref_used',   MVC_ref_global, ...
            'MVC_max_global', MVC_max_global, ...
            'maxOverP95_pct', maxOverP95_pct, ...
            'MVC_ref95',      MVC_ref95, ...
            'MVC_refMax',     MVC_refMax, ...
            'ratio_ref95_refMax', ratio95, ...
            'bip_mvc_rms',    bip_mvc_rms);

        stage3.repIdx = repQC.repIdx;

        stage3.rep = repmat(struct('is_valid',false,'CV_idx',[NaN NaN], ...
            'torque_mean_CV',NaN,'angle_mean_CV',NaN, ...
            'emgRMS_mean_CV',NaN,'emgRMS_mean_CV_normMVC',NaN, ...
            'torque_101',[],'angle_101',[], ...
            'emgRMSmean_101',[],'emgRMSmeanNorm_101',[], ...
            'emgRMS_bip_101',[],'emgRMS_bipNorm_101',[]), nReps, 1);

        %% --- per-repetition extraction and time normalization ---
        x_new = linspace(0,1,CFG.nTimeNorm);

        for r = 1:nReps
            stage3.rep(r).is_valid = logical(is_used_rep(r));
            stage3.rep(r).CV_idx   = CV_idx(r,:);
            if ~stage3.rep(r).is_valid, continue; end

            i0 = CV_idx(r,1); i1 = CV_idx(r,2);
            if ~isfinite(i0) || ~isfinite(i1) || i1 <= i0, continue; end
            idx = i0:i1;

            stage3.rep(r).torque_mean_CV = mean(torque_active(idx),'omitnan');
            stage3.rep(r).angle_mean_CV  = mean(angle(idx),'omitnan');

            emgMean  = mean(bip_trial_rms(idx,:),      2, 'omitnan');
            emgMeanN = mean(bip_trial_rms_norm(idx,:), 2, 'omitnan');

            stage3.rep(r).emgRMS_mean_CV         = mean(emgMean,'omitnan');
            stage3.rep(r).emgRMS_mean_CV_normMVC = mean(emgMeanN,'omitnan');

            x_old = (0:numel(idx)-1) ./ (numel(idx)-1);

            stage3.rep(r).torque_101 = interp1(x_old, torque_active(idx), x_new,'linear').';
            stage3.rep(r).angle_101  = interp1(x_old, angle(idx),         x_new,'linear').';
            stage3.rep(r).emgRMSmean_101     = interp1(x_old, emgMean,  x_new,'linear').';
            stage3.rep(r).emgRMSmeanNorm_101 = interp1(x_old, emgMeanN, x_new,'linear').';
            stage3.rep(r).emgRMS_bip_101     = interp1(x_old, bip_trial_rms(idx,:),      x_new,'linear').';
            stage3.rep(r).emgRMS_bipNorm_101 = interp1(x_old, bip_trial_rms_norm(idx,:), x_new,'linear').';
        end

        %% --- repetition by time matrices, ready for Stage 6A ---
        torque_101_all         = nan(nReps,CFG.nTimeNorm);
        angle_101_all          = nan(nReps,CFG.nTimeNorm);
        emgRMSmean_101_all     = nan(nReps,CFG.nTimeNorm);
        emgRMSmeanNorm_101_all = nan(nReps,CFG.nTimeNorm);

        for r = 1:nReps
            if stage3.rep(r).is_valid && ~isempty(stage3.rep(r).angle_101)
                torque_101_all(r,:)         = stage3.rep(r).torque_101(:).';
                angle_101_all(r,:)          = stage3.rep(r).angle_101(:).';
                emgRMSmean_101_all(r,:)     = stage3.rep(r).emgRMSmean_101(:).';
                emgRMSmeanNorm_101_all(r,:) = stage3.rep(r).emgRMSmeanNorm_101(:).';
            end
        end

        stage3.groupReady = struct( ...
            'torque_101_all', torque_101_all, ...
            'angle_101_all',  angle_101_all, ...
            'meanAngle_101',  mean(angle_101_all,1,'omitnan'), ...
            'emgRMSmean_101_all',     emgRMSmean_101_all, ...
            'emgRMSmeanNorm_101_all', emgRMSmeanNorm_101_all, ...
            'x_norm_101', linspace(0,100,CFG.nTimeNorm));

        %% --- QC values ---
        nValid = nnz(arrayfun(@(s) s.is_valid && ~isempty(s.angle_101), stage3.rep));

        allNorm = cat(3, stage3.rep(arrayfun(@(s) ~isempty(s.emgRMS_bipNorm_101), ...
                                             stage3.rep)).emgRMS_bipNorm_101);
        if isempty(allNorm)
            meanPct = NaN; maxPct = NaN; over100 = NaN;
        else
            meanPct = mean(allNorm(:),'omitnan');
            maxPct  = max(allNorm(:),[],'omitnan');
            over100 = 100 * nnz(allNorm > 100) / numel(allNorm);
        end

        %% --- QC figure ---
        if CFG.saveFigures
            plot_stage3(stage3, bipolar_rowcol, tag, qcDir);
        end

        %% --- save ---
        [~, baseName] = fileparts(F(k).name);
        outName = fullfile(F(k).folder, sprintf('%s_stage3_EMGnorm.mat', baseName));
        save(outName, 'stage3', 'repQC', '-v7.3');

        res(k,:) = {string(partName), string(condName), string(setName), ...
                    nReps, nValid, MVC_ref_global, maxOverP95_pct, meanPct, ...
                    maxPct, over100, median(ratio95,'omitnan'), "ok"};

        fprintf(['[%3d/%3d] %-26s %d/%d reps | ref %.3g | mean %5.1f %%MVA | ' ...
                 'max %6.1f | >100%%: %4.2f%%\n'], ...
            k, nF, tag, nValid, nReps, MVC_ref_global, meanPct, maxPct, over100);

    catch ME
        res(k,:) = {string(partName), string(condName), string(setName), ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, string(ME.message)};
        fprintf('[%3d/%3d] %-26s FAILED: %s\n', k, nF, tag, ME.message);
    end
end

%% ---------------------------
%  (5) Summary
%% ---------------------------
ok = res.status == "ok";
fprintf('\n===== STAGE 3 SUMMARY =====\n');
fprintf('Processed %d | ok %d | failed %d\n', nF, nnz(ok), nnz(~ok));

if any(ok)
    fprintf('\nMean %%MVA by condition:\n');
    for c = ["CON_75","CON_90","ECC_75","ECC_90"]
        sel = ok & res.condition == c;
        if any(sel)
            fprintf('  %-7s %5.1f +/- %4.1f %%MVA  (n = %d)\n', c, ...
                mean(res.meanPctMVA(sel)), std(res.meanPctMVA(sel)), nnz(sel));
        end
    end

    fprintf('\nSamples above 100%% MVA: median %.2f%% of samples (max %.2f%%)\n', ...
        median(res.pctSamplesOver100(ok)), max(res.pctSamplesOver100(ok)));

    % The margin the percentile reference avoids, reported per participant in
    % the manuscript.
    P = unique(res.participant(ok));
    m = nan(numel(P),1);
    for i = 1:numel(P)
        sel = ok & res.participant == P(i);
        m(i) = median(res.MVC_maxOverP95_pct(sel));
    end
    fprintf('MVC maximum exceeded the %dth percentile by %.0f to %.0f%% across participants.\n', ...
        CFG.mvc_prc, min(m), max(m));

    % One reference per participant, so every file of theirs must agree
    for i = 1:numel(P)
        sel = ok & res.participant == P(i);
        if numel(unique(round(res.MVC_ref_global(sel),12))) > 1
            fprintf('WARNING: participant %s has more than one MVC reference.\n', P(i));
        end
    end

    miss = res(ok & res.nValid < res.nReps, :);
    if ~isempty(miss)
        fprintf('\nFiles with repetitions not extracted:\n');
        disp(miss(:, {'participant','condition','set','nValid','nReps'}));
    end
end

if any(~ok)
    fprintf('\nFailures:\n');
    disp(res(~ok, {'participant','condition','set','status'}));
end

writetable(res, fullfile(ROOT,'stage3_summary.csv'));
fprintf('\nSaved summary: %s\n', fullfile(ROOT,'stage3_summary.csv'));
if CFG.saveFigures, fprintf('QC figures   : %s\n', qcDir); end

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
% selection fails immediately with a useful message rather than partway through
% a batch.
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

        % A repository containing only the processed group file is a valid root for
    % the Stage 6 analyses, which read nothing else.
    hasGroupFile = isfile(fullfile(ROOT,'groupData_stage6A.mat'));

    assert(hasParticipant || hasGroupFile, ...
        ['The selected folder does not look like the expected data root.\n' ...
         'Expected either the raw layout\n' ...
         '  <ROOT>/<participant>/<CON_75|CON_90|ECC_75|ECC_90>/<Set_N>/\n' ...
         'or a folder containing groupData_stage6A.mat\n' ...
         'Selected: %s'], ROOT);
end

function bip_rms = bipolar_bp_rms(emg32, pairs, b_bp, a_bp, win)
% Monopolar grid to bipolar derivations, band-pass filtered in both directions,
% then RMS-enveloped over a centred window. Gaps are interpolated before
% filtering and restored afterwards, so a short dropout cannot spread through
% the filter into neighbouring samples.
    N  = size(emg32,1);
    nB = size(pairs,1);

    bip_raw = nan(N,nB);
    for k = 1:nB
        bip_raw(:,k) = emg32(:,pairs(k,2)) - emg32(:,pairs(k,1));
    end

    bip_bp = nan(size(bip_raw));
    for k = 1:nB
        x = bip_raw(:,k);
        nanMask = isnan(x);
        if any(nanMask)
            x = fillmissing(x,'linear','EndValues','nearest');
        end
        y = filtfilt(b_bp, a_bp, x);
        if any(nanMask), y(nanMask) = NaN; end
        bip_bp(:,k) = y;
    end

    bip_rms = sqrt(movmean(bip_bp.^2, win, 1, 'omitnan'));
end

function plot_stage3(stage3, rowcol, tag, qcDir)
% Three views of the first valid repetition: the channel by time map, the
% time-averaged spatial map, and the grid-mean amplitude across repetitions.
    r = find(arrayfun(@(s) s.is_valid && ~isempty(s.emgRMS_bipNorm_101), stage3.rep), 1);
    if isempty(r), return; end

    f = figure('Color','w','Visible','off','Position',[100 100 1100 700]);
    tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

    a1 = nexttile;
    imagesc(a1, stage3.rep(r).emgRMS_bipNorm_101);
    colorbar(a1); xlabel(a1,'%CV'); ylabel(a1,'Bipolar channel');
    title(a1, sprintf('Rep %d, normalized RMS (%%MVA)', r), 'FontWeight','bold');

    a2 = nexttile;
    m = mean(stage3.rep(r).emgRMS_bipNorm_101, 2, 'omitnan');
    map = nan(max(rowcol(:,1)), max(rowcol(:,2)));
    for c = 1:numel(m), map(rowcol(c,1), rowcol(c,2)) = m(c); end
    imagesc(a2, map); colorbar(a2); axis(a2,'image');
    xlabel(a2,'Column (1 = lateral)'); ylabel(a2,'Row (1 = proximal)');
    title(a2,'Time-averaged spatial map','FontWeight','bold');

    a3 = nexttile([1 2]); hold(a3,'on');
    x = stage3.groupReady.x_norm_101;
    plot(a3, x, stage3.groupReady.emgRMSmeanNorm_101_all.', 'Color',[.7 .7 .7]);
    plot(a3, x, mean(stage3.groupReady.emgRMSmeanNorm_101_all,1,'omitnan'), ...
        'k','LineWidth',2);
    xlabel(a3,'%CV','FontWeight','bold'); ylabel(a3,'%MVA','FontWeight','bold');
    title(a3,'Grid-mean amplitude, all repetitions','FontWeight','bold');
    set(a3,'Box','off','TickDir','out');

    sgtitle(f, tag, 'FontWeight','bold','Interpreter','none');
    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + ".png"), ...
        'Resolution',110);
    close(f);
end
