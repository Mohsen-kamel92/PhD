%% ===== Stage 02A: Passive torque correction and CV window indices =====
%
% PURPOSE
%   Remove the passive contribution from the recorded torque, so that later
%   stages work with active torque, and convert the constant velocity windows
%   detected at Stage 1 into sample indices.
%
% SUPPORTS
%   The active torque definition in the Data Analysis section of the
%   manuscript.
%
% INPUT
%   <set folder>/<tag>_sync_rawVolts_stage1_matchQC.mat        (Stage 01)
%
% OUTPUT
%   <set folder>/<P>_<COND>_<Set>_ActiveTorqueCorrected.mat
%   <ROOT>/QC_Stage2/<tag>.png
%   <ROOT>/stage2_summary.csv
%
% METHOD
%   Torque and angle are low-pass filtered, the passive torque-angle relation
%   is fitted across the whole trial (Auto_RestTorque), and the fitted passive
%   component is subtracted to give active torque (Normalize_Torque). The
%   constant velocity windows detected at Stage 01 are carried forward and
%   converted to sample indices for Stage 03, so the window is defined once and
%   used everywhere.
%
% ON THE SAMPLING RATE
%   Fs is read from sync.Fs_emg rather than recomputed from the time vector, so
%   that every stage uses the rate Stage 00 measured.
%
% ON THE PER-REPETITION EMG COPIES
%   Setting saveRepEMGcopy false omits the per-repetition EMG copies. They
%   roughly triple the output file size and are only needed if a later stage
%   reads them directly rather than re-slicing from emg_raw.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped. A file
%   that errors is logged and the batch continues.
%
% DEPENDENCIES
%   Wfilt.m             zero-phase Butterworth wrapper         (see /functions)
%   Auto_RestTorque.m   fits the passive torque-angle relation (see /functions)
%   Normalize_Torque.m  subtracts the fitted passive component (see /functions)

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'CutOff_torque',  20);     % Hz
CFG = set_default(CFG, 'CutOff_angle',   6);      % Hz
CFG = set_default(CFG, 'saveRepEMGcopy', true);   % see header
CFG = set_default(CFG, 'saveFigures',    true);

%% ---------------------------
%  (2) Data root and file list
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 2: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

for dep = {'Wfilt','Auto_RestTorque','Normalize_Torque'}
    assert(~isempty(which(dep{1})), ...
        'Stage 2: %s not found on the MATLAB path. Add the /functions folder.', dep{1});
end

qcDir = fullfile(ROOT,'QC_Stage2');
if CFG.saveFigures && ~exist(qcDir,'dir'), mkdir(qcDir); end

F = dir(fullfile(ROOT, '**', '*_stage1_matchQC.mat'));
inSet = arrayfun(@(x) ~isempty(regexp(x.folder,'[\\/]Set_\d+$','once')), F);
F = F(inSet);
assert(~isempty(F), 'Stage 2: no Stage 1 files found under %s', ROOT);

fprintf('\n===== Stage 2: passive torque correction =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Found %d Stage 1 files.\n\n', numel(F));

%% ---------------------------
%  (3) Batch loop
%% ---------------------------
nF  = numel(F);
res = table('Size',[nF 10], ...
    'VariableTypes',{'string','string','string','double','double','double', ...
                     'double','double','double','string'}, ...
    'VariableNames',{'participant','condition','set','nReps','nWindows', ...
                     'meanActiveTorque_Nm','meanTorque_20to80_Nm', ...
                     'passiveRange_Nm','Fs_Hz','status'});

for k = 1:nF
    fp  = fullfile(F(k).folder, F(k).name);
    tag = "";
    partName = "unknown"; condName = ""; setName = "";

    try
        L = load(fp);
        assert(isfield(L,'sync'),  'No "sync" struct.');
        assert(isfield(L,'repQC'), 'No "repQC" struct.');
        sync  = L.sync;
        repQC = L.repQC;

        %% --- labels ---
        [condPath, setName]  = fileparts(F(k).folder);
        [~,        condName] = fileparts(condPath);

        % Stage 0 records the names it used, which take precedence over the
        % folder names in case a folder was renamed after synchronization.
        if isfield(sync,'meta')
            if isfield(sync.meta,'participantName'), partName = char(sync.meta.participantName); end
            if isfield(sync.meta,'conditionName'),   condName = char(sync.meta.conditionName);   end
            if isfield(sync.meta,'setName'),         setName  = char(sync.meta.setName);         end
        end

        participantName = string(partName);
        conditionName   = string(condName);
        setName         = string(setName);
        tag = participantName + "_" + conditionName + "_" + setName;

        %% --- core signals ---
        time       = sync.time(:);
        torque_raw = sync.torque_raw(:);
        angle_raw  = sync.angle(:);

        assert(isfield(sync,'emg') && ~isempty(sync.emg), ...
            'sync.emg missing. Stage 3 requires emg_raw.');
        emg_raw = sync.emg;

        if isfield(sync,'AUX'),  AUX  = sync.AUX(:);  else, AUX  = []; end
        if isfield(sync,'dac1'), dac1 = sync.dac1(:); else, dac1 = []; end

        if isfield(sync,'Fs_emg') && isfinite(sync.Fs_emg)
            Fs = sync.Fs_emg;
        else
            Fs = 1 / median(diff(time));
        end

        %% --- validate repQC ---
        assert(isfield(repQC,'repIdx') && ~isempty(repQC.repIdx), 'repQC.repIdx missing.');
        repIdx = repQC.repIdx;
        nReps  = size(repIdx,1);
        assert(size(repIdx,2) == 2, 'repQC.repIdx must be [nReps x 2].');
        used_rep_idx = (1:nReps).';

        %% --- filter and subtract passive torque ---
        torque_filt = real(Wfilt(torque_raw, CFG.CutOff_torque, 'low', Fs));
        angle_filt  = real(Wfilt(angle_raw,  CFG.CutOff_angle,  'low', Fs));

        passive_p = Auto_RestTorque(torque_filt, angle_filt, Fs);
        [tau_pass_hat, tau_active] = Normalize_Torque(torque_filt, angle_filt, passive_p);

        torque_active  = tau_active;
        angle          = angle_filt;
        angle_velocity = [0; diff(angle) ./ diff(time)];

        %% --- carry the Stage 1 CV windows forward as indices ---
        if isfield(repQC,'CVwin_globalTime') && ~isempty(repQC.CVwin_globalTime)
            CVwin_globalTime = repQC.CVwin_globalTime;
        elseif isfield(repQC,'CVwin_time') && ~isempty(repQC.CVwin_time)
            CVwin_globalTime = nan(nReps,2);
            for r = 1:nReps
                CVwin_globalTime(r,:) = time(repIdx(r,1)) + repQC.CVwin_time(r,:);
            end
        else
            error('No Stage 1 CV window found in repQC.');
        end

        CVwin_idx  = nan(nReps,2);
        CVwin_time = nan(nReps,2);

        for r = 1:nReps
            if all(isfinite(CVwin_globalTime(r,:)))
                CVwin_idx(r,1)  = nearest_idx(time, CVwin_globalTime(r,1));
                CVwin_idx(r,2)  = nearest_idx(time, CVwin_globalTime(r,2));
                CVwin_time(r,:) = time(CVwin_idx(r,:)) - time(repIdx(r,1));
            end
        end

        % activeWin_* holds the same values as CVwin_* under an older name.
        % Kept so that a script written against either name still runs.
        activeWin_idx        = CVwin_idx;
        activeWin_globalTime = CVwin_globalTime;
        activeWin_time       = CVwin_time;

        repQC.activeWin_idx        = activeWin_idx;
        repQC.activeWin_globalTime = activeWin_globalTime;
        repQC.activeWin_time       = activeWin_time;
        repQC.CVwin_idx            = CVwin_idx;
        repQC.CVwin_globalTime     = CVwin_globalTime;
        repQC.CVwin_time           = CVwin_time;

        repQC.detector.method = 'CV window carried forward from Stage 1';
        repQC.detector.source = 'repQC.CVwin_globalTime';

        nWindows = nnz(all(isfinite(CVwin_idx),2));

        %% --- per-repetition storage ---
        time_CV_raw_stage2          = cell(nReps,1);
        angle_CV_raw_stage2         = cell(nReps,1);
        torque_active_CV_raw_stage2 = cell(nReps,1);
        emg_CV_raw_stage2           = cell(nReps,1);

        time_matched   = cell(nReps,1); angle_matched       = cell(nReps,1);
        torque_matched = cell(nReps,1); emg_raw_matched     = cell(nReps,1);
        torque_raw_matched = cell(nReps,1); torque_filt_matched = cell(nReps,1);
        angle_raw_matched  = cell(nReps,1); angle_filt_matched  = cell(nReps,1);

        for r = 1:nReps
            if all(isfinite(CVwin_idx(r,:))) && CVwin_idx(r,2) > CVwin_idx(r,1)
                idx = CVwin_idx(r,1):CVwin_idx(r,2);
                time_CV_raw_stage2{r}          = time(idx);
                angle_CV_raw_stage2{r}         = angle(idx);
                torque_active_CV_raw_stage2{r} = torque_active(idx);
                if CFG.saveRepEMGcopy, emg_CV_raw_stage2{r} = emg_raw(idx,:); end
            end

            i0 = repIdx(r,1); i1 = repIdx(r,2);
            if ~isfinite(i0) || ~isfinite(i1) || i1 <= i0 || i0 < 1 || i1 > numel(time)
                continue
            end
            idx = i0:i1;
            time_matched{r}        = time(idx);
            angle_matched{r}       = angle(idx);
            torque_matched{r}      = torque_active(idx);
            torque_raw_matched{r}  = torque_raw(idx);
            torque_filt_matched{r} = torque_filt(idx);
            angle_raw_matched{r}   = angle_raw(idx);
            angle_filt_matched{r}  = angle(idx);
            if CFG.saveRepEMGcopy, emg_raw_matched{r} = emg_raw(idx,:); end
        end

        repQC.time_CV_raw_stage2          = time_CV_raw_stage2;
        repQC.angle_CV_raw_stage2         = angle_CV_raw_stage2;
        repQC.torque_active_CV_raw_stage2 = torque_active_CV_raw_stage2;
        repQC.emg_CV_raw_stage2           = emg_CV_raw_stage2;

        %% --- time-normalized active torque, for the summary only ---
        nPts = 101;
        torque_101 = nan(nReps, nPts);
        for r = 1:nReps
            i0 = CVwin_idx(r,1); i1 = CVwin_idx(r,2);
            if ~isfinite(i0) || ~isfinite(i1) || i1 <= i0, continue; end
            idx   = i0:i1;
            t_seg = time(idx);
            torque_101(r,:) = interp1(t_seg, torque_active(idx), ...
                linspace(t_seg(1), t_seg(end), nPts), 'linear');
        end
        torque_101 = torque_101(all(isfinite(torque_101),2), :);

        if isempty(torque_101)
            meanTau = NaN; meanTauMid = NaN;
        else
            mean_tau   = mean(torque_101,1);
            cvAxis     = linspace(0,100,nPts);
            meanTau    = mean(mean_tau);
            meanTauMid = mean(mean_tau(cvAxis >= 20 & cvAxis <= 80));
        end

        passiveRange = max(tau_pass_hat) - min(tau_pass_hat);

        %% --- provenance ---
        angle_is_filtered  = true;
        torque_is_filtered = true;
        angle_filter_Hz    = CFG.CutOff_angle;
        torque_filter_Hz   = CFG.CutOff_torque;
        filter_method      = "Wfilt low-pass";
        method = ['trial passive fit (Auto_RestTorque and Normalize_Torque) ' ...
                  'with the CV window carried forward from Stage 1'];
        sourceFile = fp;
        createdOn  = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
        config     = CFG;

        %% --- QC figure ---
        if CFG.saveFigures
            vthr = 14;
            if isfield(repQC,'params') && isfield(repQC.params,'velocity_thresh')
                vthr = repQC.params.velocity_thresh;
            end
            plot_stage2(time, angle, angle_velocity, torque_filt, torque_active, ...
                repIdx, CVwin_idx, vthr, tag, qcDir);
        end

        %% --- save ---
        outName = fullfile(F(k).folder, sprintf('%s_ActiveTorqueCorrected.mat', ...
            regexprep(tag,'[^\w-]','_')));

        save(outName, ...
            'participantName','conditionName','setName','method','sourceFile', ...
            'createdOn','config', ...
            'time','Fs', ...
            'angle_raw','angle_filt','angle','angle_velocity', ...
            'torque_raw','torque_filt','torque_active', ...
            'passive_p','tau_pass_hat', ...
            'emg_raw','AUX','dac1', ...
            'time_matched','angle_matched','torque_matched','emg_raw_matched', ...
            'torque_raw_matched','torque_filt_matched','angle_raw_matched','angle_filt_matched', ...
            'repQC','used_rep_idx', ...
            'angle_is_filtered','torque_is_filtered','angle_filter_Hz', ...
            'torque_filter_Hz','filter_method', ...
            '-v7.3');

        res(k,:) = {participantName, conditionName, setName, nReps, nWindows, ...
                    meanTau, meanTauMid, passiveRange, Fs, "ok"};

        fprintf(['[%3d/%3d] %-26s %d/%d windows | active torque %6.2f Nm ' ...
                 '(20-80%%: %6.2f) | passive range %5.2f Nm\n'], ...
            k, nF, tag, nWindows, nReps, meanTau, meanTauMid, passiveRange);

    catch ME
        res(k,:) = {string(partName), string(condName), string(setName), ...
                    NaN, NaN, NaN, NaN, NaN, NaN, string(ME.message)};
        fprintf('[%3d/%3d] %-26s FAILED: %s\n', k, nF, tag, ME.message);
    end
end

%% ---------------------------
%  (4) Summary
%% ---------------------------
ok = res.status == "ok";
fprintf('\n===== STAGE 2 SUMMARY =====\n');
fprintf('Processed %d | ok %d | failed %d\n', nF, nnz(ok), nnz(~ok));

if any(ok)
    miss = res(ok & res.nWindows < res.nReps, :);
    if isempty(miss)
        fprintf('All repetitions have a valid CV window.\n');
    else
        fprintf('\nFiles with repetitions missing a CV window:\n');
        disp(miss(:, {'participant','condition','set','nWindows','nReps'}));
    end

    fprintf('\nMean active torque over the 20-80%% CV window, by condition:\n');
    for c = ["CON_75","CON_90","ECC_75","ECC_90"]
        sel = ok & res.condition == c;
        if any(sel)
            fprintf('  %-7s %6.2f +/- %5.2f Nm  (n = %d)\n', c, ...
                mean(res.meanTorque_20to80_Nm(sel)), ...
                std(res.meanTorque_20to80_Nm(sel)), nnz(sel));
        end
    end

    fprintf('\nPassive torque range: median %.2f Nm (max %.2f Nm)\n', ...
        median(res.passiveRange_Nm(ok)), max(res.passiveRange_Nm(ok)));

    if numel(unique(round(res.Fs_Hz(ok),3))) > 1
        fprintf('WARNING: more than one sampling rate across files.\n');
    end
end

if any(~ok)
    fprintf('\nFailures:\n');
    disp(res(~ok, {'participant','condition','set','status'}));
end

writetable(res, fullfile(ROOT,'stage2_summary.csv'));
fprintf('\nSaved summary: %s\n', fullfile(ROOT,'stage2_summary.csv'));
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

function idx = nearest_idx(t, q)
% Index of the sample in t closest to the scalar time q.
    [~, idx] = min(abs(t - q));
end

function plot_stage2(time, angle, vel, torque_filt, torque_active, ...
                     repIdx, CVwin_idx, vthr, tag, qcDir)
% Joint angle, angular velocity and torque before and after passive
% subtraction, with repetition onsets and CV windows marked.

    f = figure('Color','w','Visible','off','Position',[100 100 1200 800]);
    tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');

    a1 = nexttile; hold(a1,'on');
    plot(a1, time, angle, 'Color',[0 .55 0], 'LineWidth',1.5);
    ylabel(a1,'Angle (deg)','FontWeight','bold');
    set(a1,'XTickLabel',[]);
    title(a1, tag + " | passive subtraction and CV windows", ...
        'Interpreter','none','FontWeight','bold');

    a2 = nexttile; hold(a2,'on');
    plot(a2, time, vel, 'Color',[.25 .25 .25], 'LineWidth',1.2);
    yline(a2,  vthr, 'm--'); yline(a2, -vthr, 'm--');
    ylabel(a2,'Vel (deg/s)','FontWeight','bold');
    set(a2,'XTickLabel',[]);

    a3 = nexttile; hold(a3,'on');
    hRaw = plot(a3, time, torque_filt,   'LineWidth',1.0, 'Color',[.6 .6 .6]);
    hNet = plot(a3, time, torque_active, 'LineWidth',1.5, 'Color',[.5 0 .5]);
    ylabel(a3,'Torque (Nm)','FontWeight','bold');
    xlabel(a3,'Time (s)','FontWeight','bold');

    for r = 1:size(repIdx,1)
        if all(isfinite(repIdx(r,:)))
            xline(a2, time(repIdx(r,1)), 'k--','LineWidth',0.8);
            xline(a3, time(repIdx(r,1)), 'k--','LineWidth',0.8);
        end
        if all(isfinite(CVwin_idx(r,:))) && CVwin_idx(r,2) > CVwin_idx(r,1)
            xline(a2, time(CVwin_idx(r,1)), 'r--','LineWidth',1.2);
            xline(a2, time(CVwin_idx(r,2)), 'r--','LineWidth',1.2);
            xline(a3, time(CVwin_idx(r,1)), 'r--','LineWidth',1.2);
            xline(a3, time(CVwin_idx(r,2)), 'r--','LineWidth',1.2);
        end
    end

    legend(a3, [hRaw hNet], {'Filtered torque','Active torque'}, ...
        'Location','best','Box','off');

    linkaxes([a1 a2 a3],'x');
    set([a1 a2 a3],'Box','off','TickDir','out','FontWeight','bold');

    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + ".png"), ...
        'Resolution',110);
    close(f);
end
