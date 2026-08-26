%% ===== Stage 02B (MVC): Passive torque correction for the reference contractions =====
%
% PURPOSE
%   Apply the same passive torque subtraction used for the trials to the
%   maximal voluntary contractions, and locate the contraction within each
%   recording.
%
% SUPPORTS
%   The reference values described in the Experimental Protocol section of the
%   manuscript. The MVC defines the peak active torque from which the 75% and
%   90% targets were derived, and the EMG amplitude used to express activity as
%   %MVA. Both are reported per file below so they can be checked rather than
%   assumed.
%
% INPUT
%   <ROOT>/<P>/MVC/<P>_MVC_<name>_sync_rawVolts.mat   (Stage 00)
%
%   MVC recordings stop at Stage 00. They have no prescribed ramp, no
%   repetitions and no matching check, so there is no Stage 1 output for them
%   and nothing to carry forward.
%
% OUTPUT
%   <ROOT>/<P>/MVC/<tag>_MVC_ActiveTorqueCorrected.mat
%   <ROOT>/QC_Stage2_MVC/<tag>.png
%   <ROOT>/stage2_mvc_summary.csv
%
% METHOD
%   Torque and angle are low-pass filtered, the passive torque-angle relation
%   is fitted across the recording and subtracted, and the contraction window
%   runs between the first and the last sample where angular velocity is more
%   negative than the threshold. The concentric rotation runs from
%   plantarflexion to dorsiflexion, so velocity is negative under the sign
%   convention used from Stage 1 onward.
%
% ON THE SEGMENT COUNT
%   The window spans the first to the last threshold crossing, so a recording
%   holding more than one attempt would produce a window covering the gap
%   between them. Contiguous runs above threshold separated by more than
%   CFG.minGap_s are counted as separate segments and reported, so any such
%   recording can be found and inspected rather than silently averaged.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped. A file
%   that errors is logged and the batch continues.
%
% DEPENDENCIES
%   Wfilt.m            zero-phase Butterworth wrapper         (see /functions)
%   Auto_RestTorque.m  fits the passive torque-angle relation (see /functions)
%   Normalize_Torque.m subtracts the fitted passive component (see /functions)

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'velocity_thresh', 14);   % deg/s, same as Stage 1
CFG = set_default(CFG, 'CutOff_torque',   20);   % Hz
CFG = set_default(CFG, 'CutOff_angle',    6);    % Hz
CFG = set_default(CFG, 'minGap_s',        0.5);  % gap splitting velocity segments
CFG = set_default(CFG, 'saveFigures',     true);

%% ---------------------------
%  (2) Data root and file list
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 2 (MVC): no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

for dep = {'Wfilt','Auto_RestTorque','Normalize_Torque'}
    assert(~isempty(which(dep{1})), ...
        'Stage 2 (MVC): %s not found on the MATLAB path. Add the /functions folder.', dep{1});
end

qcDir = fullfile(ROOT,'QC_Stage2_MVC');
if CFG.saveFigures && ~exist(qcDir,'dir'), mkdir(qcDir); end

F = dir(fullfile(ROOT, '**', '*_sync_rawVolts.mat'));

% MVC recordings only. Trials sit in Set_* folders and are handled by the trial
% version of this stage.
isMVC = arrayfun(@(x) ~isempty(regexp(x.folder,'[\\/]MVC$','once')), F);
F = F(isMVC);
assert(~isempty(F), 'Stage 2 (MVC): no MVC sync files found under %s', ROOT);

fprintf('\n===== Stage 2 (MVC): passive torque correction =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Found %d MVC recordings.\n\n', numel(F));

%% ---------------------------
%  (3) Batch loop
%% ---------------------------
nF  = numel(F);
res = table('Size',[nF 10], ...
    'VariableTypes',{'string','string','double','double','double','double', ...
                     'double','double','double','string'}, ...
    'VariableNames',{'participant','mvcFile','peakTorque_Nm','meanTorque_Nm', ...
                     'peakAt_pctWindow','windowDur_s','windowROM_deg', ...
                     'nSegments','passiveRange_Nm','status'});

for k = 1:nF
    fp  = fullfile(F(k).folder, F(k).name);
    tag = "";
    participantName = ""; mvcNameStr = "";

    try
        L = load(fp,'sync');
        assert(isfield(L,'sync'), 'No "sync" struct.');
        sync = L.sync;

        %% --- labels ---
        [partPath, ~] = fileparts(F(k).folder);
        [~, partName] = fileparts(partPath);
        mvcName = 'MVC';

        % Stage 0 records the names it used, which take precedence over the
        % folder names in case a folder was renamed after synchronization.
        if isfield(sync,'meta')
            if isfield(sync.meta,'participantName'), partName = char(sync.meta.participantName); end
            if isfield(sync.meta,'setName'),         mvcName  = char(sync.meta.setName);         end
        end
        participantName = string(partName);
        conditionName   = "MVC";
        mvcNameStr      = string(mvcName);
        tag = string(sprintf('%s_MVC_%s', partName, mvcName));

        %% --- core signals ---
        time       = sync.time(:);
        torque_raw = sync.torque_raw(:);
        angle_raw  = sync.angle(:);

        if isfield(sync,'emg'), emg_raw = sync.emg; else, emg_raw = []; end
        assert(~isempty(emg_raw), ...
            'sync.emg missing. Stage 3 needs the MVC EMG to normalize activity.');

        if isfield(sync,'Fs_emg') && isfinite(sync.Fs_emg)
            Fs = sync.Fs_emg;
        else
            Fs = 1 / median(diff(time));
        end

        %% --- filter and subtract the passive torque ---
        torque_filt = real(Wfilt(torque_raw, CFG.CutOff_torque, 'low', Fs));
        angle_filt  = real(Wfilt(angle_raw,  CFG.CutOff_angle,  'low', Fs));

        passive_p = Auto_RestTorque(torque_filt, angle_filt, Fs);
        [tau_pass_hat, tau_active] = Normalize_Torque(torque_filt, angle_filt, passive_p);

        torque_active  = tau_active;
        angle          = angle_filt;
        angle_velocity = [0; diff(angle) ./ diff(time)];

        %% --- contraction window ---
        above = angle_velocity < -CFG.velocity_thresh;
        assert(any(above), 'No samples exceed the velocity threshold.');

        idx0 = find(above, 1, 'first');
        idx1 = find(above, 1, 'last');
        assert(idx1 > idx0, 'Invalid velocity window.');

        % Counted, not merged. The window itself still spans the full range, so
        % a count above one means the recording needs inspecting.
        nSegments = count_segments(above, time, CFG.minGap_s);

        mvc_t0 = time(idx0);
        mvc_t1 = time(idx1);
        idx    = idx0:idx1;

        %% --- reference values ---
        [peakTorque, iPeak] = max(torque_active(idx));
        meanTorque   = mean(torque_active(idx));
        peakAtPct    = (iPeak - 1) / (numel(idx) - 1) * 100;
        windowDur    = mvc_t1 - mvc_t0;
        windowROM    = abs(angle(idx1) - angle(idx0));
        passiveRange = max(tau_pass_hat) - min(tau_pass_hat);

        %% --- matched segments ---
        time_matched    = time(idx);
        angle_matched   = angle(idx);
        torque_matched  = torque_active(idx);
        emg_raw_matched = emg_raw(idx,:);

        % Duplicate names, kept so that downstream code referencing either
        % continues to work. See the same note in the trial version.
        time_MVC_raw          = time_matched;
        angle_MVC_raw         = angle_matched;
        torque_active_MVC_raw = torque_matched;
        emg_MVC_raw           = emg_raw_matched;

        %% --- provenance ---
        angle_is_filtered  = true;
        torque_is_filtered = true;
        angle_filter_Hz    = CFG.CutOff_angle;
        torque_filter_Hz   = CFG.CutOff_torque;
        filter_method      = "Wfilt low-pass";
        velocity_thresh    = CFG.velocity_thresh;
        method = ['MVC passive fit (Auto_RestTorque and Normalize_Torque) ' ...
                  'and velocity-based contraction window'];
        sourceFile   = fp;
        createdOn    = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
        mvcRecording = mvcNameStr;

        %% --- QC figure ---
        if CFG.saveFigures
            plot_mvc(time, angle, angle_velocity, torque_filt, torque_active, ...
                idx0, idx1, CFG.velocity_thresh, peakTorque, tag, qcDir);
        end

        %% --- save ---
        outName = fullfile(F(k).folder, sprintf('%s_MVC_ActiveTorqueCorrected.mat', ...
            regexprep(tag,'[^\w-]','_')));

        save(outName, ...
            'participantName','conditionName','mvcRecording','method','sourceFile','createdOn', ...
            'time','Fs', ...
            'angle_raw','angle_filt','angle','angle_velocity', ...
            'torque_raw','torque_filt','torque_active', ...
            'passive_p','tau_pass_hat', ...
            'emg_raw', ...
            'time_matched','angle_matched','torque_matched','emg_raw_matched', ...
            'time_MVC_raw','angle_MVC_raw','torque_active_MVC_raw','emg_MVC_raw', ...
            'angle_is_filtered','torque_is_filtered','angle_filter_Hz', ...
            'torque_filter_Hz','filter_method', ...
            'mvc_t0','mvc_t1','idx0','idx1','velocity_thresh','nSegments', ...
            '-v7.3');

        res(k,:) = {participantName, mvcNameStr, peakTorque, meanTorque, ...
                    peakAtPct, windowDur, windowROM, nSegments, passiveRange, "ok"};

        fprintf(['[%2d/%2d] %-24s peak %6.2f Nm (at %4.1f%%) | mean %6.2f | ' ...
                 '%.2f s | ROM %5.2f deg | %d seg\n'], ...
            k, nF, tag, peakTorque, peakAtPct, meanTorque, windowDur, ...
            windowROM, nSegments);

    catch ME
        res(k,:) = {participantName, mvcNameStr, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, string(ME.message)};
        fprintf('[%2d/%2d] %-24s FAILED: %s\n', k, nF, tag, ME.message);
    end
end

%% ---------------------------
%  (4) Summary
%% ---------------------------
ok = res.status == "ok";
fprintf('\n===== STAGE 2 (MVC) SUMMARY =====\n');
fprintf('Processed %d | ok %d | failed %d\n', nF, nnz(ok), nnz(~ok));

if any(ok)
    fprintf('\nPeak MVC torque : %.2f +/- %.2f Nm (range %.2f to %.2f)\n', ...
        mean(res.peakTorque_Nm(ok)), std(res.peakTorque_Nm(ok)), ...
        min(res.peakTorque_Nm(ok)), max(res.peakTorque_Nm(ok)));
    fprintf('Window duration : %.2f +/- %.2f s\n', ...
        mean(res.windowDur_s(ok)), std(res.windowDur_s(ok)));
    fprintf('Window ROM      : %.2f +/- %.2f deg\n', ...
        mean(res.windowROM_deg(ok)), std(res.windowROM_deg(ok)));

    multi = res(ok & res.nSegments > 1, :);
    if isempty(multi)
        fprintf('Every recording held a single contraction.\n');
    else
        fprintf('\nRecordings with more than one velocity segment, check the QC figure:\n');
        disp(multi(:, {'participant','mvcFile','nSegments','windowDur_s'}));
    end

    fprintf('\nPeak torque per participant:\n');
    disp(res(ok, {'participant','peakTorque_Nm','peakAt_pctWindow','windowROM_deg'}));
end

if any(~ok)
    fprintf('\nFailures:\n');
    disp(res(~ok, {'participant','mvcFile','status'}));
end

writetable(res, fullfile(ROOT,'stage2_mvc_summary.csv'));
fprintf('\nSaved summary: %s\n', fullfile(ROOT,'stage2_mvc_summary.csv'));
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

function n = count_segments(above, time, minGap_s)
% Counts runs of samples above threshold, merging any two runs separated by
% less than minGap_s. A count above one means the recording holds more than one
% attempt, in which case the window spanning first to last crossing would cover
% the gap between them.
    d  = diff([false; above(:); false]);
    s_ = find(d ==  1);
    e_ = find(d == -1) - 1;

    if isempty(s_), n = 0; return; end

    n = 1;
    lastEnd = e_(1);
    for j = 2:numel(s_)
        if time(s_(j)) - time(lastEnd) >= minGap_s
            n = n + 1;
        end
        lastEnd = e_(j);
    end
end

function plot_mvc(time, angle, vel, torque_filt, torque_active, ...
                  idx0, idx1, vthr, peakTorque, tag, qcDir)
% Joint angle, angular velocity and torque before and after passive
% subtraction, with the detected contraction window marked on all three.

    f = figure('Color','w','Visible','off','Position',[100 100 1100 800]);
    tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');

    a1 = nexttile; hold(a1,'on');
    plot(a1, time, angle, 'Color',[0 .55 0], 'LineWidth',1.5);
    ylabel(a1,'Angle (deg)','FontWeight','bold');
    set(a1,'XTickLabel',[]);
    title(a1, tag + sprintf("  |  peak active torque %.2f Nm", peakTorque), ...
        'Interpreter','none','FontWeight','bold');

    a2 = nexttile; hold(a2,'on');
    plot(a2, time, vel, 'Color',[.25 .25 .25], 'LineWidth',1.2);
    yline(a2, -vthr, 'm--','LineWidth',1.0);
    ylabel(a2,'Vel (deg/s)','FontWeight','bold');
    set(a2,'XTickLabel',[]);

    a3 = nexttile; hold(a3,'on');
    hRaw = plot(a3, time, torque_filt,   'LineWidth',1.0, 'Color',[.6 .6 .6]);
    hNet = plot(a3, time, torque_active, 'LineWidth',1.5, 'Color',[.5 0 .5]);
    ylabel(a3,'Torque (Nm)','FontWeight','bold');
    xlabel(a3,'Time (s)','FontWeight','bold');

    for ax = [a1 a2 a3]
        xline(ax, time(idx0), 'k--','LineWidth',1.3);
        xline(ax, time(idx1), 'k--','LineWidth',1.3);
    end

    legend(a3, [hRaw hNet], {'Filtered torque','Active torque'}, ...
        'Location','best','Box','off');

    linkaxes([a1 a2 a3],'x');
    set([a1 a2 a3],'Box','off','TickDir','out','FontWeight','bold');

    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + ".png"), ...
        'Resolution',110);
    close(f);
end
