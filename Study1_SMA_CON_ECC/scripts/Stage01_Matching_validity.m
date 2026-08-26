%% ===== Stage 01: Torque-to-ramp matching and constant velocity window detection =====
%
% PURPOSE
%   For every retained set, locate each repetition, find the span of constant
%   angular velocity within it, and quantify how closely the produced torque
%   followed the prescribed ramp inside that span.
%
% SUPPORTS
%   The constant velocity phase definition and the offline torque matching
%   check reported in the Data Analysis section of the manuscript.
%
% INPUT
%   <ROOT>/<P>/<COND>/<Set_N>/<P>_<COND>_<Set_N>_sync_rawVolts.mat   (Stage 00)
%   <ROOT>/<P>/Ramps/{Con,Ecc}Match_*_{75,90}.txt
%
%   The ramp file is matched on contraction mode and intensity only. The "P1"
%   in the ramp filenames does not track the participant folder, so it is
%   ignored.
%
%   MVC recordings are skipped. They sit outside Set_* folders and have no ramp.
%
% OUTPUT
%   <set folder>/<name>_stage1_matchQC.mat            (sync + repQC)
%   <ROOT>/QC_Stage1/<tag>_matching.png
%   <ROOT>/QC_Stage1/<tag>_cvwindow.png
%   <ROOT>/QC_Stage1/SFig_CV_threshold.{png,pdf}      (supplementary figure)
%   <ROOT>/stage1_summary.csv
%
% METHOD
%   Repetition onsets come from sync.rep_on, the native-rate DAC1 edges. Within
%   each repetition the constant velocity window runs between the first and the
%   last sample where angular velocity exceeds the threshold in the direction of
%   the rotation, negative for CON dorsiflexion and positive for ECC
%   plantarflexion. Matching error is the median symmetrized percent
%   differences (SPDs) between measured torque and the prescribed ramp inside that window.
%
% ON THE VELOCITY THRESHOLD
%   The threshold sits slightly below the commanded 15 deg/s because velocity
%   oscillated around the commanded value throughout the plateau, so a threshold
%   at exactly 15 deg/s truncated the end of the detected phase. Section (7)
%   draws the supplementary figure that demonstrates this >> (S. Fig.S1)
%
% ON THE MATCHING CRITERION
%   Sets were accepted or repeated during testing on the basis of visual
%   inspection against the prescribed ramp. The offline percent difference
%   computed here is a verification statistic rather than an exclusion rule, and
%   no set is dropped on the basis of it.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set ROOT (and optionally CFG) beforehand and the prompts
%   are skipped. A file that errors is logged and the batch continues.
%
% DEPENDENCIES
%   Wfilt.m        zero-phase Butterworth wrapper              (see /functions)

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'velocity_thresh',     14);    % deg/s, CV window detection
CFG = set_default(CFG, 'validity_thresh_pct', 10);    % SPD verification criterion
CFG = set_default(CFG, 'fc_angle',            6);     % Hz
CFG = set_default(CFG, 'fc_torque',           20);    % Hz
CFG = set_default(CFG, 'eps_den',             1e-6);  % guard for the SPD denominator
CFG = set_default(CFG, 'saveFigures',         true);

% Supplementary figure, section (7)
CFG = set_default(CFG, 'makeThresholdFigure', true);
CFG = set_default(CFG, 'demoParticipant',     '08');
CFG = set_default(CFG, 'demoCondition',       'CON_90');
CFG = set_default(CFG, 'demoSet',             'Set_2');
CFG = set_default(CFG, 'demoRep',             2);
CFG = set_default(CFG, 'demoThresholds',      [15 14]);   % panel A, panel B

%% ---------------------------
%  (2) Data root and file list
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 1: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

assert(~isempty(which('Wfilt')), ...
    'Stage 1: Wfilt not found on the MATLAB path. Add the /functions folder.');

qcDir = fullfile(ROOT,'QC_Stage1');
if CFG.saveFigures && ~exist(qcDir,'dir'), mkdir(qcDir); end

F = dir(fullfile(ROOT, '**', '*_sync_rawVolts.mat'));

% Trials only. MVC recordings and stray saves sit outside Set_* folders.
inSet = arrayfun(@(x) ~isempty(regexp(x.folder,'[\\/]Set_\d+$','once')), F);
F = F(inSet);
assert(~isempty(F), 'Stage 1: no trial sync files found under %s', ROOT);

% Stage 0 writes <P>_<COND>_<Set>_sync_rawVolts.mat. Anything else in a Set
% folder is left over from an earlier run and would otherwise be processed
% twice, so it is reported and skipped.
expected = arrayfun(@(x) matches_stage0_name(x), F);
stale = F(~expected);
F = F(expected);

fprintf('\n===== Stage 1: torque matching and CV windows =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Found %d trial sync files.\n', numel(F));

if ~isempty(stale)
    fprintf('\nWARNING: %d file(s) do not match the Stage 0 naming and were skipped.\n', ...
        numel(stale));
    fprintf('Delete these before relying on the results:\n');
    for k = 1:numel(stale)
        fprintf('   %s\n', fullfile(stale(k).folder, stale(k).name));
    end
    fprintf('\n');
end

%% ---------------------------
%  (3) Batch loop
%% ---------------------------
nF  = numel(F);
res = table('Size',[nF 11], ...
    'VariableTypes',{'string','string','string','string','double','double', ...
                     'double','logical','double','double','string'}, ...
    'VariableNames',{'participant','condition','set','rampFile','nReps', ...
                     'nRepsValid','medianErr_pct','pass','meanCVdur_s', ...
                     'meanCVrom_deg','status'});

for k = 1:nF
    fp  = fullfile(F(k).folder, F(k).name);
    tag = "";
    partName = ""; condName = ""; setName = "";

    try
        L = load(fp,'sync');
        assert(isfield(L,'sync'), 'No "sync" struct in file.');
        sync = L.sync;

        [condPath, setName]  = fileparts(F(k).folder);
        [partPath, condName] = fileparts(condPath);
        [~,        partName] = fileparts(partPath);

        % Stage 0 records the names it used, which take precedence over the
        % folder names in case a folder was renamed after synchronization.
        if isfield(sync,'meta')
            if isfield(sync.meta,'participantName'), partName = char(sync.meta.participantName); end
            if isfield(sync.meta,'conditionName'),   condName = char(sync.meta.conditionName);   end
            if isfield(sync.meta,'setName'),         setName  = char(sync.meta.setName);         end
        end
        tag = string(sprintf('%s_%s_%s', partName, condName, setName));

        tok = regexp(upper(condName),'^(CON|ECC)_(\d+)$','tokens','once');
        assert(~isempty(tok), 'Cannot parse condition "%s".', condName);
        condition = tok{1};
        intensity = tok{2};

        %% --- ramp file ---
        rampDir = fullfile(partPath,'Ramps');
        assert(isfolder(rampDir), 'No Ramps folder at %s', rampDir);

        R   = dir(fullfile(rampDir,'*.txt'));
        hit = arrayfun(@(x) contains(lower(x.name), lower(condition)) && ...
                            contains(x.name, ['_' intensity]), R);
        assert(nnz(hit) == 1, 'Expected 1 ramp for %s_%s in %s, found %d.', ...
            condition, intensity, rampDir, nnz(hit));

        ramp = readmatrix(fullfile(rampDir, R(hit).name), 'Delimiter', ',');
        t_ramp = ramp(:,1); y_ramp = ramp(:,2);
        [t_ramp, iu] = unique(t_ramp,'stable');
        y_ramp = y_ramp(iu);
        t_ramp = t_ramp - t_ramp(1);
        ramp_duration = t_ramp(end);

        %% --- signals ---
        time = sync.time(:);
        Fs   = sync.Fs_emg;

        angle_smooth  = real(Wfilt(sync.angle(:),      CFG.fc_angle,  'low', Fs));
        torque_smooth = real(Wfilt(sync.torque_raw(:), CFG.fc_torque, 'low', Fs));

        assert(isfield(sync,'rep_on') && ~isempty(sync.rep_on), 'No rep_on markers.');
        onset_times = sync.rep_on(:);
        onset_idx   = nearest_idx(time, onset_times);
        nReps       = numel(onset_idx);

        angle_velocity = [0; diff(angle_smooth)./diff(time)];

        %% --- repetition windows and CV detection ---
        repIdx = nan(nReps,2);
        CVwin  = nan(nReps,2);

        for r = 1:nReps
            i0 = onset_idx(r);
            iRamp = find(time >= time(i0) + ramp_duration, 1, 'first');
            if isempty(iRamp), iRamp = numel(time); end
            if r < nReps, i1 = min(iRamp, onset_idx(r+1)-1); else, i1 = iRamp; end

            repIdx(r,:) = [i0 i1];
            t_rel = time(i0:i1) - time(i0);
            v     = angle_velocity(i0:i1);

            % Direction of the rotation: dorsiflexion is negative during CON,
            % plantarflexion positive during ECC.
            if strcmp(condition,'ECC'), above = v >  CFG.velocity_thresh;
            else,                       above = v < -CFG.velocity_thresh;
            end

            a = find(above,1,'first'); b = find(above,1,'last');
            if ~isempty(a) && ~isempty(b) && b > a
                CVwin(r,:) = [t_rel(a) t_rel(b)];
            end
        end

        %% --- matching error inside each CV window ---
        matching_errors = nan(nReps,1);
        n_valid_points  = nan(nReps,1);
        CVdur = nan(nReps,1); CVrom = nan(nReps,1);
        CVg   = nan(nReps,2);

        time_CV_raw = cell(nReps,1); angle_CV_raw = cell(nReps,1);
        torque_CV_raw = cell(nReps,1); emg_CV_raw = cell(nReps,1);

        for r = 1:nReps
            if ~all(isfinite(CVwin(r,:))), continue; end

            i0 = repIdx(r,1); i1 = repIdx(r,2);
            t_rep = time(i0:i1);
            t_rel = t_rep - t_rep(1);
            CVg(r,:) = t_rep(1) + CVwin(r,:);

            sel = find(t_rel >= CVwin(r,1) & t_rel <= CVwin(r,2));
            if numel(sel) < 5, continue; end
            idx = i0 + sel - 1;

            tqv = torque_smooth(idx);
            rav = interp1(t_ramp, y_ramp, t_rel(sel), 'linear', NaN);
            ok  = isfinite(tqv) & isfinite(rav);
            if nnz(ok) < 5, continue; end

            den = max(abs(tqv(ok)) + abs(rav(ok)), CFG.eps_den);
            matching_errors(r) = median(abs(tqv(ok)-rav(ok))./den*100);
            n_valid_points(r)  = nnz(ok);

            CVdur(r) = CVwin(r,2) - CVwin(r,1);
            CVrom(r) = abs(angle_smooth(idx(end)) - angle_smooth(idx(1)));

            time_CV_raw{r}   = time(idx);
            angle_CV_raw{r}  = angle_smooth(idx);
            torque_CV_raw{r} = torque_smooth(idx);
            if isfield(sync,'emg') && ~isempty(sync.emg)
                emg_CV_raw{r} = sync.emg(idx,:);
            end
        end

        median_error_all = median(matching_errors,'omitnan');
        nValid = nnz(isfinite(matching_errors));

        % Verification only, see the header. Rounded to one decimal so that the
        % flag agrees with the value printed in the summary.
        pass = round(median_error_all,1) <= CFG.validity_thresh_pct;

        %% --- QC figures ---
        if CFG.saveFigures
            plot_matching(time, angle_smooth, torque_smooth, t_ramp, y_ramp, ...
                repIdx, CVg, tag, median_error_all, CFG.validity_thresh_pct, ...
                pass, qcDir);
            plot_cv_windows(time, angle_velocity, repIdx, CVwin, ...
                CFG.velocity_thresh, tag, qcDir);
        end

        %% --- package and save ---
        repQC = struct();
        repQC.onset_times      = onset_times;
        repQC.onset_idx        = onset_idx;
        repQC.repIdx           = repIdx;
        repQC.CVwin_time       = CVwin;
        repQC.CVwin_globalTime = CVg;
        repQC.time_CV_raw      = time_CV_raw;
        repQC.angle_CV_raw     = angle_CV_raw;
        repQC.torque_CV_raw    = torque_CV_raw;
        repQC.emg_CV_raw       = emg_CV_raw;
        repQC.matchErr_pct_med    = matching_errors;
        repQC.median_matchErr_pct = median_error_all;
        repQC.mean_matchErr_pct   = mean(matching_errors,'omitnan');
        repQC.overall_match_pass  = pass;
        repQC.n_valid_points      = n_valid_points;
        repQC.CVdur_s             = CVdur;
        repQC.CVrom_deg           = CVrom;

        repQC.params = struct('condition',condition, ...
            'intensity',str2double(intensity), ...
            'velocity_thresh',CFG.velocity_thresh, ...
            'validity_thresh_pct',CFG.validity_thresh_pct, ...
            'eps_den',CFG.eps_den,'fc_angle',CFG.fc_angle,'fc_torque',CFG.fc_torque, ...
            'ramp_duration',ramp_duration,'rampFile',string(R(hit).name), ...
            'Fs_used',Fs, ...
            'createdOn',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'error_metric','median symmetrized percent difference within CV window');

        [~, base] = fileparts(F(k).name);
        save(fullfile(F(k).folder, sprintf('%s_stage1_matchQC.mat', base)), ...
             'sync','repQC','-v7.3');

        res(k,:) = {string(partName), string(condName), string(setName), ...
                    string(R(hit).name), nReps, nValid, median_error_all, ...
                    pass, mean(CVdur,'omitnan'), mean(CVrom,'omitnan'), "ok"};

        fprintf('[%3d/%3d] %-26s %d/%d reps | SPD %5.1f%% | pass %d | ROM %5.2f deg\n', ...
            k, nF, tag, nValid, nReps, median_error_all, pass, mean(CVrom,'omitnan'));

    catch ME
        res(k,:) = {string(partName), string(condName), string(setName), "", ...
                    NaN, NaN, NaN, false, NaN, NaN, string(ME.message)};
        fprintf('[%3d/%3d] %-26s FAILED: %s\n', k, nF, tag, ME.message);
    end
end

%% ---------------------------
%  (4) Summary
%% ---------------------------
ok = res.status == "ok";
fprintf('\n===== STAGE 1 SUMMARY =====\n');
fprintf('Processed %d | ok %d | failed %d\n', nF, nnz(ok), nnz(~ok));

if any(ok)
    fprintf('Median SPD across files : %.2f%% (range %.2f to %.2f)\n', ...
        median(res.medianErr_pct(ok)), min(res.medianErr_pct(ok)), max(res.medianErr_pct(ok)));

    bad = res(ok & ~res.pass, :);
    if isempty(bad)
        fprintf('All files met the %.0f%% verification criterion.\n', CFG.validity_thresh_pct);
    else
        fprintf('\nFiles above the verification criterion:\n');
        disp(bad(:, {'participant','condition','set','medianErr_pct'}));
    end

    drop = res(ok & res.nRepsValid < res.nReps, :);
    if ~isempty(drop)
        fprintf('\nFiles with repetitions lacking a valid CV window:\n');
        disp(drop(:, {'participant','condition','set','nRepsValid','nReps'}));
    end

    % Reported because a difference in angular excursion between modes would
    % confound any comparison of mechanical work.
    isC = ok & startsWith(res.condition,"CON");
    isE = ok & startsWith(res.condition,"ECC");
    fprintf('\nCV ROM  CON %.2f +/- %.2f deg | ECC %.2f +/- %.2f deg\n', ...
        mean(res.meanCVrom_deg(isC)), std(res.meanCVrom_deg(isC)), ...
        mean(res.meanCVrom_deg(isE)), std(res.meanCVrom_deg(isE)));
    fprintf('CV dur  CON %.3f +/- %.3f s   | ECC %.3f +/- %.3f s\n', ...
        mean(res.meanCVdur_s(isC)), std(res.meanCVdur_s(isC)), ...
        mean(res.meanCVdur_s(isE)), std(res.meanCVdur_s(isE)));
end

if any(~ok)
    fprintf('\nFailures:\n');
    disp(res(~ok, {'participant','condition','set','status'}));
end

writetable(res, fullfile(ROOT,'stage1_summary.csv'));
fprintf('\nSaved summary: %s\n', fullfile(ROOT,'stage1_summary.csv'));
if CFG.saveFigures, fprintf('QC figures   : %s\n', qcDir); end

%% ---------------------------
%  (5) Supplementary figure: velocity threshold demonstration
%% ---------------------------
%  One representative repetition, showing the constant velocity phase detected
%  at the commanded 15 deg/s (A) and at the threshold used here (B). Velocity
%  oscillated around the commanded value throughout the plateau, so the higher
%  threshold truncated the end of the detected phase.

if CFG.makeThresholdFigure
    demoFile = dir(fullfile(ROOT, CFG.demoParticipant, CFG.demoCondition, ...
        CFG.demoSet, '*_stage1_matchQC.mat'));

    if isempty(demoFile)
        warning('Stage 1: demo file not found for %s %s %s. Threshold figure skipped.', ...
            CFG.demoParticipant, CFG.demoCondition, CFG.demoSet);
    else
        Ld = load(fullfile(demoFile(1).folder, demoFile(1).name), 'sync', 'repQC');

        tD   = Ld.sync.time(:);
        angD = real(Wfilt(Ld.sync.angle(:), CFG.fc_angle, 'low', Ld.sync.Fs_emg));
        velD = [0; diff(angD)./diff(tD)];

        i0  = Ld.repQC.repIdx(CFG.demoRep,1);
        i1  = Ld.repQC.repIdx(CFG.demoRep,2);
        idx = i0:i1;

        tRel  = tD(idx) - tD(i0);
        v     = velD(idx);
        isECC = strcmpi(Ld.repQC.params.condition, 'ECC');

        figV = figure('Color','w','Name','Supplementary: CV threshold', ...
                      'Position',[100 100 800 700]);
        tlV = tiledlayout(figV,2,1,'TileSpacing','compact','Padding','compact');

        letters = {'A','B'};
        fprintf('\nCV threshold demonstration (%s %s %s, repetition %d):\n', ...
            CFG.demoParticipant, CFG.demoCondition, CFG.demoSet, CFG.demoRep);

        for kk = 1:2
            thr = CFG.demoThresholds(kk);
            axV = nexttile(tlV,kk); hold(axV,'on');

            if isECC, above = v > thr; else, above = v < -thr; end

            % Vertical lines mark the first and the last crossing, which is how
            % the window is defined in section (3).
            a_ = find(above, 1, 'first');
            b_ = find(above, 1, 'last');
            if ~isempty(a_) && ~isempty(b_) && b_ > a_
                xline(axV, tRel(a_), '--', 'Color',[0 0.4470 0.7410], ...
                    'LineWidth',2, 'HandleVisibility','off');
                xline(axV, tRel(b_), '--', 'Color',[0 0.4470 0.7410], ...
                    'LineWidth',2, 'HandleVisibility','off');
            end

            plot(axV, tRel, v, 'k', 'LineWidth',2.0);
            yline(axV, 0, ':k', 'LineWidth',1);
            if isECC
                yline(axV,  thr, '--r', 'LineWidth',2);
            else
                yline(axV, -thr, '--r', 'LineWidth',2);
            end

            ylabel(axV, 'Angular velocity [\circ/s]', ...
                'FontWeight','bold','FontSize',12);
            ylim(axV, [-25 40]);
            yticks(axV, -20:10:40);
            xlim(axV, [0 tRel(end)]);
            box(axV,'off');
            set(axV,'TickDir','out','LineWidth',2.5, ...
                'FontWeight','bold','FontSize',12);

            if kk == 1
                set(axV,'XTick',[]);
            else
                xlabel(axV,'Contraction window [s]', ...
                    'FontWeight','bold','FontSize',12);
            end

            text(axV, 0.01, 0.97, letters{kk}, 'Units','normalized', ...
                'FontSize',16, 'FontWeight','bold', ...
                'HorizontalAlignment','left','VerticalAlignment','top');

            text(axV, 0.98, 0.97, sprintf('%d\\circ/s', thr), ...
                'Units','normalized','FontSize',12,'FontWeight','bold', ...
                'Color','r','HorizontalAlignment','right', ...
                'VerticalAlignment','top');

            if ~isempty(a_) && ~isempty(b_)
                fprintf('  %2d deg/s: window %.3f to %.3f s (duration %.3f s)\n', ...
                    thr, tRel(a_), tRel(b_), tRel(b_) - tRel(a_));
            end
        end

        if CFG.saveFigures
            exportgraphics(figV, fullfile(qcDir,'SFig_CV_threshold.png'), ...
                'Resolution', 300);
            exportgraphics(figV, fullfile(qcDir,'SFig_CV_threshold.pdf'), ...
                'ContentType','vector');
            fprintf('Threshold figure saved to %s\n', qcDir);
        end
    end
end

%% =========================================================
%  Local functions
%% =========================================================
function S = set_default(S, fieldName, value)
% Fills a config field only when the caller has not already supplied it.
    if ~isfield(S, fieldName) || isempty(S.(fieldName))
        S.(fieldName) = value;
    end
end

function tf = matches_stage0_name(f)
% Stage 0 writes <P>_<COND>_<Set>_sync_rawVolts.mat into the set folder.
    [~, base] = fileparts(f.name);
    [condPath, setName]  = fileparts(f.folder);
    [~,        condName] = fileparts(condPath);
    tf = endsWith(base, sprintf('%s_%s_sync_rawVolts', condName, setName));
end

function idx = nearest_idx(t, q)
    idx = zeros(numel(q),1);
    for i = 1:numel(q)
        [~, idx(i)] = min(abs(t - q(i)));
    end
end

function plot_matching(time, angS, tqS, t_ramp, y_ramp, repIdx, CVg, ...
                       tag, medErr, crit, pass, qcDir)
% Measured torque against the prescribed ramp, with the CV window of each
% repetition marked and the difference between the two shaded.
    f = figure('Color','w','Visible','off','Position',[100 100 1100 600]);
    tl = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');

    a1 = nexttile(tl,1);
    plot(a1, time, angS, 'g-', 'LineWidth',1.3);
    ylabel(a1,'Angle (deg)','FontWeight','bold');
    set(a1,'Box','off','XColor','none','TickDir','out','FontWeight','bold');
    title(a1, sprintf('%s | median SPD %.1f%% (criterion %.0f%%) | pass %d', ...
        tag, medErr, crit, pass), 'FontWeight','bold','Interpreter','none');

    a2 = nexttile(tl,2); hold(a2,'on');
    plot(a2, time, tqS, 'Color',[.6 .6 .6], 'LineWidth',0.5);

    for r = 1:size(repIdx,1)
        if ~all(isfinite(repIdx(r,:))), continue; end
        idx   = repIdx(r,1):repIdx(r,2);
        t_rep = time(idx);
        t_rel = t_rep - t_rep(1);
        rampR = interp1(t_ramp, y_ramp, t_rel, 'linear', NaN);

        plot(a2, t_rep, rampR, 'b-', 'LineWidth',1.2);
        plot(a2, t_rep, tqS(idx), '-', 'Color',[.5 0 .5], 'LineWidth',1.2);

        if all(isfinite(CVg(r,:)))
            m = t_rep >= CVg(r,1) & t_rep <= CVg(r,2);
            if nnz(m) > 2
                tv = t_rep(m); rv = rampR(m); qv = tqS(idx); qv = qv(m);
                g = isfinite(rv) & isfinite(qv);
                if nnz(g) > 2
                    fill(a2, [tv(g); flipud(tv(g))], [rv(g); flipud(qv(g))], ...
                        'r', 'FaceAlpha',0.15, 'EdgeColor','none');
                end
            end
            xline(a2, CVg(r,1), 'k--'); xline(a2, CVg(r,2), 'k--');
        end
    end

    xlabel(a2,'Time (s)','FontWeight','bold');
    ylabel(a2,'Torque (Nm)','FontWeight','bold');
    set(a2,'Box','off','TickDir','out','FontWeight','bold');

    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + "_matching.png"), ...
        'Resolution',110);
    close(f);
end

function plot_cv_windows(time, vel, repIdx, CVwin, vthr, tag, qcDir)
% One row per repetition, so a repetition with no valid window is obvious.
    n = size(repIdx,1);
    f = figure('Color','w','Visible','off','Position',[100 100 900 max(400,120*n)]);
    tiledlayout(f, n, 1, 'TileSpacing','compact','Padding','compact');

    for r = 1:n
        ax = nexttile; hold(ax,'on');
        if ~all(isfinite(repIdx(r,:))), continue; end
        idx   = repIdx(r,1):repIdx(r,2);
        t_rel = time(idx) - time(idx(1));

        plot(ax, t_rel, vel(idx), 'k', 'LineWidth',1.2);
        yline(ax, vthr,'--r'); yline(ax,-vthr,'--r');
        if all(isfinite(CVwin(r,:)))
            xline(ax, CVwin(r,1), '--b','LineWidth',1.5);
            xline(ax, CVwin(r,2), '--b','LineWidth',1.5);
        end
        ylabel(ax,'\omega'); title(ax, sprintf('Rep %02d', r));
        set(ax,'Box','off','TickDir','out');
    end
    xlabel('Rep time (s)');
    sgtitle(f, tag + " | CV windows", 'FontWeight','bold','Interpreter','none');

    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + "_cvwindow.png"), ...
        'Resolution',110);
    close(f);
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