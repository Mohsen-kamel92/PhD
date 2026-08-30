%% ===== Stage 00: Synchronize SPIKE torque and angle with HDsEMG =====
%
% PURPOSE
%   Bring the dynamometer recording (SPIKE, 2 kHz) and the high-density EMG
%   recording (MEACS, ~2048 Hz) onto a common time base, so that every later
%   stage can index torque, joint angle and EMG at the same sample.
%
% SUPPORTS
%   The synchronization described in the Data Collection section of the
%   manuscript, including the three recordings in which the auxiliary channel
%   missed the opening trigger pulses.
%
% INPUT
%   Raw SPIKE .mat files and MEACS .sig folders, laid out as below.
%
% OUTPUT
%   <set folder>/<P>_<COND>_<Set>_sync_rawVolts.mat   (struct "sync")
%   <ROOT>/QC_Stage0/<tag>.png
%   <ROOT>/stage0_summary.csv
%
% EXPECTED LAYOUT
%   <ROOT>/<P>/<CON_75|CON_90|ECC_75|ECC_90>/<Set_N>/<spike>.mat
%   <ROOT>/<P>/zEMG/<CONDITION>/<Set_N>/                  (*.sig files)
%   <ROOT>/<P>/MVC/<mvcname>.mat                          (optional)
%   <ROOT>/<P>/zEMG/MVC/<mvcname>/                        (optional)
%
%   The EMG folder mirrors the SPIKE folder. For MVC recordings the SPIKE file
%   sits directly in the MVC folder and the EMG folder carries the same
%   basename as that file.
%
%   Every Set_* folder present is treated as retained data. Failed attempts
%   were deleted during collection rather than kept, so no exclusion logic is
%   applied here.
%
% TWO TIMING CORRECTIONS
%   Both systems are cropped to their own first trigger and re-zeroed, so
%   alignment rests on the assumption that the first AUX edge and the first
%   DAC edge are the same physical event. Two ways that can fail:
%
%   1) CONSTANT OFFSET. AUX misses the opening pulse or pulses, so EMG sample
%      zero corresponds to a later DAC time. The AUX-vs-DAC residual is then a
%      large, tightly clustered non-zero value with no trend. Detected from the
%      median residual (excluding the first onset, which is forced to time zero
%      when AUX starts high) and removed by shifting the EMG.
%
%   2) LINEAR DRIFT. The two systems run on independent clocks, so a small rate
%      difference makes the time bases diverge over a recording. The residual
%      then shows a linear trend, and the EMG is resampled using both slope and
%      intercept, t_DAC = a*t_AUX + b.
%
%   The constant offset is tested first, because a large step would otherwise
%   distort the trend fit.
%
% ON SAMPLING RATE
%   sync.Fs_emg holds the rate measured from the recording itself, which is the
%   value every later stage uses and the value reported in the manuscript.
%   sync.Fs_emg_nominal holds the manufacturer figure of 2048 Hz and is stored
%   for reference only.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set ROOT (and optionally CFG) beforehand and the prompts
%   are skipped. A recording that errors is logged and the batch continues.
%
% DEPENDENCIES
%   MAECS_read.m   reads a MEACS .sig folder into a struct     (see /functions)
%   Wfilt.m        zero-phase Butterworth wrapper              (see /functions)

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'th_dac',        0.85);   % V, DAC1 high threshold

CFG = set_default(CFG, 'offset_mode',   "auto"); % "auto" | "off"
CFG = set_default(CFG, 'offset_min_s',  0.020);  % median residual counting as real
CFG = set_default(CFG, 'offset_max_sd', 0.050);  % refuse above this scatter

CFG = set_default(CFG, 'drift_mode',    "auto"); % "auto" | "force" | "off"
CFG = set_default(CFG, 'drift_r2_min',  0.80);
CFG = set_default(CFG, 'drift_min_ms',  3.0);
CFG = set_default(CFG, 'drift_max_ppm', 2000);

CFG = set_default(CFG, 'includeMVC',    true);   % also sync the MVC recordings
CFG = set_default(CFG, 'saveFigures',   true);
CFG = set_default(CFG, 'fc_angle_qc',   6);      % Hz, QC display only
CFG = set_default(CFG, 'fc_torque_qc',  20);     % Hz, QC display only

%% ---------------------------
%  (2) Data root
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 0: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);
qcDir = fullfile(ROOT,'QC_Stage0');
if CFG.saveFigures && ~exist(qcDir,'dir'), mkdir(qcDir); end

for dep = {'MAECS_read','Wfilt'}
    assert(~isempty(which(dep{1})), ...
        'Stage 0: %s not found on the MATLAB path. Add the /functions folder.', dep{1});
end

fprintf('\n===== Stage 0: synchronization =====\n');
fprintf('Root: %s\n', ROOT);

%% ---------------------------
%  (3) Build the job list
%% ---------------------------
jobs = struct('spikeFile',{},'emgDir',{},'participant',{},'condition',{},'set',{});

P = dir(ROOT);
P = P([P.isdir] & ~startsWith({P.name},'.'));

for ip = 1:numel(P)
    pName = P(ip).name;
    pPath = fullfile(ROOT, pName);
    if ~isfolder(fullfile(pPath,'zEMG')), continue; end

    % ---- trial conditions
    C = dir(pPath);
    C = C([C.isdir]);
    for ic = 1:numel(C)
        cName = C(ic).name;
        if isempty(regexp(cName,'^(CON|ECC)_\d+$','once')), continue; end

        Sdir = dir(fullfile(pPath, cName, 'Set_*'));
        Sdir = Sdir([Sdir.isdir]);

        for is = 1:numel(Sdir)
            sName    = Sdir(is).name;
            setPath  = fullfile(pPath, cName, sName);
            emgDir   = fullfile(pPath, 'zEMG', cName, sName);
            spikeMat = find_spike_mat(setPath);

            if isempty(spikeMat) || ~isfolder(emgDir), continue; end

            jobs(end+1) = struct('spikeFile',spikeMat,'emgDir',emgDir, ...
                'participant',pName,'condition',cName,'set',sName); %#ok<SAGROW>
        end
    end

    % ---- MVC
    if CFG.includeMVC && isfolder(fullfile(pPath,'MVC'))
        M = dir(fullfile(pPath,'MVC','*.mat'));
        for im = 1:numel(M)
            spikeMat = fullfile(M(im).folder, M(im).name);
            if ~is_spike_file(spikeMat), continue; end
            [~, base] = fileparts(M(im).name);
            emgDir = fullfile(pPath, 'zEMG', 'MVC', base);
            if ~isfolder(emgDir), continue; end

            jobs(end+1) = struct('spikeFile',spikeMat,'emgDir',emgDir, ...
                'participant',pName,'condition','MVC','set',base); %#ok<SAGROW>
        end
    end
end

nJ = numel(jobs);
assert(nJ > 0, 'Stage 0: no SPIKE and EMG pairs found under %s', ROOT);
fprintf('Found %d recordings to synchronize.\n\n', nJ);

%% ---------------------------
%  (4) Batch loop
%% ---------------------------
res = table('Size',[nJ 12], ...
    'VariableTypes',{'string','string','string','double','double','double', ...
                     'double','logical','double','double','logical','string'}, ...
    'VariableNames',{'participant','condition','set','overlap_s','nOnsets', ...
                     'resid_med_ms','resid_sd_ms','offset_applied','offset_ms', ...
                     'drift_r2','drift_applied','status'});

for k = 1:nJ
    J   = jobs(k);
    tag = string(sprintf('%s_%s_%s', J.participant, J.condition, J.set));

    try
        [sync, info] = sync_one(J, CFG);

        outName = fullfile(fileparts(J.spikeFile), ...
            sprintf('%s_sync_rawVolts.mat', regexprep(tag,'[^\w-]','_')));
        save(outName, 'sync', '-v7.3');

        if CFG.saveFigures
            make_qc_figure(sync, info, tag, qcDir, CFG.fc_angle_qc, CFG.fc_torque_qc);
        end

        res(k,:) = {string(J.participant), string(J.condition), string(J.set), ...
                    info.overlap_s, info.nOnsets, info.resid_med_ms, info.resid_sd_ms, ...
                    info.offset.applied, info.offset.shift_s*1e3, ...
                    info.drift.r2, info.drift.applied, "ok"};

        note = "";
        if info.offset.applied
            note = note + sprintf(" | OFFSET %+.0f ms", info.offset.shift_s*1e3);
        end
        if info.drift.applied
            note = note + " | DRIFT CORRECTED";
        end

         if info.nOnsets < 3
            fprintf('[%3d/%3d] %-26s %d onset | timing check not applicable%s\n', ...
                k, nJ, tag, info.nOnsets, note);
        else
            fprintf('[%3d/%3d] %-26s %d onsets | resid %+7.2f +/- %5.2f ms%s\n', ...
                k, nJ, tag, info.nOnsets, info.resid_med_ms, info.resid_sd_ms, note);
        end
    catch ME
        res(k,:) = {string(J.participant), string(J.condition), string(J.set), ...
                    NaN, NaN, NaN, NaN, false, NaN, NaN, false, string(ME.message)};
        fprintf('[%3d/%3d] %-26s FAILED: %s\n', k, nJ, tag, ME.message);
    end
end

%% ---------------------------
%  (5) Summary
%% ---------------------------
ok = res.status == "ok";
fprintf('\n===== STAGE 0 SUMMARY =====\n');
fprintf('Processed %d | ok %d | failed %d\n', nJ, nnz(ok), nnz(~ok));

if any(ok)
    fprintf('AUX-DAC residual : median %+.3f ms | SD median %.3f ms\n', ...
        median(res.resid_med_ms(ok),'omitnan'), median(res.resid_sd_ms(ok),'omitnan'));
    fprintf('Offset applied   : %d file(s)\n', nnz(res.offset_applied(ok)));
    fprintf('Drift applied    : %d file(s)\n', nnz(res.drift_applied(ok)));

    corrected = res(ok & (res.offset_applied | res.drift_applied), :);
    if ~isempty(corrected)
        fprintf('\nCorrected files:\n');
        disp(corrected(:, {'participant','condition','set','offset_ms','drift_applied'}));
    end

    noisy = res(ok & res.resid_sd_ms > 2, :);
    if ~isempty(noisy)
        fprintf('\nResidual SD above 2 ms after correction (AUX trigger quality):\n');
        disp(noisy(:, {'participant','condition','set','resid_med_ms','resid_sd_ms'}));
    end
end

if any(~ok)
    fprintf('\nFailures:\n');
    disp(res(~ok, {'participant','condition','set','status'}));
end

writetable(res, fullfile(ROOT,'stage0_summary.csv'));
fprintf('\nSaved summary: %s\n', fullfile(ROOT,'stage0_summary.csv'));
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
function tf = is_spike_file(fp)
% A SPIKE export is identified by the three channels the pipeline needs.
    tf = false;
    try
        v = who('-file', fp);
        tf = all(ismember({'DAC1','Torque','Angle'}, v));
    catch
    end
end

function fp = find_spike_mat(folder)
% Returns the first .mat in the folder that looks like a SPIKE export, so that
% outputs written by this or any later stage are never picked up as input.
    fp = '';
    M = dir(fullfile(folder,'*.mat'));
    M = M(~contains({M.name},{'sync_rawVolts','stage1','stage2','stage3', ...
                              'stage4','stage5','stage6','ActiveTorque'}));
    for i = 1:numel(M)
        cand = fullfile(M(i).folder, M(i).name);
        if is_spike_file(cand), fp = cand; return; end
    end
end

function [sync, info] = sync_one(J, CFG)

    %% --- SPIKE ---
    S = load(J.spikeFile);
    for k = {'DAC1','Torque','Angle'}
        assert(isfield(S,k{1}), 'Missing "%s" in SPIKE file.', k{1});
    end

    dv = S.DAC1.values(:);
    i1 = find(dv > CFG.th_dac, 1, 'first');
    i2 = find(dv > CFG.th_dac, 1, 'last');
    assert(~isempty(i1) && ~isempty(i2) && i2 > i1, 'No valid DAC1-high span.');

    fn = fieldnames(S);
    for i = 1:numel(fn)
        f = fn{i};
        if isstruct(S.(f)) && isfield(S.(f),'times') && isfield(S.(f),'values')
            S.(f).times  = S.(f).times(i1:i2);
            S.(f).values = S.(f).values(i1:i2);
            S.(f).times  = S.(f).times - S.(f).times(1);
        end
    end

    %% --- HDsEMG ---
    HDEMG  = MAECS_read(J.emgDir);
    gNames = fieldnames(HDEMG);
    assert(isscalar(gNames), 'Expected one grid, found %d in %s', numel(gNames), J.emgDir);
    G = HDEMG.(gNames{1});

    timeEMG = G.Times(:);
    aux     = G.AUX(:);
    emg     = G.EMG;

    assert(numel(aux) == numel(timeEMG) && size(emg,1) == numel(timeEMG), ...
        'EMG, AUX and time length mismatch.');

    Fs_emg_actual = 1 / median(diff(timeEMG));

    aux(~isfinite(aux)) = 0;
    j1 = find(aux > 0, 1, 'first');
    j2 = find(aux > 0, 1, 'last');
    assert(~isempty(j1) && ~isempty(j2) && j2 > j1, 'No valid AUX-high span.');

    emg_vals = emg(j1:j2,:);
    emg_time = timeEMG(j1:j2) - timeEMG(j1);
    AUX_cut  = aux(j1:j2);

    %% --- align and interpolate ---
    % Both records were cropped to their own first trigger and re-zeroed, so
    % dt_align is zero by construction. It is retained only for provenance.
    tAux0    = emg_time(find(AUX_cut > 0, 1, 'first'));
    tDAC0    = S.DAC1.times(find(S.DAC1.values > CFG.th_dac, 1, 'first'));
    dt_align = tDAC0 - tAux0;

    tTorque = S.Torque.times - dt_align;
    tAngle  = S.Angle.times  - dt_align;
    tDAC    = S.DAC1.times   - dt_align;

    torque_i = interp1(tTorque, S.Torque.values, emg_time, 'linear');
    angle_i  = interp1(tAngle,  S.Angle.values,  emg_time, 'linear');
    dac1_i   = interp1(tDAC,    S.DAC1.values,   emg_time, 'linear');

    valid = isfinite(torque_i) & isfinite(angle_i) & isfinite(dac1_i);
    a = find(valid,1,'first'); b = find(valid,1,'last');
    assert(~isempty(a) && b > a, 'No valid overlap after interpolation.');

    keep     = a:b;
    time     = emg_time(keep);
    emg_vals = emg_vals(keep,:);
    AUX_cut  = AUX_cut(keep);
    torque_i = torque_i(keep);
    angle_i  = angle_i(keep);
    dac1_i   = dac1_i(keep);

    %% --- onsets at the native DAC rate ---
    dh = S.DAC1.values(:) > CFG.th_dac;
    ei = find(diff(dh) == 1) + 1;
    don = tDAC(ei);
    if dh(1), don = [tDAC(1); don(:)]; end

    dac_on = don(don >= time(1) & don <= time(end));
    dac_on = dac_on(:);
    assert(~isempty(dac_on), 'No onsets inside the retained window.');

    aux_on = aux_onsets(AUX_cut, time);

    %% --- 1. constant offset -----------------------------------------
    % The first AUX onset is forced to time(1) when AUX starts high, so it is
    % excluded from the estimate. A large median with tight scatter means the
    % whole EMG record sits at the wrong time relative to the DAC.
    offs = struct('applied',false,'shift_s',0,'median_s',NaN,'sd_s',NaN, ...
                  'mode',string(CFG.offset_mode),'reason',"");

    nO = min(numel(aux_on), numel(dac_on));
    resid_med_ms = NaN; resid_sd_ms = NaN;

    if nO >= 3
        r  = aux_on(2:nO) - dac_on(2:nO);
        md = median(r); sd = std(r);
        offs.median_s = md; offs.sd_s = sd;
        resid_med_ms = md*1e3; resid_sd_ms = sd*1e3;

        switch lower(string(CFG.offset_mode))
            case "off"
                offs.reason = "offset correction disabled";
            case "auto"
                if abs(md) > CFG.offset_min_s && sd < CFG.offset_max_sd
                    [time, emg_vals, AUX_cut, torque_i, angle_i, dac1_i] = ...
                        shift_emg(md, time, emg_vals, AUX_cut, torque_i, angle_i, dac1_i);

                    dac_on = dac_on(dac_on >= time(1) & dac_on <= time(end));
                    aux_on = aux_onsets(AUX_cut, time);

                    offs.applied = true;
                    offs.shift_s = md;
                    offs.reason  = "AUX start offset relative to DAC";

                    n2 = min(numel(aux_on), numel(dac_on));
                    if n2 >= 3
                        r2v = aux_on(2:n2) - dac_on(2:n2);
                        resid_med_ms = median(r2v)*1e3;
                        resid_sd_ms  = std(r2v)*1e3;
                    end
                elseif abs(md) > CFG.offset_min_s
                    offs.reason = sprintf('offset %.3f s but scatter too high (SD %.1f ms)', ...
                        md, sd*1e3);
                else
                    offs.reason = "no constant offset above threshold";
                end
        end
        elseif nO == 1
        offs.reason = "single onset, timing check not applicable";
    else
        offs.reason = sprintf('only %d comparable onsets', nO);
    end

    %% --- 2. linear drift --------------------------------------------
    drift = struct('applied',false,'slope',1,'intercept',0,'r2',NaN,'ppm',NaN, ...
                   'total_ms',NaN,'resid_pre',[],'resid_post',[], ...
                   'mode',string(CFG.drift_mode),'reason',"");

    nCmp = min(numel(aux_on), numel(dac_on));

    if nCmp >= 4
        A = aux_on(1:nCmp); D = dac_on(1:nCmp);
        r_pre = A - D;

        p  = polyfit(A, D, 1);
        q  = polyfit(A, r_pre, 1);
        rf = polyval(q, A);
        r2 = 1 - sum((r_pre-rf).^2) / max(sum((r_pre-mean(r_pre)).^2), eps);

        drift.slope = p(1); drift.intercept = p(2); drift.r2 = r2;
        drift.ppm = abs(p(1)-1)*1e6;
        drift.total_ms = abs(q(1))*(A(end)-A(1))*1e3;
        drift.resid_pre = r_pre;

        doIt = false;
        switch lower(string(CFG.drift_mode))
            case "off",   drift.reason = "disabled";
            case "force", doIt = true; drift.reason = "forced";
            case "auto"
                if r2 >= CFG.drift_r2_min && drift.total_ms >= CFG.drift_min_ms && ...
                        drift.ppm <= CFG.drift_max_ppm
                    doIt = true; drift.reason = "systematic linear trend detected";
                elseif drift.ppm > CFG.drift_max_ppm
                    drift.reason = sprintf('rate error implausible (%.0f ppm)', drift.ppm);
                else
                    drift.reason = "no systematic trend above threshold";
                end
        end

        if doIt
            tc = polyval(p, time);
            ec = interp1(tc, emg_vals, time, 'linear', NaN);
            ac = interp1(tc, AUX_cut,  time, 'nearest', NaN);
            g  = all(isfinite(ec),2) & isfinite(ac);

            time = time(g); emg_vals = ec(g,:); AUX_cut = ac(g);
            torque_i = torque_i(g); angle_i = angle_i(g); dac1_i = dac1_i(g);

            dac_on = dac_on(dac_on >= time(1) & dac_on <= time(end));
            aux_on = aux_onsets(AUX_cut, time);

            n2 = min(numel(aux_on), numel(dac_on));
            drift.resid_post = aux_on(1:n2) - dac_on(1:n2);
            drift.applied = true;

            if n2 >= 3
                r2v = aux_on(2:n2) - dac_on(2:n2);
                resid_med_ms = median(r2v)*1e3;
                resid_sd_ms  = std(r2v)*1e3;
            end
        end
        elseif nCmp == 1
        drift.reason = "single onset, timing check not applicable";
    else
        drift.reason = sprintf('only %d comparable onsets', nCmp);
    end

    %% --- pack ---
    sync = struct();
    sync.time   = time;
    sync.Fs_emg = Fs_emg_actual;      % measured; used by every later stage
    sync.Fs_emg_nominal = 2048;       % manufacturer figure, reference only
    sync.Fs_spike = 2000;             % metadata only

    sync.emg = emg_vals;  sync.AUX = AUX_cut;
    sync.torque_raw = torque_i;  sync.angle = angle_i;  sync.dac1 = dac1_i;
    sync.aux_on = aux_on;  sync.dac_on = dac_on;  sync.rep_on = dac_on;

    sync.meta = struct();
    sync.meta.method = ['crop to first trigger, interpolate to the EMG clock, ' ...
                        'crop to overlap; onsets taken from the native-rate DAC1; ' ...
                        'constant AUX offset and linear clock drift corrected ' ...
                        'where detected'];
    sync.meta.thr_dac         = CFG.th_dac;
    sync.meta.dt_align_s      = dt_align;
    sync.meta.participantName = J.participant;
    sync.meta.conditionName   = J.condition;
    sync.meta.setName         = J.set;
    sync.meta.spike_file      = J.spikeFile;
    sync.meta.emg_dir         = J.emgDir;
    sync.meta.createdOn       = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    sync.meta.onset_detection = "native DAC1 (2 kHz) rising edge, mapped to the EMG clock";
    sync.meta.const_offset    = offs;
    sync.meta.clock_drift     = drift;
    sync.meta.crop = struct('spike_i1',i1,'spike_i2',i2,'emg_j1',j1,'emg_j2',j2);
    sync.meta.config = CFG;

    info = struct('overlap_s',time(end)-time(1),'nOnsets',numel(dac_on), ...
                  'resid_med_ms',resid_med_ms,'resid_sd_ms',resid_sd_ms, ...
                  'offset',offs,'drift',drift);
end

function on = aux_onsets(AUX_cut, time)
    high = AUX_cut > 0;
    on = time(find(diff(high) == 1) + 1);
    if high(1), on = [time(1); on(:)]; end
    on = on(:);
end

function [t, e, x, tq, an, dc] = shift_emg(shift_s, t, e, x, tq, an, dc)
% Corrected EMG at time t is the original EMG at time t + shift_s.
    e2 = interp1(t, e, t + shift_s, 'linear', NaN);
    x2 = interp1(t, x, t + shift_s, 'nearest', NaN);
    g  = all(isfinite(e2),2) & isfinite(x2);
    t = t(g); e = e2(g,:); x = x2(g);
    tq = tq(g); an = an(g); dc = dc(g);
end

function make_qc_figure(sync, info, tag, qcDir, fcA, fcT)
% One page per recording, for visual confirmation that the two systems line up.
    Fs = sync.Fs_emg;
    angS = Wfilt(sync.angle,      fcA, 'low', Fs);
    tqS  = Wfilt(sync.torque_raw, fcT, 'low', Fs);
    w    = max(1, ceil(0.050*Fs));
    env  = mean(sqrt(movmean(sync.emg.^2, w, 1, 'Endpoints','shrink')), 2, 'omitnan');

    f = figure('Color','w','Visible','off','Position',[100 100 1200 800]);
    tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');

    a1 = nexttile;
    plot(a1, sync.time, angS, 'Color',[0 .55 0], 'LineWidth',1.5); hold(a1,'on');
    ylabel(a1,'Angle (deg)','FontWeight','bold'); set(a1,'XTickLabel',[]);
    for i = 1:numel(sync.dac_on), xline(a1, sync.dac_on(i), '--k'); end

    ttl = tag + sprintf("  |  resid %+.2f +/- %.2f ms", info.resid_med_ms, info.resid_sd_ms);
    if info.offset.applied
        ttl = ttl + sprintf("  |  offset %+.0f ms applied", info.offset.shift_s*1e3);
    end
    if info.drift.applied
        ttl = ttl + "  |  drift corrected";
    end
    title(a1, ttl, 'Interpreter','none','FontWeight','bold');

    a2 = nexttile;
    yyaxis(a2,'left');  plot(a2, sync.time, env, 'r', 'LineWidth',1.5);
    ylabel(a2,'EMG RMS (V)','FontWeight','bold'); a2.YAxis(1).Color = [1 0 0];
    yyaxis(a2,'right'); plot(a2, sync.time, tqS, 'Color',[.5 0 .5], 'LineWidth',1.5);
    ylabel(a2,'Torque (Nm)','FontWeight','bold'); a2.YAxis(2).Color = [.5 0 .5];
    set(a2,'XTickLabel',[]);

    a3 = nexttile; hold(a3,'on');
    n = min(numel(sync.aux_on), numel(sync.dac_on));
    if n >= 2
        plot(a3, sync.dac_on(2:n), (sync.aux_on(2:n)-sync.dac_on(2:n))*1e3, ...
             'o-','LineWidth',1.4);
    end
    yline(a3, 0, '--');
    ylabel(a3,'AUX - DAC (ms)','FontWeight','bold');
    xlabel(a3,'Time (s)','FontWeight','bold');

    linkaxes([a1 a2 a3],'x');
    set([a1 a2 a3],'Box','off','TickDir','out','FontWeight','bold');
    exportgraphics(f, fullfile(qcDir, regexprep(tag,'[^\w-]','_') + ".png"), 'Resolution',110);
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
