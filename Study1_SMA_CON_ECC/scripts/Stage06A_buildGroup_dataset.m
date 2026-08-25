%% ===== Stage 06A: Build the group dataset =====
%
% PURPOSE
%   Read every condition-level file, index it by participant and condition, and
%   reshape the outcomes into the arrays the statistics stages expect.
%
% SUPPORTS
%   Every result in the manuscript. This stage produces the single file that
%   Stages 6B, 6C1, 6C2 and 6C4 all read.
%
% INPUT
%   <ROOT>/Stage5B/<P>_<COND>_stage5B_condAgg.mat        (Stage 05B)
%
% OUTPUT
%   <ROOT>/groupData_stage6A.mat
%   <ROOT>/stage6A_check.csv
%
% DESIGN
%   2 x 2 repeated measures: contraction type (CON, ECC) by intensity (75, 90).
%   Conditions are indexed in a fixed order that every later stage relies on:
%       1 = CON_75   2 = ECC_75   3 = CON_90   4 = ECC_90
%
% WHAT THIS STAGE IS FOR
%   Stage 05B leaves one struct per participant and condition. The ANOVA and SPM
%   stages need every outcome as [nSub x nCond x nTime], so the reshaping is the
%   purpose. Participant and condition are read from each file's metadata rather
%   than from the order the files happen to be opened in, so the assignment
%   cannot silently go wrong.
%
% ON THE TWO SCALAR SETS
%   groupData.scalarCVP holds means over the analysis window and is what the
%   statistics use. groupData.scalar holds full-cycle means and is kept for
%   reference only.
%
% ON SAVING
%   MATLAB string arrays do not survive -v7.3 reliably, so participant names and
%   source paths are stored as cell arrays of char.
%
% RUNNING THIS STAGE
%   Standalone, run the script and choose the data root when prompted. From the
%   pipeline driver, set CFG.ROOT beforehand and the prompt is skipped. Unlike
%   the earlier stages this one fails rather than continuing, because a missing
%   or duplicated file would make the group arrays wrong rather than incomplete.
%
% DEPENDENCIES
%   None beyond base MATLAB.

close all;

%% ---------------------------
%  (1) Options
%% ---------------------------
% The driver may set any of these in CFG before calling. Fields are only filled
% in here when absent, so a one-button run and a standalone run behave the same.

if ~exist('CFG','var') || ~isstruct(CFG), CFG = struct(); end

CFG = set_default(CFG, 'requiredConditions', {'CON_75','ECC_75','CON_90','ECC_90'});
CFG = set_default(CFG, 'binsToUse',          {'all','early','late'});
CFG = set_default(CFG, 'cvp_lo',             20);   % analysis window, % of the CV phase
CFG = set_default(CFG, 'cvp_hi',             80);

requiredConditions = CFG.requiredConditions;
binsToUse          = CFG.binsToUse;
nCond              = numel(requiredConditions);

%% ---------------------------
%  (2) Data root and input files
%% ---------------------------
% The driver sets CFG.ROOT. Standalone, the folder is always chosen through the
% dialog, so a variable left over from an earlier run cannot silently redirect
% the analysis to the wrong data.
if ~isfield(CFG,'ROOT') || ~(ischar(CFG.ROOT) || isstring(CFG.ROOT)) || strlength(string(CFG.ROOT)) == 0
    picked = uigetdir(pwd, 'Select the ROOT data folder (contains 01, 02, 03, ...)');
    if isequal(picked,0)
        error('Stage 6A: no data root selected.');
    end
    CFG.ROOT = picked;
end
ROOT = char(CFG.ROOT);
validate_root(ROOT);

inDir = fullfile(ROOT,'Stage5B');
assert(isfolder(inDir), ...
    'Stage 6A: %s not found. Run Stage 5B first.', inDir);

F = dir(fullfile(inDir, '**', '*_stage5B_condAgg.mat'));
assert(~isempty(F), 'Stage 6A: no Stage 5B files found under %s', inDir);

fprintf('\n===== Stage 6A: build the group dataset =====\n');
fprintf('Root: %s\n', ROOT);
fprintf('Found %d Stage 5B files.\n', numel(F));

%% ---------------------------
%  (3) Read and index every file
%% ---------------------------
rec = struct('file',{},'participant',{},'condition',{},'condIdx',{},'condAgg',{});

for k = 1:numel(F)
    fp = fullfile(F(k).folder, F(k).name);
    L  = load(fp);
    assert(isfield(L,'condAgg'), 'File missing condAgg: %s', fp);
    A = L.condAgg;

    assert(isfield(A,'meta') && isfield(A.meta,'participantName') && ...
           isfield(A.meta,'conditionName'), 'File missing meta labels: %s', fp);

    pName = char(strtrim(string(A.meta.participantName)));
    cName = upper(char(strtrim(string(A.meta.conditionName))));

    ci = find(strcmp(requiredConditions, cName), 1);
    assert(~isempty(ci), 'Unrecognised condition "%s" in %s', cName, fp);

    rec(end+1) = struct('file', fp, 'participant', pName, ...
        'condition', cName, 'condIdx', ci, 'condAgg', A); %#ok<SAGROW>
end

participants = sort(unique({rec.participant}));
nSub = numel(participants);

fprintf('Participants: %d (%s)\n', nSub, strjoin(participants, ', '));

% Every participant must have all four conditions, exactly once. A missing or
% duplicated file would silently misalign the group arrays, so this fails hard.
for i = 1:nSub
    for c = 1:nCond
        n = nnz(strcmp({rec.participant}, participants{i}) & [rec.condIdx] == c);
        assert(n == 1, 'Participant %s has %d files for %s, expected 1.', ...
            participants{i}, n, requiredConditions{c});
    end
end
assert(numel(rec) == nSub*nCond, ...
    'Expected %d files, found %d.', nSub*nCond, numel(rec));

%% ---------------------------
%  (4) Template and consistency
%% ---------------------------
tmpl = rec(1).condAgg;

if isfield(tmpl,'params') && isfield(tmpl.params,'nTime')
    nT = tmpl.params.nTime;
else
    nT = numel(tmpl.curve.all.torque_101);
end

% Channel count taken from the data rather than assumed
if isfield(tmpl,'params') && isfield(tmpl.params,'nBipolar')
    nB = tmpl.params.nBipolar;
elseif isfield(tmpl,'mapInfo') && isfield(tmpl.mapInfo,'bipolar_rowcol')
    nB = size(tmpl.mapInfo.bipolar_rowcol,1);
else
    nB = 28;
end

if isfield(tmpl.params,'time_axis_values')
    xAxis = tmpl.params.time_axis_values;
else
    xAxis = linspace(0,100,nT);
end

for r = 1:numel(rec)
    assert(isfield(rec(r).condAgg,'params') && ...
           rec(r).condAgg.params.nTime == nT, ...
           'nTime mismatch in %s', rec(r).file);
end

%% ---------------------------
%  (5) Build the structure
%% ---------------------------
groupData = struct();

groupData.meta = struct( ...
    'createdOn',     char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'script',        mfilename, ...
    'nParticipants', nSub, ...
    'nConditions',   nCond, ...
    'nFilesUsed',    numel(rec), ...
    'qc_note',       ['No amplitude clipping applied. Isolated high %MVA ' ...
                      'values are retained and reduced through set, condition ' ...
                      'and group averaging.']);

groupData.design = struct( ...
    'factor1_name',   'contractionType', ...
    'factor1_levels', {{'CON','ECC'}}, ...
    'factor2_name',   'intensity', ...
    'factor2_levels', {{'75','90'}}, ...
    'conditionOrder', {requiredConditions});

groupData.participants    = participants(:);     % cell of char, survives -v7.3
groupData.conditions      = requiredConditions(:);
groupData.nSub            = nSub;
groupData.nCond           = nCond;
groupData.nTime           = nT;
groupData.nBipolar        = nB;
groupData.xAxis           = xAxis;
groupData.xAxis_label     = '%CV';
groupData.meanAngle_label = 'Angle (deg)';

if isfield(tmpl,'mapInfo'),         groupData.mapInfo         = tmpl.mapInfo;         else, groupData.mapInfo         = struct(); end
if isfield(tmpl,'coords'),          groupData.coords          = tmpl.coords;          else, groupData.coords          = struct(); end
if isfield(tmpl,'centroid_params'), groupData.centroid_params = tmpl.centroid_params; else, groupData.centroid_params = struct(); end

groupData.sourceFiles = cell(nSub, nCond);

for b = 1:numel(binsToUse)
    bin = binsToUse{b};

    % Full-cycle means, kept for reference
    groupData.scalar.(bin).meanTorque    = nan(nSub, nCond);
    groupData.scalar.(bin).meanAngle     = nan(nSub, nCond);
    groupData.scalar.(bin).meanEMG       = nan(nSub, nCond);
    groupData.scalar.(bin).meanCentroidX = nan(nSub, nCond);
    groupData.scalar.(bin).meanCentroidY = nan(nSub, nCond);

    % Means over the analysis window, which the statistics use
    groupData.scalarCVP.(bin).meanTorque    = nan(nSub, nCond);
    groupData.scalarCVP.(bin).meanEMG       = nan(nSub, nCond);
    groupData.scalarCVP.(bin).meanCentroidX = nan(nSub, nCond);
    groupData.scalarCVP.(bin).meanCentroidY = nan(nSub, nCond);

    groupData.timeseries.(bin).torque_101      = nan(nSub, nCond, nT);
    groupData.timeseries.(bin).angle_101       = nan(nSub, nCond, nT);
    groupData.timeseries.(bin).emgMeanNorm_101 = nan(nSub, nCond, nT);
    groupData.timeseries.(bin).centroidXc_101  = nan(nSub, nCond, nT);
    groupData.timeseries.(bin).centroidYc_101  = nan(nSub, nCond, nT);

    groupData.bipolar.(bin).emgRMS_bipNorm_101 = nan(nSub, nCond, nB, nT);
    groupData.map.(bin).emgRMS_bipNorm_map     = nan(nSub, nCond, nB);
end

groupData.romCV = nan(nSub, nCond);

%% ---------------------------
%  (6) Fill
%% ---------------------------
inCVP    = xAxis >= CFG.cvp_lo & xAxis <= CFG.cvp_hi;
reqCurve = {'torque_101','angle_101','emgMeanNorm_101','centroidXc_101','centroidYc_101'};

for r = 1:numel(rec)
    p = find(strcmp(participants, rec(r).participant), 1);
    c = rec(r).condIdx;
    A = rec(r).condAgg;

    groupData.sourceFiles{p,c} = rec(r).file;

    for b = 1:numel(binsToUse)
        bin = binsToUse{b};

        assert(isfield(A,'curve') && isfield(A.curve,bin), ...
            'Missing curve.%s in %s', bin, rec(r).file);

        for f = 1:numel(reqCurve)
            fn = reqCurve{f};
            assert(isfield(A.curve.(bin),fn) && ~isempty(A.curve.(bin).(fn)) && ...
                   numel(A.curve.(bin).(fn)) == nT, ...
                   'Bad curve.%s.%s in %s', bin, fn, rec(r).file);
        end

        bipTS  = A.bipolar.(bin).emgRMS_bipNorm_101;
        mapVec = A.map.(bin).emgRMS_bipNorm_map(:);
        assert(isequal(size(bipTS),[nB nT]), 'Bad bipolar size in %s', rec(r).file);
        assert(numel(mapVec) == nB,          'Bad map size in %s',     rec(r).file);

        tq = reshape(A.curve.(bin).torque_101,      1, []);
        an = reshape(A.curve.(bin).angle_101,       1, []);
        em = reshape(A.curve.(bin).emgMeanNorm_101, 1, []);
        cx = reshape(A.curve.(bin).centroidXc_101,  1, []);
        cy = reshape(A.curve.(bin).centroidYc_101,  1, []);

        groupData.timeseries.(bin).torque_101(p,c,:)      = tq;
        groupData.timeseries.(bin).angle_101(p,c,:)       = an;
        groupData.timeseries.(bin).emgMeanNorm_101(p,c,:) = em;
        groupData.timeseries.(bin).centroidXc_101(p,c,:)  = cx;
        groupData.timeseries.(bin).centroidYc_101(p,c,:)  = cy;

        groupData.bipolar.(bin).emgRMS_bipNorm_101(p,c,:,:) = bipTS;
        groupData.map.(bin).emgRMS_bipNorm_map(p,c,:)       = mapVec;

        groupData.scalar.(bin).meanTorque(p,c)    = mean(tq,'omitnan');
        groupData.scalar.(bin).meanAngle(p,c)     = mean(an,'omitnan');
        groupData.scalar.(bin).meanEMG(p,c)       = mean(em,'omitnan');
        groupData.scalar.(bin).meanCentroidX(p,c) = mean(cx,'omitnan');
        groupData.scalar.(bin).meanCentroidY(p,c) = mean(cy,'omitnan');

        groupData.scalarCVP.(bin).meanTorque(p,c)    = mean(tq(inCVP),'omitnan');
        groupData.scalarCVP.(bin).meanEMG(p,c)       = mean(em(inCVP),'omitnan');
        groupData.scalarCVP.(bin).meanCentroidX(p,c) = mean(cx(inCVP),'omitnan');
        groupData.scalarCVP.(bin).meanCentroidY(p,c) = mean(cy(inCVP),'omitnan');
    end

    % Angular excursion over the whole constant velocity phase
    anAll = reshape(A.curve.all.angle_101, 1, []);
    groupData.romCV(p,c) = abs(anAll(end) - anAll(1));
end

groupData.meanAngle_101 = squeeze(mean(mean( ...
    groupData.timeseries.all.angle_101,1,'omitnan'),2,'omitnan')).';

groupData.params = struct('cvp_lo', CFG.cvp_lo, 'cvp_hi', CFG.cvp_hi, ...
    'cvp_note', ["scalarCVP holds means over the analysis window; " ...
                 "scalar holds full-cycle means"]);

%% ---------------------------
%  (7) Size checks
%% ---------------------------
for b = 1:numel(binsToUse)
    bin = binsToUse{b};
    assert(isequal(size(groupData.scalar.(bin).meanEMG), [nSub nCond]));
    assert(isequal(size(groupData.timeseries.(bin).emgMeanNorm_101), [nSub nCond nT]));
    assert(isequal(size(groupData.bipolar.(bin).emgRMS_bipNorm_101), [nSub nCond nB nT]));
    assert(isequal(size(groupData.map.(bin).emgRMS_bipNorm_map),     [nSub nCond nB]));
end
assert(~any(cellfun(@isempty, groupData.sourceFiles(:))), ...
    'Some participant and condition cells were never filled.');

%% ---------------------------
%  (8) Report
%% ---------------------------
fprintf('\n===== STAGE 6A GROUP DATA =====\n');
fprintf('%d participants x %d conditions x %d time points x %d channels\n', ...
    nSub, nCond, nT, nB);
fprintf('Analysis window: %d-%d%% CV\n\n', CFG.cvp_lo, CFG.cvp_hi);

S = groupData.scalarCVP.all;
fprintf('Condition means over the analysis window:\n');
fprintf('  %-8s %13s %13s %13s %13s\n','cond','EMG (%MVA)','X (mm)','Y (mm)','torque (Nm)');
for c = 1:nCond
    fprintf('  %-8s %6.1f +/-%5.1f %6.2f +/-%5.2f %6.2f +/-%5.2f %6.1f +/-%4.1f\n', ...
        requiredConditions{c}, ...
        mean(S.meanEMG(:,c)),       std(S.meanEMG(:,c)), ...
        mean(S.meanCentroidX(:,c)), std(S.meanCentroidX(:,c)), ...
        mean(S.meanCentroidY(:,c)), std(S.meanCentroidY(:,c)), ...
        mean(S.meanTorque(:,c)),    std(S.meanTorque(:,c)));
end

% Reported because a difference in angular excursion between modes would
% confound any comparison of mechanical work.
fprintf('\nCV ROM by condition:\n');
for c = 1:nCond
    fprintf('  %-8s %5.2f +/- %4.2f deg\n', requiredConditions{c}, ...
        mean(groupData.romCV(:,c)), std(groupData.romCV(:,c)));
end

%% ---------------------------
%  (9) Save
%% ---------------------------
outName = fullfile(ROOT, 'groupData_stage6A.mat');
save(outName, 'groupData', '-v7.3');

% Flat table, for checking the arrays against the earlier stage summaries
chk = table('Size',[nSub*nCond 8], ...
    'VariableTypes',{'string','string','double','double','double','double','double','double'}, ...
    'VariableNames',{'participant','condition','meanEMG','meanCentroidX', ...
                     'meanCentroidY','meanTorque','romCV','dEMG_earlyLate'});
i = 0;
for p = 1:nSub
    for c = 1:nCond
        i = i + 1;
        chk(i,:) = {string(participants{p}), string(requiredConditions{c}), ...
            S.meanEMG(p,c), S.meanCentroidX(p,c), S.meanCentroidY(p,c), ...
            S.meanTorque(p,c), groupData.romCV(p,c), ...
            groupData.scalarCVP.late.meanEMG(p,c) - groupData.scalarCVP.early.meanEMG(p,c)};
    end
end
writetable(chk, fullfile(ROOT,'stage6A_check.csv'));

fprintf('\nSaved group data : %s\n', outName);
fprintf('Saved check table: %s\n', fullfile(ROOT,'stage6A_check.csv'));

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

    assert(hasParticipant, ...
        ['The selected folder does not look like the expected data root.\n' ...
         'Expected <ROOT>/<participant>/<CON_75|CON_90|ECC_75|ECC_90>/<Set_N>/\n' ...
         'Selected: %s'], ROOT);
end
