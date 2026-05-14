%% ===== Stage 6A: Build Group Dataset from Stage5B Files (Sequential per participant) =====
% PURPOSE:
%   Build one group-level dataset for the 2x2 repeated-measures design:
%       Contraction type: CON vs ECC
%       Intensity:        75 vs 90
%
% WORKFLOW:
%   For each participant, user selects 4 Stage5B files in this order:
%       1) CON_75
%       2) ECC_75
%       3) CON_90
%       4) ECC_90
%
% OUTPUT:
%   groupData_stage6A.mat
%
% MAIN OUTPUT CONTENT:
%   groupData.participants
%   groupData.conditions
%   groupData.scalar.(bin).*
%   groupData.timeseries.(bin).*
%
% CONDITIONS (fixed order):
%   1 = CON_75
%   2 = ECC_75
%   3 = CON_90
%   4 = ECC_90

clear; clc;

%% ---------------------------
%  (1) User parameters
%% ---------------------------
expected_nParticipants = 9;  % <-- SET THIS
requiredConditions = {'CON_75','ECC_75','CON_90','ECC_90'};
binsToUse = {'all','early','late'};
nCond = numel(requiredConditions);

%% ---------------------------
%  (2) Load files participant-by-participant
%% ---------------------------
records = struct( ...
    'file', "", ...
    'participant', "", ...
    'condition', "", ...
    'condIdx', NaN, ...
    'condAgg', struct());

recCount = 0;

fprintf('\n===== Stage 6A file loading =====\n');
fprintf('You will now select %d files per participant in this fixed order:\n', nCond);
fprintf('1) CON_75\n2) ECC_75\n3) CON_90\n4) ECC_90\n\n');

for p = 1:expected_nParticipants
    fprintf('\n--- Participant block %d / %d ---\n', p, expected_nParticipants);

    for c = 1:nCond
        condNameExpected = requiredConditions{c};

        [f, pathf] = uigetfile('*_stage5B_condAgg.mat', ...
            sprintf('Participant %d/%d | Select %s Stage5B file', p, expected_nParticipants, condNameExpected));

        if isequal(f,0)
            error('User cancelled while selecting participant %d, condition %s.', p, condNameExpected);
        end

        fullp = fullfile(pathf, f);
        L = load(fullp);
        assert(isfield(L,'condAgg'), 'File missing condAgg: %s', fullp);

        A = L.condAgg;

        % --- participant name
        assert(isfield(A,'meta') && isfield(A.meta,'participantName') && ~isempty(A.meta.participantName), ...
            'Missing meta.participantName in file: %s', fullp);
        participant = string(strtrim(A.meta.participantName));

        % --- condition name
        assert(isfield(A,'meta') && isfield(A.meta,'conditionName') && ~isempty(A.meta.conditionName), ...
            'Missing meta.conditionName in file: %s', fullp);
        condNameFound = upper(string(strtrim(A.meta.conditionName)));

        % --- enforce expected condition
        if condNameFound ~= string(condNameExpected)
            error('Wrong file selected for participant block %d.\nExpected: %s\nFound:    %s\nFile: %s', ...
                p, condNameExpected, condNameFound, fullp);
        end

        recCount = recCount + 1;
        records(recCount).file        = string(fullp);
        records(recCount).participant = participant;
        records(recCount).condition   = condNameFound;
        records(recCount).condIdx     = c;
        records(recCount).condAgg     = A;
    end
end

fprintf('\n✅ Loaded %d participant-condition files.\n', recCount);

%% ---------------------------
%  (3) Participant consistency check within each 4-file block
%% ---------------------------
participantsInOrder = strings(expected_nParticipants,1);

for p = 1:expected_nParticipants
    idx0 = (p-1)*nCond + 1;
    idx1 = p*nCond;

    blockParticipants = string({records(idx0:idx1).participant});
    if numel(unique(blockParticipants)) ~= 1
        error('Participant mismatch inside participant block %d. Selected files do not belong to the same participant.', p);
    end

    participantsInOrder(p) = blockParticipants(1);
end

if numel(unique(participantsInOrder)) ~= expected_nParticipants
    warning('Some participant names are duplicated across blocks. Please verify the selected files.');
end

%% ---------------------------
%  (4) Determine common nTime from first file
%% ---------------------------
tmpl = records(1).condAgg;

if isfield(tmpl,'params') && isfield(tmpl.params,'nTime')
    nT = tmpl.params.nTime;
else
    assert(isfield(tmpl,'curve') && isfield(tmpl.curve,'all') && ...
           isfield(tmpl.curve.all,'torque_101') && ~isempty(tmpl.curve.all.torque_101), ...
           'Cannot determine nTime from template condAgg.');
    nT = numel(tmpl.curve.all.torque_101);
end

xAxis = linspace(0,100,nT);

%% ---------------------------
%  (5) Check all files have same nTime
%% ---------------------------
for r = 1:numel(records)
    A = records(r).condAgg;

    if ~isfield(A,'params') || ~isfield(A.params,'nTime')
        error('File missing params.nTime: %s', records(r).file);
    end

    if A.params.nTime ~= nT
        error('nTime mismatch in file: %s', records(r).file);
    end
end

%% ---------------------------
%  (6) Build groupData structure
%% ---------------------------
nSub = expected_nParticipants;

groupData = struct();

groupData.meta = struct();
groupData.meta.createdOn = datestr(now);
groupData.meta.script = mfilename;
groupData.meta.nParticipants = nSub;
groupData.meta.nConditions = nCond;
groupData.meta.nFilesUsed = numel(records);

groupData.design = struct();
groupData.design.factor1_name = 'contractionType';
groupData.design.factor1_levels = {'CON','ECC'};
groupData.design.factor2_name = 'intensity';
groupData.design.factor2_levels = {'75','90'};
groupData.design.conditionOrder = requiredConditions;

groupData.participants = participantsInOrder(:);
groupData.conditions = requiredConditions(:);
groupData.nSub = nSub;
groupData.nCond = nCond;
groupData.nTime = nT;
groupData.xAxis = xAxis;

% Keep spatial definitions from template
if isfield(tmpl,'mapInfo'), groupData.mapInfo = tmpl.mapInfo; else, groupData.mapInfo = struct(); end
if isfield(tmpl,'coords'), groupData.coords = tmpl.coords; else, groupData.coords = struct(); end
if isfield(tmpl,'centroid_params'), groupData.centroid_params = tmpl.centroid_params; else, groupData.centroid_params = struct(); end

groupData.sourceFiles = strings(nSub, nCond);

for b = 1:numel(binsToUse)
    bin = binsToUse{b};

    % Scalar outcomes: [nSub x nCond]
    groupData.scalar.(bin).meanTorque    = nan(nSub, nCond);
    groupData.scalar.(bin).meanAngle     = nan(nSub, nCond);
    groupData.scalar.(bin).meanEMG       = nan(nSub, nCond);
    groupData.scalar.(bin).meanCentroidX = nan(nSub, nCond);
    groupData.scalar.(bin).meanCentroidY = nan(nSub, nCond);

    % Time-series outcomes: [nSub x nCond x nT]
    groupData.timeseries.(bin).torque_101      = nan(nSub, nCond, nT);
    groupData.timeseries.(bin).angle_101       = nan(nSub, nCond, nT);
    groupData.timeseries.(bin).emgMeanNorm_101 = nan(nSub, nCond, nT);
    groupData.timeseries.(bin).centroidXc_101  = nan(nSub, nCond, nT);
    groupData.timeseries.(bin).centroidYc_101  = nan(nSub, nCond, nT);

    % Spatial channel outcomes
    groupData.bipolar.(bin).emgRMS_bipNorm_101 = nan(nSub, nCond, 28, nT); % [nSub x nCond x 28 x 101]
    groupData.map.(bin).emgRMS_bipNorm_map28   = nan(nSub, nCond, 28);     % [nSub x nCond x 28]
end

%% ---------------------------
%  (7) Fill groupData
%% ---------------------------
for p = 1:nSub
    for c = 1:nCond
        recIdx = (p-1)*nCond + c;
        A = records(recIdx).condAgg;

        groupData.sourceFiles(p,c) = records(recIdx).file;

        for b = 1:numel(binsToUse)
            bin = binsToUse{b};

            assert(isfield(A,'curve') && isfield(A.curve,bin), ...
                'Missing curve.%s in file: %s', bin, records(recIdx).file);

            reqFields = {'torque_101','angle_101','emgMeanNorm_101','centroidXc_101','centroidYc_101'};
            for rf = 1:numel(reqFields)
                fn = reqFields{rf};
                assert(isfield(A.curve.(bin), fn) && ~isempty(A.curve.(bin).(fn)), ...
                    'Missing curve.%s.%s in file: %s', bin, fn, records(recIdx).file);
                assert(numel(A.curve.(bin).(fn)) == nT, ...
                    'Wrong length in curve.%s.%s in file: %s', bin, fn, records(recIdx).file);
            end
                        % --- spatial channel data ---
            assert(isfield(A,'bipolar') && isfield(A.bipolar,bin) && ...
                   isfield(A.bipolar.(bin),'emgRMS_bipNorm_101') && ...
                   ~isempty(A.bipolar.(bin).emgRMS_bipNorm_101), ...
                   'Missing bipolar.%s.emgRMS_bipNorm_101 in file: %s', ...
                   bin, records(recIdx).file);

            assert(isfield(A,'map') && isfield(A.map,bin) && ...
                   isfield(A.map.(bin),'emgRMS_bipNorm_map28') && ...
                   ~isempty(A.map.(bin).emgRMS_bipNorm_map28), ...
                   'Missing map.%s.emgRMS_bipNorm_map28 in file: %s', ...
                   bin, records(recIdx).file);

            bipTS = A.bipolar.(bin).emgRMS_bipNorm_101;   % [28 x nT]
            map28 = A.map.(bin).emgRMS_bipNorm_map28;     % [28 x 1] or [1 x 28]

            assert(ismatrix(bipTS) && all(size(bipTS) == [28 nT]), ...
                   'Wrong size in bipolar.%s.emgRMS_bipNorm_101 in file: %s', ...
                   bin, records(recIdx).file);

            map28 = map28(:);
            assert(numel(map28) == 28, ...
                   'Wrong size in map.%s.emgRMS_bipNorm_map28 in file: %s', ...
                   bin, records(recIdx).file);

            groupData.bipolar.(bin).emgRMS_bipNorm_101(p,c,:,:) = bipTS;
            groupData.map.(bin).emgRMS_bipNorm_map28(p,c,:)     = map28;
            torqueTS = reshape(A.curve.(bin).torque_101,      1, []);
            angleTS  = reshape(A.curve.(bin).angle_101,       1, []);
            emgTS    = reshape(A.curve.(bin).emgMeanNorm_101, 1, []);
            cxTS     = reshape(A.curve.(bin).centroidXc_101,  1, []);
            cyTS     = reshape(A.curve.(bin).centroidYc_101,  1, []);

            groupData.timeseries.(bin).torque_101(p,c,:)      = torqueTS;
            groupData.timeseries.(bin).angle_101(p,c,:)       = angleTS;
            groupData.timeseries.(bin).emgMeanNorm_101(p,c,:) = emgTS;
            groupData.timeseries.(bin).centroidXc_101(p,c,:)  = cxTS;
            groupData.timeseries.(bin).centroidYc_101(p,c,:)  = cyTS;

            groupData.scalar.(bin).meanTorque(p,c)    = mean(torqueTS, 'omitnan');
            groupData.scalar.(bin).meanAngle(p,c)     = mean(angleTS,  'omitnan');
            groupData.scalar.(bin).meanEMG(p,c)       = mean(emgTS,    'omitnan');
            groupData.scalar.(bin).meanCentroidX(p,c) = mean(cxTS,     'omitnan');
            groupData.scalar.(bin).meanCentroidY(p,c) = mean(cyTS,     'omitnan');
        end
    end
end

%% ---------------------------
%  (8) QC summary
%% ---------------------------
fprintf('\n===== STAGE 6A GROUP DATA SUMMARY =====\n');
fprintf('Participants included: %d\n', nSub);
fprintf('Conditions: %s\n', strjoin(string(requiredConditions), ', '));
fprintf('nTime = %d\n', nT);
fprintf('bipolar emgRMS_bipNorm_101 size: '); disp(size(groupData.bipolar.(bin).emgRMS_bipNorm_101))
fprintf('map emgRMS_bipNorm_map28 size:   '); disp(size(groupData.map.(bin).emgRMS_bipNorm_map28))
for b = 1:numel(binsToUse)
    bin = binsToUse{b};

    fprintf('\n--- BIN: %s ---\n', upper(bin));
    fprintf('meanTorque size:    '); disp(size(groupData.scalar.(bin).meanTorque))
    fprintf('meanAngle size:     '); disp(size(groupData.scalar.(bin).meanAngle))
    fprintf('meanEMG size:       '); disp(size(groupData.scalar.(bin).meanEMG))
    fprintf('meanCentroidX size: '); disp(size(groupData.scalar.(bin).meanCentroidX))
    fprintf('meanCentroidY size: '); disp(size(groupData.scalar.(bin).meanCentroidY))

    fprintf('torque_101 size:      '); disp(size(groupData.timeseries.(bin).torque_101))
    fprintf('angle_101 size:       '); disp(size(groupData.timeseries.(bin).angle_101))
    fprintf('emgMeanNorm_101 size: '); disp(size(groupData.timeseries.(bin).emgMeanNorm_101))
    fprintf('centroidXc_101 size:  '); disp(size(groupData.timeseries.(bin).centroidXc_101))
    fprintf('centroidYc_101 size:  '); disp(size(groupData.timeseries.(bin).centroidYc_101))
end

for b = 1:numel(binsToUse)
    bin = binsToUse{b};

    assert(all(size(groupData.scalar.(bin).meanTorque)    == [nSub nCond]));
    assert(all(size(groupData.scalar.(bin).meanAngle)     == [nSub nCond]));
    assert(all(size(groupData.scalar.(bin).meanEMG)       == [nSub nCond]));
    assert(all(size(groupData.scalar.(bin).meanCentroidX) == [nSub nCond]));
    assert(all(size(groupData.scalar.(bin).meanCentroidY) == [nSub nCond]));

    assert(all(size(groupData.timeseries.(bin).torque_101)      == [nSub nCond nT]));
    assert(all(size(groupData.timeseries.(bin).angle_101)       == [nSub nCond nT]));
    assert(all(size(groupData.timeseries.(bin).emgMeanNorm_101) == [nSub nCond nT]));
    assert(all(size(groupData.bipolar.(bin).emgRMS_bipNorm_101) == [nSub nCond 28 nT]));
    assert(all(size(groupData.map.(bin).emgRMS_bipNorm_map28)   == [nSub nCond 28]));
    assert(all(size(groupData.timeseries.(bin).centroidXc_101)  == [nSub nCond nT]));
    assert(all(size(groupData.timeseries.(bin).centroidYc_101)  == [nSub nCond nT]));
end

fprintf('\n✅ Stage 6A groupData built successfully.\n');

%% ---------------------------
%  (9) Save groupData
%% ---------------------------
outFolder = uigetdir(pwd, 'Select output folder for Stage6A groupData');
if isequal(outFolder,0)
    outFolder = pwd;
end

outName = fullfile(outFolder, 'groupData_stage6A.mat');
save(outName, 'groupData', '-v7.3');

fprintf('✅ Saved Stage 6A group dataset:\n%s\n', outName);