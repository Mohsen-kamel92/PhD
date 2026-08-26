%% ===== Pipeline driver =====
%
% PURPOSE
%   Run the analysis stages in order from a single data root, so the whole
%   pipeline can be reproduced without opening each script.
%
% USAGE
%   Open this file and press Run. Choose the data root when prompted, or set
%   PIPE.ROOT below to skip the dialog.
%
%   To run part of the pipeline, edit PIPE.stages. The stage names are listed
%   in full below and any subset can be selected, provided the inputs each
%   stage needs already exist.
%
% WHAT EACH STAGE NEEDS
%   Stages 0 to 5B walk the raw participant folders and need the full dataset.
%   Stage 6A builds groupData_stage6A.mat from the Stage 5B output.
%   Stages 6B onward read only groupData_stage6A.mat, so they run against a
%   repository that ships that file alone.
%
% ON CFG
%   Each stage fills its own settings through set_default, which only writes a
%   field that is absent. Two stages use the same field name for different
%   purposes, so CFG is rebuilt from PIPE.shared before every stage rather than
%   carried forward. Anything placed in PIPE.shared reaches every stage; a
%   setting for one stage only should be edited in that stage.
%
% DEPENDENCIES
%   MATLAB R2023a or later (xregion in the QC figures)
%   Signal Processing Toolbox
%   Statistics and Machine Learning Toolbox
%   spm1d for MATLAB, for Stages 6C1, 6C2 and 6C4 (https://spm1d.org)
%   The functions folder shipped with this repository
%   MAECS_read, for Stage 0 only, available from the system developers

close all;

%% ---------------------------
%  (1) What to run
%% ---------------------------
PIPE = struct();

% Leave empty to be prompted. Set to a path to run without any dialog.
PIPE.ROOT = '';

% Path to spm1d. Leave empty if it is already on the MATLAB path.
PIPE.spm1dPath = '';

% Stages to run, in order. Comment out any that should be skipped.
PIPE.stages = { ...
    'Stage00_sync.m'
    'Stage01_Matching_validity.m'
    'Stage02A_active_torque.m'
    'Stage02B_active_torqueMVC.m'
    'Stage03_bipolar_SignalNormalization.m'
    'Stage04_WeightedCentroid_processing.m'
    'Stage05A_setLevel_aggregation_acrossR....m'
    'Stage05B_ConditionLevel_aggregation_a....m'
    'Stage06A_buildGroup_dataset.m'
    'Stage06B_scalar_RM_ANOVA.m'
    'Stage06C1_SPM_TimeSeries_inference.m'
    'Stage06C2_Early_vs_Late.m'
    'Stage06C3_heat_maps.m'
    'Stage06C4_entropy_spm.m'
    };

% Settings passed to every stage. Anything not listed here is left to each
% stage's own defaults.
PIPE.shared = struct( ...
    'saveFigures', true, ...
    'showTitles',  true);

%% ---------------------------
%  (2) Paths
%% ---------------------------
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'functions'));
addpath(fullfile(thisDir,'scripts'));

if ~isempty(PIPE.spm1dPath)
    addpath(genpath(PIPE.spm1dPath));
end

%% ---------------------------
%  (3) Data root
%% ---------------------------
if isempty(PIPE.ROOT)
    picked = uigetdir(pwd, 'Select the data root folder');
    if isequal(picked,0)
        error('Pipeline: no data root selected.');
    end
    PIPE.ROOT = picked;
end
assert(isfolder(PIPE.ROOT), 'Pipeline: data root not found at %s', PIPE.ROOT);

PIPE.shared.ROOT = PIPE.ROOT;

fprintf('\n========================================\n');
fprintf('  Pipeline\n');
fprintf('  Root   : %s\n', PIPE.ROOT);
fprintf('  Stages : %d\n', numel(PIPE.stages));
fprintf('  Started: %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf('========================================\n');

%% ---------------------------
%  (4) Run
%% ---------------------------
PIPE.log = table('Size',[numel(PIPE.stages) 4], ...
    'VariableTypes',{'string','double','string','string'}, ...
    'VariableNames',{'stage','seconds','status','message'});

for kStage = 1:numel(PIPE.stages)
    stageFile = PIPE.stages{kStage};
    stagePath = fullfile(thisDir,'scripts',stageFile);

    fprintf('\n\n>>> [%d/%d] %s\n', kStage, numel(PIPE.stages), stageFile);

    if ~isfile(stagePath)
        PIPE.log(kStage,:) = {string(stageFile), NaN, "missing", ...
                              string(sprintf('not found at %s', stagePath))};
        fprintf('    SKIPPED: file not found.\n');
        continue
    end

    % Rebuilt each time, so a field one stage defines cannot bind the value
    % another stage would have chosen for the same name.
    CFG = PIPE.shared; %#ok<NASGU>

    tStage = tic;
    try
        run(stagePath);
        PIPE.log(kStage,:) = {string(stageFile), toc(tStage), "ok", ""};
        fprintf('\n    Completed in %.1f s.\n', toc(tStage));
    catch ME
        PIPE.log(kStage,:) = {string(stageFile), toc(tStage), "failed", ...
                              string(ME.message)};
        fprintf('\n    FAILED after %.1f s: %s\n', toc(tStage), ME.message);
        fprintf('    Later stages may depend on this one.\n');
    end

    % Each stage leaves its working variables behind, and some names are reused
    % with different meanings, so the workspace is cleared between stages.
    clearvars -except PIPE thisDir kStage
    close all;
end

%% ---------------------------
%  (5) Report
%% ---------------------------
fprintf('\n\n========================================\n');
fprintf('  Pipeline finished %s\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf('========================================\n\n');

disp(PIPE.log);

nOK = nnz(PIPE.log.status == "ok");
fprintf('\n%d of %d stages completed. Total %.1f min.\n', ...
    nOK, height(PIPE.log), sum(PIPE.log.seconds,'omitnan')/60);

if nOK < height(PIPE.log)
    fprintf('\nStages that did not complete:\n');
    disp(PIPE.log(PIPE.log.status ~= "ok", {'stage','status','message'}));
end

writetable(PIPE.log, fullfile(PIPE.ROOT,'pipeline_log.csv'));
fprintf('\nLog saved to %s\n', fullfile(PIPE.ROOT,'pipeline_log.csv'));
