%% ===== Stage 4: Weighted Centroid Tracking (PHYSICAL COORDS, mm) =====
% Uses Stage-3 output:
%   stage3.rep(r).emgRMS_bipNorm_101  (28 x 101)  -> weights
%   stage3.map.bipolar_rowcol         (28 x 2)    -> [row col] for each bipolar
%
% Coordinate system (documented, reviewer-friendly):
% - Grid has 8 mono rows x 4 cols (IED = 10 mm typical).
% - After bipolar derivation along the proximal–distal axis, there are 7 bipolar rows x 4 cols.
% - Bipolar channel is located at MIDPOINT between two mono electrodes:
%     y(row_bp) = (row_bp-1)*IED + IED/2  -> 5,15,...,65 mm
% - Columns keep their electrode centers:
%     x(col)    = (col-1)*IED            -> 0,10,20,30 mm
%
% Axes:
% - x: lateral -> medial (mm)   (col 1 = lateral at x=0)
% - y: proximal -> distal (mm)  (bp row 1 at y=5 mm, more distal increases y)
%
% Output:
%   stage4.rep(r).centroidX_mm_101, centroidY_mm_101, sumW_101
%   plus optional centered versions relative to grid geometric center.

clc; clear; close all;
%% --- Load Stage 3 file ---
[file3, path3] = uigetfile('*_stage3_EMGnorm.mat', 'Select Stage-3 EMGnorm file');
if isequal(file3,0), return; end
L = load(fullfile(path3, file3));

assert(isfield(L,'stage3'), 'Loaded file missing stage3 struct.');
stage3 = L.stage3;

assert(isfield(stage3,'rep') && isfield(stage3,'map') && isfield(stage3,'params'), ...
    'stage3 missing rep/map/params.');
assert(isfield(stage3.map,'bipolar_rowcol'), 'stage3.map missing bipolar_rowcol.');

rc = stage3.map.bipolar_rowcol;   % bring the location of bp channels as [28 x 2] => [row_bp col]
nB = size(rc,1); % number of bipolar channels

nReps = numel(stage3.rep);
nT = stage3.params.timeNorm_points; % should be 101

%% --- Physical grid parameters ---
IED_mm = 10;   % inter-electrode distance in mm

% Convert bipolar row/col indices -> physical coordinates (mm)
% x: columns: 0,10,20,30 mm
x_i_mm = (rc(:,2) - 1) * IED_mm; % columns 1 - 4, but, because MATLAB indexing starts at 1, we minus 1 here to start at 0

% y: bipolar rows are midpoints between mono rows: 5,15,...,65 mm
y_i_mm = ((rc(:,1) - 1) * IED_mm) + (IED_mm/2); % same here, we subtract 1 to start at 0, but we are adding + (IED/2) = +5 as bp physically locates in between monos 

% Geometric center (for displacement reporting)
x_center_mm = mean([min(x_i_mm) max(x_i_mm)]);   % e.g., 15 mm, x_center = (0 + 30)/2 = 15 mm
y_center_mm = mean([min(y_i_mm) max(y_i_mm)]);   % e.g., 35 mm, y_center = (5 + 65)/2 = 35 mm

%% --- Build Stage 4 output struct ---
stage4 = struct();
stage4.meta   = stage3.meta;
stage4.map    = stage3.map;

stage4.params = struct( ...
    'IED_mm', IED_mm, ...
    'coord_x', 'lateral_to_medial_mm', ...
    'coord_y', 'proximal_to_distal_mm', ...
    'origin_note', 'x=0 at col1 (lateral), y=5mm at bp row1 (between mono row1&2)', ...
    'x_center_mm', x_center_mm, ...
    'y_center_mm', y_center_mm, ...
    'weights', 'emgRMS_bipNorm_101', ...
    'eps_denom', 1e-12 );

stage4.coords = struct( ...
    'x_i_mm', x_i_mm, ...
    'y_i_mm', y_i_mm, ...
    'x_center_mm', x_center_mm, ...
    'y_center_mm', y_center_mm );

stage4.rep = repmat(struct( ...
    'is_valid', false, ...
    'centroidX_mm_101', nan(1,nT), ... % pre-fills with NaNs so missing/invalid reps are obvious
    'centroidY_mm_101', nan(1,nT), ...
    'centroidX_centered_mm_101', nan(1,nT), ...
    'centroidY_centered_mm_101', nan(1,nT), ...
    'sumW_101', nan(1,nT) ...
), nReps, 1);

eps_denom = stage4.params.eps_denom;

%% --- Compute centroid per rep over time ---
for r = 1:nReps
    stage4.rep(r).is_valid = logical(stage3.rep(r).is_valid);

    if ~stage4.rep(r).is_valid
        continue
    end

    W = stage3.rep(r).emgRMS_bipNorm_101;  % [28 x 101], pull the weight matrix
    if isempty(W) || ~ismatrix(W) || size(W,1) ~= nB || size(W,2) ~= nT
        warning('Rep %d: invalid weight matrix size. Expected %dx%d. Skipping.', r, nB, nT);
        stage4.rep(r).is_valid = false;
        continue
    end

    % RMS weights should be non-negative; guard just in case
    W(~isfinite(W)) = NaN; %
    W(W < 0) = 0; % RMS should be ≥ 0 always

    sumW = sum(W, 1, 'omitnan');        % [1 x 101], means Sum the total activation at each time point
    sumW_safe = max(sumW, eps_denom);   % avoid divide-by-zero

    % Weighted centroid in mm - WHERE is the centre of the muscle activity?
    % Cx=distance relative to origin.*amplitude at t devided by the sume of
    % EMG activity for all chs at t, so the mass of the weithed muscle
    % activity at this time point
    cx = (x_i_mm.' * W) ./ sumW_safe;   % [1 x 101], we transpose the x_i' to be [1x28] to fit w's matrix [28x101], then multiplication is valid as 1x101
    cy = (y_i_mm.' * W) ./ sumW_safe;   % [1 x 101], same here

    stage4.rep(r).centroidX_mm_101 = cx;
    stage4.rep(r).centroidY_mm_101 = cy;

    % Centered (displacement relative to geometric grid center)
    stage4.rep(r).centroidX_centered_mm_101 = cx - x_center_mm; % just in case we want to plot geometric SMA
    stage4.rep(r).centroidY_centered_mm_101 = cy - y_center_mm;

    stage4.rep(r).sumW_101 = sumW;
end

%% --- Per-rep QC figures (same layout as cumulative) ---
doPlotEachRep = true;

if doPlotEachRep
    valid = arrayfun(@(s) s.is_valid, stage4.rep);
    repList = find(valid);

    xAxis = linspace(0,100,nT);

    for rr = 1:numel(repList)
        r = repList(rr);

        cx  = stage4.rep(r).centroidX_mm_101;
        cy  = stage4.rep(r).centroidY_mm_101;
        cxc = stage4.rep(r).centroidX_centered_mm_101;
        cyc = stage4.rep(r).centroidY_centered_mm_101;

        figName = sprintf('Centroid Rep %02d (mm)', r);
        figure('Color','w','Name',figName);

        tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

        % --- Top: Absolute ---
        ax1 = nexttile(tl,1); hold(ax1,'on');
        plot(ax1, xAxis, cx, 'LineWidth', 1.8);
        plot(ax1, xAxis, cy, 'LineWidth', 1.8);
        ylabel(ax1,'Centroid (mm)','FontWeight','bold');
        title(ax1, sprintf('Rep %02d | Absolute centroid (mm): X=lateral→medial, Y=proximal→distal', r), ...
            'FontWeight','bold','Interpreter','none');
        legend(ax1, {'X (mm)','Y (mm)'}, 'Location','best', 'Box','off');
        set(ax1,'Box','off','TickDir','out','FontWeight','bold');

        % --- Bottom: Centered ---
        ax2 = nexttile(tl,2); hold(ax2,'on');
        plot(ax2, xAxis, cxc, 'LineWidth', 1.8);
        plot(ax2, xAxis, cyc, 'LineWidth', 1.8);
        xlabel(ax2,'% CV','FontWeight','bold');
        ylabel(ax2,'Centroid displacement (mm)','FontWeight','bold');
        title(ax2, 'Centered centroid relative to geometric grid center', ...
            'FontWeight','bold','Interpreter','none');
        legend(ax2, {'X-centered','Y-centered'}, 'Location','best', 'Box','off');
        set(ax2,'Box','off','TickDir','out','FontWeight','bold');

        ax3 = nexttile(tl,3); hold(ax3,'on');

      meanW = stage4.rep(r).sumW_101 / nB;   % nB = number of bipolar channels (28)

      plot(ax3, xAxis, meanW, 'LineWidth', 1.8);

      xlabel(ax3,'% CV','FontWeight','bold');
      ylabel(ax3,'Mean EMG (%MVC per channel)','FontWeight','bold');

     title(ax3,'Mean EMG activity across grid', ...
     'FontWeight','bold','Interpreter','none');

     set(ax3,'Box','off','TickDir','out','FontWeight','bold');
    end
end

%% --- Save Stage 4 ---
[~, baseName, ~] = fileparts(file3);

% Remove the Stage3 suffix if it exists
baseName = strrep(baseName, '_stage3_EMGnorm', '');

outName = fullfile(path3, sprintf('%s_stage4_centroid_mm.mat', baseName));

save(outName, 'stage4', 'stage3', '-v7.3');

fprintf('\n✅ Saved Stage 4 centroid (mm) output:\n%s\n', outName);



% 
% 
% %% ===== Stage 4 STRUCT CHECK =====
% 
% fprintf('\n===== Stage 4 STRUCT CHECK =====\n')
% 
% % --- Basic dimensions ---
% nReps = numel(stage4.rep);
% nB    = numel(stage4.coords.x_i_mm);
% nT    = stage3.params.timeNorm_points;
% 
% fprintf('Repetitions: %d\n', nReps);
% fprintf('Bipolar channels: %d\n', nB);
% fprintf('Time-normalized samples: %d\n', nT);
% 
% % --- Coordinate checks ---
% fprintf('\n--- Coordinate vectors ---\n');
% disp(size(stage4.coords.x_i_mm))
% disp(size(stage4.coords.y_i_mm))
% 
% fprintf('X range: %.2f → %.2f mm\n', ...
%     min(stage4.coords.x_i_mm), max(stage4.coords.x_i_mm));
% 
% fprintf('Y range: %.2f → %.2f mm\n', ...
%     min(stage4.coords.y_i_mm), max(stage4.coords.y_i_mm));
% 
% fprintf('Grid center: (%.2f , %.2f) mm\n', ...
%     stage4.coords.x_center_mm, stage4.coords.y_center_mm);
% 
% 
% % --- Rep-level checks ---
% validCount = 0;
% 
% for r = 1:nReps
% 
%     if stage4.rep(r).is_valid
% 
%         validCount = validCount + 1;
% 
%         fprintf('\nRep %d\n', r);
% 
%         disp(size(stage4.rep(r).centroidX_mm_101))
%         disp(size(stage4.rep(r).centroidY_mm_101))
%         disp(size(stage4.rep(r).centroidX_centered_mm_101))
%         disp(size(stage4.rep(r).centroidY_centered_mm_101))
%         disp(size(stage4.rep(r).sumW_101))
% 
%         % quick centroid sanity check
%         cx = stage4.rep(r).centroidX_mm_101;
%         cy = stage4.rep(r).centroidY_mm_101;
% 
%         fprintf('X centroid range: %.2f → %.2f mm\n', min(cx), max(cx));
%         fprintf('Y centroid range: %.2f → %.2f mm\n', min(cy), max(cy));
% 
%     end
% 
% end
% 
% fprintf('\nValid reps: %d / %d\n', validCount, nReps);
% 
% 
% fprintf('\n--- sumW QC ---\n')
% for r = 1:nReps
%     if ~stage4.rep(r).is_valid, continue; end
%     sw = stage4.rep(r).sumW_101;
%     fprintf('Rep %d: sumW min=%.2f, median=%.2f, max=%.2f\n', r, ...
%         min(sw,[],'omitnan'), median(sw,'omitnan'), max(sw,[],'omitnan'));
% end