%% ===== Stage 6C3 MID 20-80%: Spatial Heatmaps in CENTERED GRID COORDINATES =====
% Uses Stage6A groupData
%
% PURPOSE:
%   Plot spatial EMG heatmaps for:
%       CON_75 vs ECC_75
%       CON_90 vs ECC_90
%   using the MID 20-80% CV phase.
%
% INPUTS USED:
%   groupData.bipolar.(bin).emgRMS_bipNorm_101 -> [nSub x 4 x 28 x 101]
%   groupData.timeseries.(bin).centroidXc_101  -> [nSub x 4 x 101]
%   groupData.timeseries.(bin).centroidYc_101  -> [nSub x 4 x 101]
%   groupData.mapInfo.bipolar_rowcol           -> [28 x 2]
%   groupData.coords.x_center_mm
%   groupData.coords.y_center_mm
%
% IMPORTANT:
%   - Heatmap values are already normalized to %MVC in Stage 3
%   - Colorbar remains in %MVC units
%   - Coordinates are plotted in CENTERED grid coordinates:
%         lateral  = negative X
%         medial   = positive X
%         proximal = negative Y
%         distal   = positive Y
%   - Centroid is already stored as centered in mm, so it is plotted directly
%   - MID phase here = 20-80% CV
%
% NOTE:
%   - ECC flipping is NOT needed here because maps are averaged across the
%     full 20-80% window, and averaging is order-independent.

clear; clc; close all;

%% -------------------------------------------------------
% Settings
%% -------------------------------------------------------
binName = 'all';                 % 'all' | 'early' | 'late'
midWin = [20 80];                % target CV window

interpMethod = 'natural';
extrapMethod = 'nearest';

nInterpX = 121;
nInterpY = 241;

% Display range in %MVC
% Keep units as %MVC, but use a tighter display range for stronger contrast
climMVC = [20 80];

% Colormaps
mapColormap  = jet(256);         % stronger, publication-like topography
diffColormap = parula(256);      % for difference maps

% Plot styling
electrodeMarkerSize = 95;
electrodeLineWidth  = 1.4;
centroidMarkerSize  = 16;
centroidLineWidth   = 2.2;

%% -------------------------------------------------------
% Load group dataset
%% -------------------------------------------------------
[f,p] = uigetfile('groupData_stage6A.mat','Select Stage6A groupData file');
if isequal(f,0), return; end

L = load(fullfile(p,f));
assert(isfield(L,'groupData'), 'Selected file does not contain variable "groupData".');
groupData = L.groupData;

fprintf('\n===== Stage 6C3 MID loaded =====\n');
fprintf('Participants: %d\n', groupData.nSub);
fprintf('Conditions: %s\n', strjoin(string(groupData.conditions), ', '));
fprintf('Bin: %s\n', binName);
fprintf('Target CV window: %d-%d %%\n', midWin(1), midWin(2));

%% -------------------------------------------------------
% Condition indices
%% -------------------------------------------------------
condNames = string(groupData.conditions);

idx.CON75 = find(condNames=="CON_75",1);
idx.ECC75 = find(condNames=="ECC_75",1);
idx.CON90 = find(condNames=="CON_90",1);
idx.ECC90 = find(condNames=="ECC_90",1);

assert(~isempty(idx.CON75) && ~isempty(idx.ECC75) && ~isempty(idx.CON90) && ~isempty(idx.ECC90), ...
    'Could not find all required conditions in groupData.conditions.');

%% -------------------------------------------------------
% Required fields check
%% -------------------------------------------------------
assert(isfield(groupData,'bipolar') && isfield(groupData.bipolar,binName), ...
    'groupData.bipolar.%s missing.', binName);
assert(isfield(groupData.bipolar.(binName),'emgRMS_bipNorm_101'), ...
    'groupData.bipolar.%s.emgRMS_bipNorm_101 missing.', binName);

assert(isfield(groupData,'timeseries') && isfield(groupData.timeseries,binName), ...
    'groupData.timeseries.%s missing.', binName);
assert(isfield(groupData.timeseries.(binName),'centroidXc_101') && ...
       isfield(groupData.timeseries.(binName),'centroidYc_101'), ...
    'Centroid timeseries missing in groupData.timeseries.%s.', binName);

assert(isfield(groupData,'coords') && ...
       isfield(groupData.coords,'x_center_mm') && ...
       isfield(groupData.coords,'y_center_mm'), ...
    'groupData.coords.x_center_mm / y_center_mm missing.');

assert(isfield(groupData,'mapInfo') && isfield(groupData.mapInfo,'bipolar_rowcol'), ...
    'groupData.mapInfo.bipolar_rowcol missing.');

assert(isfield(groupData,'xAxis'), 'groupData.xAxis missing.');

%% -------------------------------------------------------
% Mid-phase mask
%% -------------------------------------------------------
xAxis_full = groupData.xAxis(:).';
midMask = xAxis_full >= midWin(1) & xAxis_full <= midWin(2);
assert(any(midMask), 'No samples found inside the requested mid window.');

xAxis_mid = xAxis_full(midMask); %#ok<NASGU>
fprintf('MID mask points: %d / %d\n', nnz(midMask), numel(xAxis_full));

%% -------------------------------------------------------
% Extract true MID 20-80% bipolar maps from [nSub x 4 x 28 x 101]
%% -------------------------------------------------------
bipTS_all = groupData.bipolar.(binName).emgRMS_bipNorm_101;   % [nSub x 4 x 28 x 101]

assert(ndims(bipTS_all)==4, 'bipolar emgRMS_bipNorm_101 must be [nSub x 4 x 28 x nT].');
assert(size(bipTS_all,1)==groupData.nSub, 'Unexpected size in dim 1.');
assert(size(bipTS_all,2)==4, 'Unexpected size in dim 2.');
assert(size(bipTS_all,3)==28, 'Unexpected size in dim 3.');
assert(size(bipTS_all,4)==numel(xAxis_full), 'Unexpected size in time dimension.');

% Average over MID 20-80% time window -> [nSub x 4 x 28]
map28_mid_all = mean(bipTS_all(:,:,:,midMask), 4, 'omitnan');

fprintf('Heatmap source: groupData.bipolar.%s.emgRMS_bipNorm_101 averaged over MID 20-80%%\n', binName);

map28_CON75 = squeeze(mean(map28_mid_all(:,idx.CON75,:), 1, 'omitnan'));
map28_ECC75 = squeeze(mean(map28_mid_all(:,idx.ECC75,:), 1, 'omitnan'));
map28_CON90 = squeeze(mean(map28_mid_all(:,idx.CON90,:), 1, 'omitnan'));
map28_ECC90 = squeeze(mean(map28_mid_all(:,idx.ECC90,:), 1, 'omitnan'));

map28_CON75 = map28_CON75(:);
map28_ECC75 = map28_ECC75(:);
map28_CON90 = map28_CON90(:);
map28_ECC90 = map28_ECC90(:);

%% -------------------------------------------------------
% Extract centered centroid traces and average over MID 20-80%
%% -------------------------------------------------------
Cx = groupData.timeseries.(binName).centroidXc_101;   % [nSub x 4 x 101]
Cy = groupData.timeseries.(binName).centroidYc_101;   % [nSub x 4 x 101]

assert(ndims(Cx)==3 && ndims(Cy)==3, ...
    'centroidXc_101 / centroidYc_101 must be [nSub x nCond x nT].');

assert(size(Cx,1)==groupData.nSub && size(Cx,2)==4, 'Unexpected size for centroidXc_101.');
assert(size(Cy,1)==groupData.nSub && size(Cy,2)==4, 'Unexpected size for centroidYc_101.');
assert(size(Cx,3)==numel(xAxis_full) && size(Cy,3)==numel(xAxis_full), ...
    'Centroid time dimension does not match groupData.xAxis.');

% Subject mean centroid over MID contraction cycle: [nSub x 4]
Cx_mean_mid = squeeze(mean(Cx(:,:,midMask), 3, 'omitnan'));
Cy_mean_mid = squeeze(mean(Cy(:,:,midMask), 3, 'omitnan'));

% Group mean centered centroid positions
xc_CON75 = mean(Cx_mean_mid(:,idx.CON75), 'omitnan');
yc_CON75 = mean(Cy_mean_mid(:,idx.CON75), 'omitnan');

xc_ECC75 = mean(Cx_mean_mid(:,idx.ECC75), 'omitnan');
yc_ECC75 = mean(Cy_mean_mid(:,idx.ECC75), 'omitnan');

xc_CON90 = mean(Cx_mean_mid(:,idx.CON90), 'omitnan');
yc_CON90 = mean(Cy_mean_mid(:,idx.CON90), 'omitnan');

xc_ECC90 = mean(Cx_mean_mid(:,idx.ECC90), 'omitnan');
yc_ECC90 = mean(Cy_mean_mid(:,idx.ECC90), 'omitnan');

%% -------------------------------------------------------
% Centroid statistics (MID 20-80%, centered mm)
%% -------------------------------------------------------
Cx_CON75 = Cx_mean_mid(:,idx.CON75);
Cx_ECC75 = Cx_mean_mid(:,idx.ECC75);
Cy_CON75 = Cy_mean_mid(:,idx.CON75);
Cy_ECC75 = Cy_mean_mid(:,idx.ECC75);

Cx_CON90 = Cx_mean_mid(:,idx.CON90);
Cx_ECC90 = Cx_mean_mid(:,idx.ECC90);
Cy_CON90 = Cy_mean_mid(:,idx.CON90);
Cy_ECC90 = Cy_mean_mid(:,idx.ECC90);

[~,pX75,~,stX75] = ttest(Cx_CON75,Cx_ECC75);
[~,pY75,~,stY75] = ttest(Cy_CON75,Cy_ECC75);

[~,pX90,~,stX90] = ttest(Cx_CON90,Cx_ECC90);
[~,pY90,~,stY90] = ttest(Cy_CON90,Cy_ECC90);

fprintf('\n===== Centroid Results | MID 20-80%% | 75%% MVC =====\n');
fprintf('\nMedial-Lateral (X centered, mm):\n');
fprintf('CON: %.3f ± %.3f mm\n', mean(Cx_CON75,'omitnan'), std(Cx_CON75,'omitnan'));
fprintf('ECC: %.3f ± %.3f mm\n', mean(Cx_ECC75,'omitnan'), std(Cx_ECC75,'omitnan'));
fprintf('t(%d)=%.3f , p=%.5f\n', stX75.df, stX75.tstat, pX75);

fprintf('\nProximal-Distal (Y centered, mm):\n');
fprintf('CON: %.3f ± %.3f mm\n', mean(Cy_CON75,'omitnan'), std(Cy_CON75,'omitnan'));
fprintf('ECC: %.3f ± %.3f mm\n', mean(Cy_ECC75,'omitnan'), std(Cy_ECC75,'omitnan'));
fprintf('t(%d)=%.3f , p=%.5f\n', stY75.df, stY75.tstat, pY75);

fprintf('\n===== Centroid Results | MID 20-80%% | 90%% MVC =====\n');
fprintf('\nMedial-Lateral (X centered, mm):\n');
fprintf('CON: %.3f ± %.3f mm\n', mean(Cx_CON90,'omitnan'), std(Cx_CON90,'omitnan'));
fprintf('ECC: %.3f ± %.3f mm\n', mean(Cx_ECC90,'omitnan'), std(Cx_ECC90,'omitnan'));
fprintf('t(%d)=%.3f , p=%.5f\n', stX90.df, stX90.tstat, pX90);

fprintf('\nProximal-Distal (Y centered, mm):\n');
fprintf('CON: %.3f ± %.3f mm\n', mean(Cy_CON90,'omitnan'), std(Cy_CON90,'omitnan'));
fprintf('ECC: %.3f ± %.3f mm\n', mean(Cy_ECC90,'omitnan'), std(Cy_ECC90,'omitnan'));
fprintf('t(%d)=%.3f , p=%.5f\n', stY90.df, stY90.tstat, pY90);

%% -------------------------------------------------------
% Save centroid subject-level table
%% -------------------------------------------------------
centroidTable = table( ...
    string(groupData.participants), ...
    Cx_CON75, Cy_CON75, Cx_ECC75, Cy_ECC75, ...
    Cx_CON90, Cy_CON90, Cx_ECC90, Cy_ECC90, ...
    'VariableNames', { ...
    'Participant', ...
    'Cx_CON75_mid20to80_mm','Cy_CON75_mid20to80_mm', ...
    'Cx_ECC75_mid20to80_mm','Cy_ECC75_mid20to80_mm', ...
    'Cx_CON90_mid20to80_mm','Cy_CON90_mid20to80_mm', ...
    'Cx_ECC90_mid20to80_mm','Cy_ECC90_mid20to80_mm'});

%% -------------------------------------------------------
% Get bipolar coordinates in CENTERED grid system
%% -------------------------------------------------------
rc = groupData.mapInfo.bipolar_rowcol;   % [28 x 2] = [row col]
assert(size(rc,1)==28 && size(rc,2)==2, 'bipolar_rowcol must be [28 x 2].');

row_bp = rc(:,1);
col_bp = rc(:,2);

% Absolute physical coordinates:
%   x_abs = [0 10 20 30]
%   y_abs = [5 15 25 35 45 55 65]
x_abs = (col_bp - 1) * 10;
y_abs = ((row_bp - 1) * 10) + 5;

% Convert to CENTERED coordinates
x_center_mm = groupData.coords.x_center_mm;   % usually 15
y_center_mm = groupData.coords.y_center_mm;   % usually 35

x_ch = x_abs - x_center_mm;   % -> [-15 -5 5 15]
y_ch = y_abs - y_center_mm;   % -> [-30 -20 -10 0 10 20 30]

%% -------------------------------------------------------
% Interpolation grid in CENTERED coordinates
%% -------------------------------------------------------
xq = linspace(min(x_ch), max(x_ch), nInterpX);
yq = linspace(min(y_ch), max(y_ch), nInterpY);
[Xq, Yq] = meshgrid(xq, yq);

%% -------------------------------------------------------
% Interpolate condition maps
%% -------------------------------------------------------
Z_CON75 = make_interp_map(x_ch, y_ch, map28_CON75, Xq, Yq, interpMethod, extrapMethod);
Z_ECC75 = make_interp_map(x_ch, y_ch, map28_ECC75, Xq, Yq, interpMethod, extrapMethod);
Z_CON90 = make_interp_map(x_ch, y_ch, map28_CON90, Xq, Yq, interpMethod, extrapMethod);
Z_ECC90 = make_interp_map(x_ch, y_ch, map28_ECC90, Xq, Yq, interpMethod, extrapMethod);

Z_DIFF75 = Z_ECC75 - Z_CON75;
Z_DIFF90 = Z_ECC90 - Z_CON90;

diffLim75 = max(abs(Z_DIFF75(:)), [], 'omitnan');
diffLim90 = max(abs(Z_DIFF90(:)), [], 'omitnan');

fprintf('\nDifference-map limits:\n');
fprintf('75%% | max abs diff = %.3f %%MVC\n', diffLim75);
fprintf('90%% | max abs diff = %.3f %%MVC\n', diffLim90);

%% =======================================================
% ===== Plot 75% heatmaps ================================
%% =======================================================
figure('Color','w','Name','Spatial EMG Heatmaps | MID 20-80% | 75% MVC');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
plot_spatial_map_centered(Xq,Yq,Z_CON75,x_ch,y_ch,xc_CON75,yc_CON75,climMVC,'CON 75%', ...
    electrodeMarkerSize,electrodeLineWidth,centroidMarkerSize,centroidLineWidth,mapColormap);

nexttile;
plot_spatial_map_centered(Xq,Yq,Z_ECC75,x_ch,y_ch,xc_ECC75,yc_ECC75,climMVC,'ECC 75%', ...
    electrodeMarkerSize,electrodeLineWidth,centroidMarkerSize,centroidLineWidth,mapColormap);

nexttile;
axis off
colormap(mapColormap)
cb = colorbar('eastoutside');
cb.Label.String = 'RMS %MVC';
cb.FontSize = 11;
cb.TickDirection = 'out';
cb.AxisLocation = 'out';
cb.LineWidth = 1.0;
if climMVC(2) <= 60
    cb.Ticks = 0:10:climMVC(2);
else
    cb.Ticks = 0:20:climMVC(2);
end
caxis(climMVC);
title(sprintf('MID %d-%d%% color scale', midWin(1), midWin(2)), 'FontWeight','bold');

%% =======================================================
% ===== Plot 90% heatmaps ================================
%% =======================================================
figure('Color','w','Name','Spatial EMG Heatmaps | MID 20-80% | 90% MVC');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
plot_spatial_map_centered(Xq,Yq,Z_CON90,x_ch,y_ch,xc_CON90,yc_CON90,climMVC,'CON 90%', ...
    electrodeMarkerSize,electrodeLineWidth,centroidMarkerSize,centroidLineWidth,mapColormap);

nexttile;
plot_spatial_map_centered(Xq,Yq,Z_ECC90,x_ch,y_ch,xc_ECC90,yc_ECC90,climMVC,'ECC 90%', ...
    electrodeMarkerSize,electrodeLineWidth,centroidMarkerSize,centroidLineWidth,mapColormap);

nexttile;
axis off
colormap(mapColormap)
cb = colorbar('eastoutside');
cb.Label.String = 'RMS %MVC';
cb.FontSize = 11;
cb.TickDirection = 'out';
cb.AxisLocation = 'out';
cb.LineWidth = 1.0;
if climMVC(2) <= 60
    cb.Ticks = 0:10:climMVC(2);
else
    cb.Ticks = 0:20:climMVC(2);
end
caxis(climMVC);
title(sprintf('MID %d-%d%% color scale', midWin(1), midWin(2)), 'FontWeight','bold');

%% -------------------------------------------------------
% Optional difference-map figures
%% -------------------------------------------------------
figure('Color','w','Name','Spatial EMG Difference Maps | MID 20-80%');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot_spatial_map_centered(Xq,Yq,Z_DIFF75,x_ch,y_ch,NaN,NaN,[-diffLim75 diffLim75], ...
    'ECC - CON | 75%', electrodeMarkerSize,electrodeLineWidth,centroidMarkerSize,centroidLineWidth,diffColormap);

nexttile;
plot_spatial_map_centered(Xq,Yq,Z_DIFF90,x_ch,y_ch,NaN,NaN,[-diffLim90 diffLim90], ...
    'ECC - CON | 90%', electrodeMarkerSize,electrodeLineWidth,centroidMarkerSize,centroidLineWidth,diffColormap);

%% -------------------------------------------------------
% Save centroid table
%% -------------------------------------------------------
[fn,pp] = uiputfile('CentroidResults_Stage6C3_mid20to80.csv','Save centroid results table');
if ~isequal(fn,0)
    writetable(centroidTable, fullfile(pp,fn));
    fprintf('\n✅ Centroid table saved:\n%s\n', fullfile(pp,fn));
end

fprintf('\n✅ Stage 6C3 MID finished successfully.\n');

%% =======================================================
% LOCAL FUNCTIONS
%% =======================================================
function Z = make_interp_map(x, y, v, Xq, Yq, interpMethod, extrapMethod)
    x = x(:);
    y = y(:);
    v = v(:);

    valid = isfinite(x) & isfinite(y) & isfinite(v);
    x = x(valid);
    y = y(valid);
    v = v(valid);

    F = scatteredInterpolant(x, y, v, interpMethod, extrapMethod);
    Z = F(Xq, Yq);
end

function plot_spatial_map_centered(Xq,Yq,Z,x_ch,y_ch,xc,yc,climVals,titleStr,...
    electrodeMarkerSize,electrodeLineWidth,centroidMarkerSize,centroidLineWidth,mapColormap)

    imagesc(Xq(1,:), Yq(:,1), Z);
    set(gca,'YDir','normal');
    axis equal tight;
    hold on;

    colormap(gca, mapColormap);
    caxis(climVals);

    scatter(x_ch, y_ch, electrodeMarkerSize, 'o', ...
        'MarkerEdgeColor','k', ...
        'MarkerFaceColor','none', ...
        'LineWidth', electrodeLineWidth);

    scatter(0, 0, electrodeMarkerSize, ...
        'o', 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'LineWidth',1.0);

    if isfinite(xc) && isfinite(yc)
        plot(xc, yc, 'w+', 'MarkerSize', centroidMarkerSize, 'LineWidth', centroidLineWidth);
        plot(xc, yc, 'k+', 'MarkerSize', centroidMarkerSize-2, 'LineWidth', 1.2);
    end

    xlabel('Medial-Lateral distance (mm)','FontWeight','bold');
    ylabel('Proximal-Distal distance (mm)','FontWeight','bold');
    title(titleStr,'FontWeight','bold');

    xticks([-15 -5 0 5 15]);
    yticks([-30 -20 -10 0 10 20 30]);

    xlim([min(Xq(:)) max(Xq(:))]);
    ylim([min(Yq(:)) max(Yq(:))]);

    box off;
    set(gca,'TickDir','out','LineWidth',1.2,'FontWeight','bold');
end