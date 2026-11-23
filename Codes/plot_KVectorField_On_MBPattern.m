function plot_KVectorField_On_MBPattern(Kx_crop, Ky_crop, Inc_Mask_Crop, ...
    Number_Meridians, Xcrop_mm, Ycrop_mm, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_KVectorField_On_MBPattern
%
% Plot a sparse, color-coded 2D k-vector field on top of the Multimeridian
% MB-mode sampling pattern.
%
% This function is intended for visual comparison between:
%   1) The estimated 2D phase-gradient direction field, given by Kx_crop
%      and Ky_crop.
%   2) The radial multimeridian sampling geometry used by the MB-mode
%      acquisition protocol.
%
% IMPORTANT ORIENTATION CONVENTION:
%   This function follows the same visual convention used in the corrected
%   elastogram figures:
%
%       horizontal axis = X [mm]
%       vertical axis   = Y [mm]
%
%   and the displayed Y-axis is reversed:
%
%       set(gca,'YDir','reverse')
%
%   so that the figure matches the k-Wave / GT_Mask visualization convention
%   used throughout the simulation pipeline.
%
% IMPORTANT MATRIX CONVENTION:
%   Kx_crop, Ky_crop, and Inc_Mask_Crop are assumed to be stored as:
%
%       [X index x Y index]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%   - Kx_crop (2D matrix, double):
%       Cropped x-component of the estimated local wave vector [rad/m].
%
%   - Ky_crop (2D matrix, double):
%       Cropped y-component of the estimated local wave vector [rad/m].
%
%   - Inc_Mask_Crop (2D logical/double matrix):
%       Cropped binary inclusion mask aligned with Kx_crop and Ky_crop.
%
%   - Number_Meridians (scalar, integer):
%       Number of MB-mode meridians to overlay, e.g. 4, 8, 16, or 32.
%
%   - Xcrop_mm (vector, double):
%       Cropped horizontal/X spatial axis [mm].
%
%   - Ycrop_mm (vector, double):
%       Cropped vertical/Y spatial axis [mm].
%
% Optional input:
%
%   - varargin{1} = opts (struct):
%
%       opts.vectorStride (scalar, integer):
%           Subsampling stride used to display only a few vectors.
%           Higher values produce fewer arrows.
%           default = 28
%
%       opts.vectorScale_mm (scalar, double):
%           Display length of normalized vectors [mm].
%           This only affects visualization.
%           default = 0.55
%
%       opts.forceOutward (logical):
%           If true, vectors are forced to point outward from the specified
%           source origin. This is useful for divergent wave propagation.
%           default = true
%
%       opts.sourceOrigin_mm (1x2 double):
%           Source origin used for outward vector orientation:
%
%               [x0_mm, y0_mm]
%
%           default = [0 0]
%
%       opts.showMBPoints (logical):
%           If true, show gray dots along the multimeridian pattern.
%           default = true
%
%       opts.meridianLineWidth (scalar, double):
%           Line width of the MB meridian traces.
%           default = 1.0
%
%       opts.vectorLineWidth (scalar, double):
%           Line width of the k-vector arrows.
%           default = 2.0
%
%       opts.kClim ([] or 1x2 double):
%           Manual color limits for |k| [rad/m]. If empty, robust color
%           limits are computed from the displayed vectors.
%           default = []
%
%       opts.useRobustColorScale (logical):
%           If true, compute robust |k| color limits using mean ± N std,
%           excluding NaNs/Infs. This avoids outliers dominating the
%           colormap.
%           default = true
%
%       opts.robustNumStd (scalar, double):
%           Number of standard deviations used in the robust color limits.
%           default = 2
%
%       opts.maskLineColor:
%           Color used for the inclusion outline.
%           default = 'k'
%
%       opts.maskLineWidth:
%           Line width used for the inclusion outline.
%           default = 2.5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================== Defaults ===============================

opts = struct();
opts.vectorStride        = 28;
opts.vectorScale_mm      = 0.55;
opts.forceOutward        = true;
opts.sourceOrigin_mm     = [0 0];
opts.showMBPoints        = true;
opts.meridianLineWidth   = 1.0;
opts.vectorLineWidth     = 2.0;
opts.kClim               = [];
opts.useRobustColorScale = true;
opts.robustNumStd        = 2;
opts.maskLineColor       = 'k';
opts.maskLineWidth       = 2.5;

if nargin >= 7 && ~isempty(varargin{1})
    user_opts = varargin{1};
    fn = fieldnames(user_opts);

    for k = 1:numel(fn)
        opts.(fn{k}) = user_opts.(fn{k});
    end
end

%% ============================== Checks =================================

if ~isequal(size(Kx_crop), size(Ky_crop))
    error('Kx_crop and Ky_crop must have the same size.');
end

if ~isequal(size(Kx_crop), size(Inc_Mask_Crop))
    error('Inc_Mask_Crop must have the same size as Kx_crop and Ky_crop.');
end

[nX, nY] = size(Kx_crop);

if numel(Xcrop_mm) ~= nX
    error('Length of Xcrop_mm must match size(Kx_crop,1).');
end

if numel(Ycrop_mm) ~= nY
    error('Length of Ycrop_mm must match size(Kx_crop,2).');
end

Xcrop_mm = double(Xcrop_mm(:));
Ycrop_mm = double(Ycrop_mm(:));

Inc_Mask_Crop = logical(Inc_Mask_Crop);

%% ======================= Build MB Meridian Pattern ======================

x_min = min(Xcrop_mm);
x_max = max(Xcrop_mm);
y_min = min(Ycrop_mm);
y_max = max(Ycrop_mm);

x_center = mean([x_min x_max]);
y_center = mean([y_min y_max]);

Rx = (x_max - x_min) / 2;
Ry = (y_max - y_min) / 2;

r = linspace(-1, 1, nX).';

theta = linspace(0, pi, Number_Meridians + 1);
theta = theta(1:end-1);

rr = repmat(r, 1, Number_Meridians);
tt = repmat(theta, nX, 1);

% Same alternating radial ordering convention used in the acquisition
% vector generation.
rr(:,2:2:end) = rr(end:-1:1, 2:2:end);

X_mer_mm = x_center + (rr .* cos(tt)) * Rx;
Y_mer_mm = y_center + (rr .* sin(tt)) * Ry;

%% ========================== Sparse k-vector field =======================

[Xg_mm, Yg_mm] = ndgrid(Xcrop_mm, Ycrop_mm);

Kmag_crop = hypot(Kx_crop, Ky_crop);

valid = isfinite(Kx_crop) & isfinite(Ky_crop) & ...
        isfinite(Kmag_crop) & (Kmag_crop > 0);

sample_mask = false(size(valid));
sample_mask(1:opts.vectorStride:end, 1:opts.vectorStride:end) = true;

plot_mask = valid & sample_mask;

if nnz(plot_mask) < 3
    warning('Very few vectors selected. Decrease opts.vectorStride if needed.');
end

Xv = Xg_mm(plot_mask);
Yv = Yg_mm(plot_mask);

Kxv = Kx_crop(plot_mask);
Kyv = Ky_crop(plot_mask);
Kmagv = Kmag_crop(plot_mask);

if isempty(Kmagv)
    warning('No valid k-vectors found for plotting.');
    return;
end

% Force outward direction for divergent propagation.
if opts.forceOutward
    rx = Xv - opts.sourceOrigin_mm(1);
    ry = Yv - opts.sourceOrigin_mm(2);

    inward = (Kxv .* rx + Kyv .* ry) < 0;

    Kxv(inward) = -Kxv(inward);
    Kyv(inward) = -Kyv(inward);
end

% Normalize vector direction for visualization. The vector magnitude is
% represented by color only.
Knorm = hypot(Kxv, Kyv) + eps;

Ux = Kxv ./ Knorm;
Uy = Kyv ./ Knorm;

%% ========================== Robust color limits =========================

if ~isempty(opts.kClim)
    kmin = opts.kClim(1);
    kmax = opts.kClim(2);

elseif opts.useRobustColorScale
    k_valid = Kmagv(isfinite(Kmagv));

    k_mean = mean(k_valid, 'omitnan');
    k_std  = std(k_valid,  0, 'omitnan');

    kmin = k_mean - opts.robustNumStd * k_std;
    kmax = k_mean + opts.robustNumStd * k_std;

    kmin = max(kmin, 0);

    if ~isfinite(kmin) || ~isfinite(kmax) || kmax <= kmin
        kmin = min(k_valid);
        kmax = max(k_valid);
    end

else
    kmin = min(Kmagv);
    kmax = max(Kmagv);
end

if abs(kmax - kmin) < eps
    kmax = kmin + eps;
end

%% =============================== Plot ==================================

figure('Name','Sparse k-Vector Field on Multimeridian Pattern', ...
    'Units','normalized', ...
    'Position',[0.10 0.10 0.62 0.72]);

hold on;

% Plot MB meridian pattern first.
for b = 1:Number_Meridians
    plot(X_mer_mm(:,b), Y_mer_mm(:,b), '-', ...
        'Color', [0.72 0.72 0.72], ...
        'LineWidth', opts.meridianLineWidth);

    if opts.showMBPoints
        plot(X_mer_mm(:,b), Y_mer_mm(:,b), '.', ...
            'Color', [0.55 0.55 0.55], ...
            'MarkerSize', 5);
    end
end

% Plot inclusion outline.
plot_mask_outline(Xcrop_mm, Ycrop_mm, Inc_Mask_Crop, ...
    opts.maskLineColor, opts.maskLineWidth);

% Color-coded vectors: color indicates |k|, arrow length is normalized.
cmap = turbo(256);

color_idx = round(1 + 255 * (Kmagv - kmin) / (kmax - kmin));
color_idx = min(max(color_idx, 1), 256);

for i = 1:numel(Kmagv)
    thisColor = cmap(color_idx(i), :);

    quiver(Xv(i), Yv(i), ...
        opts.vectorScale_mm * Ux(i), ...
        opts.vectorScale_mm * Uy(i), ...
        0, ...
        'Color', thisColor, ...
        'LineWidth', opts.vectorLineWidth, ...
        'MaxHeadSize', 1.8);
end

axis image;
set(gca, 'YDir', 'reverse');   % match GT_Mask / elastogram orientation
grid on;

xlabel('X [mm]');
ylabel('Y [mm]');

title({ ...
    'Sparse 2D k-Vector Field on Multimeridian MB Pattern', ...
    sprintf('%d meridians | normalized arrows, color-coded by |k|', Number_Meridians)});

colormap(gca, turbo);
cb = colorbar;
ylabel(cb, '|k| [rad/m]');
clim([kmin kmax]);

hold off;

end

%% ===================================================================== %%
function plot_mask_outline(X_mm, Y_mm, Mask, lineColor, lineWidth)
% Plot the outline of a binary mask over the current axes.
%
% Matrix convention:
%   Mask is stored as [X x Y].
%
% Plotting convention:
%   horizontal axis -> X_mm
%   vertical axis   -> Y_mm
%
% Therefore, the mask is passed directly, matching the previous elastogram
% plotting convention used in the PhaseGradient2D functions.

if nargin < 4 || isempty(lineColor)
    lineColor = 'k';
end

if nargin < 5 || isempty(lineWidth)
    lineWidth = 2.5;
end

Mask = logical(Mask);

if ~any(Mask(:))
    return;
end

contour(X_mm, Y_mm, double(Mask), [0.5 0.5], ...
    'LineColor', lineColor, ...
    'LineWidth', lineWidth, ...
    'LineStyle', '-');

end