function [Bkg_Mask, Inc_Mask, bg_value, params] = Get_Bkg_ROI_Masks(GT_Mask, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 23-04-2026
% Date of last modification (DD/MM/YY): 04-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get_Bkg_ROI_Masks
%
% This function takes a ground-truth phase-speed map, GT_Mask, and creates
% two binary ROI masks:
%
%   1) Bkg_Mask:
%      Background-only ROI mask.
%
%      In this updated version, the background ROI is restricted to a
%      circular support region spanning the valid square grid. This makes
%      the background mask geometrically consistent with the circular
%      support observed in the interpolated MB/OCE elastograms.
%
%      Therefore, instead of being a square background with an inclusion
%      cutout, the background ROI becomes:
%
%          circular support
%          minus inclusion/transition region
%          minus transition dead band
%
%   2) Inc_Mask:
%      Inclusion-only ROI mask.
%
%      The inclusion-mask logic is intentionally preserved from the previous
%      version. It is still obtained from pixels sufficiently far away from
%      the dominant background value, followed by erosion to avoid transition
%      pixels.
%
% IMPORTANT:
%   Bkg_Mask and Inc_Mask are NOT complements of each other.
%
%   A dead band is left between them so that neither mask includes transition
%   pixels produced by spatial smoothing, interpolation, or finite boundary
%   thickness.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%   GT_Mask : 2D numeric matrix
%       Ground-truth phase-speed map.
%
% Optional name-value inputs:
%
%   'BgTolFrac' : scalar
%       Fraction of the GT dynamic range used to define the background core
%       around the estimated dominant background value.
%       Default: 0.05
%
%   'IncTolFrac' : scalar
%       Fraction of the GT dynamic range used to define the inclusion core
%       away from the estimated dominant background value.
%       Must be greater than BgTolFrac.
%       Default: 0.15
%
%   'MarginPx' : nonnegative integer
%       Erosion margin applied to both masks to avoid transition pixels.
%       Default: 2
%
%   'NBins' : integer
%       Number of histogram bins used to estimate the dominant background
%       value.
%       Default: 256
%
%   'UseCircularBackgroundSupport' : logical
%       If true, restrict Bkg_Mask to a circular support region.
%       Default: true
%
%   'CircleRadiusFrac' : scalar in (0,1]
%       Fraction of the maximum inscribed circle radius used for the
%       background support.
%
%       CircleRadiusFrac = 1 gives the largest circle that fits inside the
%       image grid.
%
%       Default: 1.0
%
% Outputs:
%
%   Bkg_Mask : logical matrix
%       Circular-support background ROI mask with inclusion/transition
%       region excluded.
%
%   Inc_Mask : logical matrix
%       Inclusion ROI mask.
%
%   bg_value : scalar
%       Estimated dominant background value.
%
%   params : struct
%       Structure containing all parameters and intermediate values used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse inputs ===============================

p = inputParser;

addParameter(p, 'BgTolFrac', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);

addParameter(p, 'IncTolFrac', 0.15, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);

addParameter(p, 'MarginPx', 2, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && floor(x) == x);

addParameter(p, 'NBins', 256, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 8 && floor(x) == x);

addParameter(p, 'UseCircularBackgroundSupport', true, ...
    @(x) islogical(x) && isscalar(x));

addParameter(p, 'CircleRadiusFrac', 1.0, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);

parse(p, varargin{:});

BgTolFrac = p.Results.BgTolFrac;
IncTolFrac = p.Results.IncTolFrac;
MarginPx = p.Results.MarginPx;
NBins = p.Results.NBins;
UseCircularBackgroundSupport = p.Results.UseCircularBackgroundSupport;
CircleRadiusFrac = p.Results.CircleRadiusFrac;

if IncTolFrac <= BgTolFrac
    error('IncTolFrac must be strictly greater than BgTolFrac.');
end

%% ============================ Validate input =============================

if ~isnumeric(GT_Mask) || ndims(GT_Mask) ~= 2
    error('GT_Mask must be a 2D numeric matrix.');
end

valid_vals = GT_Mask(~isnan(GT_Mask));

if isempty(valid_vals)
    error('GT_Mask contains no valid values.');
end

map_min = min(valid_vals);
map_max = max(valid_vals);
map_range = map_max - map_min;

if map_range <= eps
    error('GT_Mask appears to be constant; cannot separate background and inclusions.');
end

[Nx, Ny] = size(GT_Mask);

%% =================== Estimate dominant background value ==================

[counts, edges] = histcounts(valid_vals, NBins);
[~, idx_max] = max(counts);

bin_low  = edges(idx_max);
bin_high = edges(idx_max + 1);

vals_in_mode_bin = valid_vals(valid_vals >= bin_low & valid_vals <= bin_high);

if isempty(vals_in_mode_bin)
    bg_value = 0.5 * (bin_low + bin_high);
else
    bg_value = median(vals_in_mode_bin);
end

%% ================= Build core masks with transition dead band =============

bg_tol  = BgTolFrac  * map_range;
inc_tol = IncTolFrac * map_range;

dist_to_bg = abs(GT_Mask - bg_value);

% Background core:
% pixels close to the dominant background value.
Bkg_core = dist_to_bg <= bg_tol;

% Inclusion core:
% pixels sufficiently different from the dominant background value.
Inc_core = dist_to_bg >= inc_tol;

% Exclude invalid pixels explicitly.
Bkg_core(isnan(GT_Mask)) = false;
Inc_core(isnan(GT_Mask)) = false;

%% ====================== Preserve original inclusion logic =================
% The inclusion mask is computed exactly as before:
% thresholding away from the background value, then erosion.

if MarginPx > 0

    se = make_disk_kernel(MarginPx);

    Inc_Mask = binary_erode_custom(Inc_core, se);

else

    Inc_Mask = Inc_core;
end

%% ====================== Build circular background support =================
% This is the only intentional conceptual update:
%
%   the background ROI is restricted to a circular region that spans the
%   grid similarly to the valid elastogram support.
%
% The circle is centered at the image center and has radius equal to the
% largest inscribed circle, optionally scaled by CircleRadiusFrac.

if UseCircularBackgroundSupport

    [I, J] = ndgrid(1:Nx, 1:Ny);

    cx = (Nx + 1) / 2;
    cy = (Ny + 1) / 2;

    r_max = 0.5 * min(Nx, Ny);
    r_circle = CircleRadiusFrac * r_max;

    CircularSupport = ((I - cx).^2 + (J - cy).^2) <= r_circle^2;

else

    CircularSupport = true(Nx, Ny);
end

CircularSupport(isnan(GT_Mask)) = false;

%% ====================== Compute background mask ==========================
% The background mask keeps the original background-core logic but restricts
% it to the circular support.
%
% Then, as before, it is eroded to leave a margin away from transitions.

Bkg_core_circular = Bkg_core & CircularSupport;

if MarginPx > 0

    se = make_disk_kernel(MarginPx);

    Bkg_Mask = binary_erode_custom(Bkg_core_circular, se);

else

    Bkg_Mask = Bkg_core_circular;
end

%% ====================== Final safety: no overlap =========================

overlap_mask = Bkg_Mask & Inc_Mask;

if any(overlap_mask(:))
    Bkg_Mask(overlap_mask) = false;
    Inc_Mask(overlap_mask) = false;
end

%% ====================== Return parameters ================================

params = struct();

params.BgTolFrac = BgTolFrac;
params.IncTolFrac = IncTolFrac;
params.MarginPx = MarginPx;
params.NBins = NBins;

params.UseCircularBackgroundSupport = UseCircularBackgroundSupport;
params.CircleRadiusFrac = CircleRadiusFrac;

params.bg_tol = bg_tol;
params.inc_tol = inc_tol;

params.bg_value = bg_value;

params.map_min = map_min;
params.map_max = map_max;
params.map_range = map_range;

params.CircularSupport = CircularSupport;
params.Bkg_core = Bkg_core;
params.Bkg_core_circular = Bkg_core_circular;
params.Inc_core = Inc_core;

end

%% ===================================================================== %%
function se = make_disk_kernel(radius_px)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_disk_kernel
%
% Builds a binary disk-shaped structuring element without requiring the
% Image Processing Toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = radius_px;

[x, y] = ndgrid(-r:r, -r:r);

se = (x.^2 + y.^2) <= r^2;

end

%% ===================================================================== %%
function out = binary_erode_custom(mask, se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% binary_erode_custom
%
% Binary erosion implemented with conv2 so that the Image Processing
% Toolbox is not required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

se = logical(se);

required_sum = sum(se(:));

conv_result = conv2(double(mask), double(se), 'same');

out = (conv_result == required_sum);

end