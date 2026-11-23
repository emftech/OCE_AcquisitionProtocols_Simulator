function [Ep_Bkg, Ep_Inc, Ep_Global, ROI_Masks, Stats] = Compute_PrecisionError_Masks( ...
    Elastogram_Masked, Inc_Mask, Bkg_Mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 23-04-2026
% Date of last modification (DD/MM/YY): 13-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute_PrecisionError_Masks
%
% Computes Precision Error (Ep) inside background and inclusion ROIs using
% precomputed binary masks.
%
% Precision error is defined as:
%
%   Ep_ROI = sigma_ROI / mu_ROI * 100
%
% where:
%
%   sigma_ROI = standard deviation of the estimated phase speed in the ROI
%   mu_ROI    = mean estimated phase speed in the ROI
%
% This function also supports cylindrical/no-inclusion simulations. If the
% inclusion mask contains no valid pixels, Ep_Inc is returned as NaN and
% Ep_Global = Ep_Bkg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =============================== Input check =============================

% The function requires an elastogram and two ROI masks.
if nargin < 3
    error('Compute_PrecisionError_Masks requires Elastogram_Masked, Inc_Mask, and Bkg_Mask.');
end

%% ============================ Dimension checks ===========================

% Resize inclusion mask if it does not match the elastogram size.
if ~isequal(size(Inc_Mask), size(Elastogram_Masked))
    Inc_Mask = imresize(logical(Inc_Mask), size(Elastogram_Masked), 'nearest');
end

% Resize background mask if it does not match the elastogram size.
if ~isequal(size(Bkg_Mask), size(Elastogram_Masked))
    Bkg_Mask = imresize(logical(Bkg_Mask), size(Elastogram_Masked), 'nearest');
end

% Convert masks to logical arrays.
Inc_Mask = logical(Inc_Mask);
Bkg_Mask = logical(Bkg_Mask);

%% ============================ Valid pixel masks ==========================

% Valid elastogram pixels are finite estimated speed values.
valid_map = isfinite(Elastogram_Masked);

% Background valid mask.
mask_background = Bkg_Mask & valid_map;

% Inclusion valid mask. This may be empty for homogeneous/cylindrical cases.
mask_inclusion  = Inc_Mask & valid_map;

% Global valid mask combines all available ROI pixels.
mask_global     = (mask_background | mask_inclusion) & valid_map;

%% ============================ Extract ROI values =========================

% Extract estimated speeds inside each ROI.
vals_bg  = Elastogram_Masked(mask_background);
vals_inc = Elastogram_Masked(mask_inclusion);
vals_glb = Elastogram_Masked(mask_global);

% Remove any remaining invalid values.
vals_bg  = vals_bg(isfinite(vals_bg));
vals_inc = vals_inc(isfinite(vals_inc));
vals_glb = vals_glb(isfinite(vals_glb));

%% =============================== Safety checks ===========================

% Background is required for all experiments.
if isempty(vals_bg)
    error('Background ROI has no valid pixels after NaN/Inf removal.');
end

% Determine whether an inclusion ROI exists.
hasInclusion = ~isempty(vals_inc);

%% ============================ Background stats ===========================

% Initialize output structure.
Stats = struct();

% Store whether the experiment contains a valid inclusion ROI.
Stats.HasInclusion = hasInclusion;

% Compute background statistics.
Stats.Background.mean   = mean(vals_bg, 'omitnan');
Stats.Background.std    = std(vals_bg, 0, 'omitnan');
Stats.Background.median = median(vals_bg, 'omitnan');
Stats.Background.N      = numel(vals_bg);

% Compute background precision error.
Ep_Bkg = 100 * Stats.Background.std / Stats.Background.mean;

%% ============================ Inclusion stats ============================

if hasInclusion

    % Compute inclusion statistics.
    Stats.Inclusion.mean    = mean(vals_inc, 'omitnan');
    Stats.Inclusion.std     = std(vals_inc, 0, 'omitnan');
    Stats.Inclusion.median  = median(vals_inc, 'omitnan');
    Stats.Inclusion.N       = numel(vals_inc);

    % Compute inclusion precision error.
    Ep_Inc = 100 * Stats.Inclusion.std / Stats.Inclusion.mean;

else

    % No inclusion ROI exists.
    Stats.Inclusion.mean    = nan;
    Stats.Inclusion.std     = nan;
    Stats.Inclusion.median  = nan;
    Stats.Inclusion.N       = 0;

    Ep_Inc = nan;
end

%% ============================== Global stats =============================

% Global mask must contain at least background pixels.
if isempty(vals_glb)
    error('Global ROI has no valid pixels after NaN/Inf removal.');
end

% Compute global descriptive statistics.
Stats.Global.mean       = mean(vals_glb, 'omitnan');
Stats.Global.std        = std(vals_glb, 0, 'omitnan');
Stats.Global.median     = median(vals_glb, 'omitnan');
Stats.Global.N          = numel(vals_glb);

% If inclusion exists, global Ep is computed over the combined ROI.
% If not, global Ep reduces to background Ep.
if hasInclusion
    Ep_Global = 100 * Stats.Global.std / Stats.Global.mean;
else
    Ep_Global = Ep_Bkg;
end

% Store final Ep values.
Stats.Background.Ep = Ep_Bkg;
Stats.Inclusion.Ep  = Ep_Inc;
Stats.Global.Ep     = Ep_Global;

%% ============================== Mask output ==============================

% Store the exact masks used for the precision calculation.
ROI_Masks = struct();
ROI_Masks.background = mask_background;
ROI_Masks.inclusion  = mask_inclusion;
ROI_Masks.global     = mask_global;

%% ============================ Print summary ==============================

fprintf('\n===== Precision Error Summary =====\n');

fprintf('Background Ep = %.3f %% | mean = %.4f m/s | std = %.4f m/s | N = %d pixels\n', ...
    Ep_Bkg, Stats.Background.mean, Stats.Background.std, Stats.Background.N);

if hasInclusion
    fprintf('Inclusion  Ep = %.3f %% | mean = %.4f m/s | std = %.4f m/s | N = %d pixels\n', ...
        Ep_Inc, Stats.Inclusion.mean, Stats.Inclusion.std, Stats.Inclusion.N);
else
    fprintf('Inclusion  Ep = NaN | No inclusion ROI detected.\n');
end

fprintf('Global     Ep = %.3f %% | mean = %.4f m/s | std = %.4f m/s | N = %d pixels\n', ...
    Ep_Global, Stats.Global.mean, Stats.Global.std, Stats.Global.N);

fprintf('===================================\n\n');

end