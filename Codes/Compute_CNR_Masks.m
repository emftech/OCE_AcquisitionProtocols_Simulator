function [CNR_Inc_vs_Bkg, ROI_Masks, Stats] = Compute_CNR_Masks( ...
    Elastogram_Masked, Inc_Mask, Bkg_Mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Description:
% Compute the Contrast-to-Noise Ratio (CNR) between the inclusion ROI and
% the background ROI using precomputed binary masks.
%
% CNR = |mean(background) - mean(ROI)| / sqrt((sigma_background^2 + sigma_ROI^2)/2)
%
% Inputs:
%   - Elastogram_Masked (2D matrix, double):
%       Estimated elastogram on which the CNR metric will be computed.
%       NaN values are omitted from all regional statistics.
%
%   - Inc_Mask (2D logical/double matrix):
%       Binary mask containing all inclusion pixels to be used as ROI.
%
%   - Bkg_Mask (2D logical/double matrix):
%       Binary mask containing background pixels.
%
% Outputs:
%   - CNR_Inc_vs_Bkg (scalar, double):
%       CNR of the inclusion ROI with respect to the background.
%
%   - ROI_Masks (struct):
%       Masks used for the computation:
%           ROI_Masks.background
%           ROI_Masks.inclusion
%
%   - Stats (struct):
%       Regional statistics used to compute CNR.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input checks

if nargin < 3
    error('Compute_CNR_Masks requires Elastogram_Masked, Inc_Mask, and Bkg_Mask.');
end

if ~isequal(size(Inc_Mask), size(Elastogram_Masked))
    Inc_Mask = imresize(logical(Inc_Mask), size(Elastogram_Masked), 'nearest');
end

if ~isequal(size(Bkg_Mask), size(Elastogram_Masked))
    Bkg_Mask = imresize(logical(Bkg_Mask), size(Elastogram_Masked), 'nearest');
end

Inc_Mask = logical(Inc_Mask);
Bkg_Mask = logical(Bkg_Mask);

%% Extract valid regional values

valid_map = isfinite(Elastogram_Masked);

mask_background = Bkg_Mask & valid_map;
mask_inclusion  = Inc_Mask & valid_map;

vals_bg  = Elastogram_Masked(mask_background);
vals_inc = Elastogram_Masked(mask_inclusion);

vals_bg  = vals_bg(isfinite(vals_bg));
vals_inc = vals_inc(isfinite(vals_inc));

%% Safety checks

if isempty(vals_bg)
    error('Background ROI has no valid pixels after NaN/Inf removal.');
end

if isempty(vals_inc)
    error('Inclusion ROI has no valid pixels after NaN/Inf removal.');
end

%% Region statistics

mB = mean(vals_bg, 'omitnan');
sB = std(vals_bg, 0, 'omitnan');

mI = mean(vals_inc, 'omitnan');
sI = std(vals_inc, 0, 'omitnan');

%% Compute CNR

denominator = sqrt((sB^2 + sI^2)/2);

if denominator <= eps
    CNR_Inc_vs_Bkg = NaN;
    warning('CNR denominator is approximately zero. Returning NaN.');
else
    CNR_Inc_vs_Bkg = abs(mB - mI) / denominator;
end

%% Store outputs

ROI_Masks = struct();
ROI_Masks.background = mask_background;
ROI_Masks.inclusion  = mask_inclusion;

Stats = struct();

Stats.Background.mean   = mB;
Stats.Background.std    = sB;
Stats.Background.median = median(vals_bg, 'omitnan');
Stats.Background.N      = numel(vals_bg);

Stats.Inclusion.mean    = mI;
Stats.Inclusion.std     = sI;
Stats.Inclusion.median  = median(vals_inc, 'omitnan');
Stats.Inclusion.N       = numel(vals_inc);

Stats.CNR_Inc_vs_Bkg = CNR_Inc_vs_Bkg;

%% Print summary

fprintf('\n===== CNR Summary =====\n');
fprintf('Background: mean = %.4f m/s | std = %.4f m/s | N = %d pixels\n', ...
    Stats.Background.mean, Stats.Background.std, Stats.Background.N);
fprintf('Inclusion : mean = %.4f m/s | std = %.4f m/s | N = %d pixels\n', ...
    Stats.Inclusion.mean, Stats.Inclusion.std, Stats.Inclusion.N);
fprintf('CNR (Inclusion vs Background) = %.4f\n', CNR_Inc_vs_Bkg);
fprintf('=======================\n\n');

end