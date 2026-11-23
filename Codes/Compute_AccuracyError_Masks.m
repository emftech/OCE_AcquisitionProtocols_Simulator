function [Ea_Bkg, Ea_Inc, Ea_Global, Ea_Map, Stats] = Compute_AccuracyError_Masks( ...
    GT_Mask, Inc_Mask, Bkg_Mask, Elastogram_Masked)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Description:
% Compute Accuracy Error (Ea) between the estimated phase-speed elastogram
% and the ground-truth speed map using separated background and inclusion
% masks.
%
% Ea = |c_hat - c_GT| / c_GT * 100%
%
% Inputs:
%   - GT_Mask             : ground-truth phase speed map [m/s]
%   - Inc_Mask            : binary inclusion mask
%   - Bkg_Mask            : binary background mask
%   - Elastogram_Masked   : estimated phase-speed elastogram [m/s]
%
% Outputs:
%   - Ea_Bkg      : mean accuracy error in background ROI [%]
%   - Ea_Inc      : mean accuracy error in inclusion ROI [%]
%   - Ea_Global   : mean accuracy error over background + inclusion ROIs [%]
%   - Ea_Map      : pixel-wise accuracy error map [%]
%   - Stats       : structure with mean, std, median, and valid pixel counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check dimensions

if ~isequal(size(GT_Mask), size(Elastogram_Masked))
    GT_Mask = imresize(GT_Mask, size(Elastogram_Masked), 'nearest');
end

if ~isequal(size(Inc_Mask), size(Elastogram_Masked))
    Inc_Mask = imresize(logical(Inc_Mask), size(Elastogram_Masked), 'nearest');
end

if ~isequal(size(Bkg_Mask), size(Elastogram_Masked))
    Bkg_Mask = imresize(logical(Bkg_Mask), size(Elastogram_Masked), 'nearest');
end

Inc_Mask = logical(Inc_Mask);
Bkg_Mask = logical(Bkg_Mask);

%% Compute pixel-wise accuracy error map

Ea_Map = abs(Elastogram_Masked - GT_Mask) ./ GT_Mask * 100;

% Remove invalid values
Ea_Map(~isfinite(Ea_Map)) = nan;

%% Define valid ROI pixels

Valid_Bkg = Bkg_Mask & isfinite(Ea_Map);
Valid_Inc = Inc_Mask & isfinite(Ea_Map);

Valid_Global = (Valid_Bkg | Valid_Inc) & isfinite(Ea_Map);

%% Compute regional mean accuracy errors

Ea_Bkg = mean(Ea_Map(Valid_Bkg), 'omitnan');
Ea_Inc = mean(Ea_Map(Valid_Inc), 'omitnan');
Ea_Global = mean(Ea_Map(Valid_Global), 'omitnan');

%% Store additional statistics

Stats = struct();

Stats.Bkg.mean   = Ea_Bkg;
Stats.Bkg.std    = std(Ea_Map(Valid_Bkg), 0, 'omitnan');
Stats.Bkg.median = median(Ea_Map(Valid_Bkg), 'omitnan');
Stats.Bkg.N      = nnz(Valid_Bkg);

Stats.Inc.mean   = Ea_Inc;
Stats.Inc.std    = std(Ea_Map(Valid_Inc), 0, 'omitnan');
Stats.Inc.median = median(Ea_Map(Valid_Inc), 'omitnan');
Stats.Inc.N      = nnz(Valid_Inc);

Stats.Global.mean   = Ea_Global;
Stats.Global.std    = std(Ea_Map(Valid_Global), 0, 'omitnan');
Stats.Global.median = median(Ea_Map(Valid_Global), 'omitnan');
Stats.Global.N      = nnz(Valid_Global);

%% Print summary

fprintf('\n===== Accuracy Error Summary =====\n');
fprintf('Background Ea = %.3f %% | std = %.3f %% | N = %d pixels\n', ...
    Stats.Bkg.mean, Stats.Bkg.std, Stats.Bkg.N);
fprintf('Inclusion  Ea = %.3f %% | std = %.3f %% | N = %d pixels\n', ...
    Stats.Inc.mean, Stats.Inc.std, Stats.Inc.N);
fprintf('Global     Ea = %.3f %% | std = %.3f %% | N = %d pixels\n', ...
    Stats.Global.mean, Stats.Global.std, Stats.Global.N);
fprintf('==================================\n\n');

end