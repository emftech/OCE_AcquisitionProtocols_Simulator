function [Ea_Bkg, Ea_Inc, Ea_Global, Ea_Map, Stats] = Compute_AccuracyError_Masks_v2( ...
    GT_Mask, Inc_Mask, Bkg_Mask, Elastogram_Masked)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 04-06-2026
% Date of last modification (DD/MM/YY): 13-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute_AccuracyError_Masks_v2
%
% Computes Accuracy Error (Ea) between an estimated phase-speed elastogram
% and the ground-truth phase-speed map using separated inclusion and
% background ROI masks.
%
% The reported Ea values are computed globally per ROI:
%
%   Ea_ROI = |mean(c_hat_ROI) - mode(c_GT_ROI)| / mode(c_GT_ROI) * 100
%
% A local pixel-wise Ea map is also computed only for visualization:
%
%   Ea_Map(x,y) = |c_hat(x,y) - c_GT(x,y)| / c_GT(x,y) * 100
%
% This version also supports cylindrical/no-inclusion simulations. If the
% inclusion mask is empty, Ea_Inc is returned as NaN and Ea_Global = Ea_Bkg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Dimension checks ===========================

% Resize the ground-truth speed map if it does not match the elastogram.
% Nearest-neighbor interpolation is used to preserve discrete GT values.
if ~isequal(size(GT_Mask), size(Elastogram_Masked))
    GT_Mask = imresize(GT_Mask, size(Elastogram_Masked), 'nearest');
end

% Resize the inclusion ROI mask if needed.
if ~isequal(size(Inc_Mask), size(Elastogram_Masked))
    Inc_Mask = imresize(logical(Inc_Mask), size(Elastogram_Masked), 'nearest');
end

% Resize the background ROI mask if needed.
if ~isequal(size(Bkg_Mask), size(Elastogram_Masked))
    Bkg_Mask = imresize(logical(Bkg_Mask), size(Elastogram_Masked), 'nearest');
end

% Force numerical type for arithmetic operations.
GT_Mask = double(GT_Mask);
Elastogram_Masked = double(Elastogram_Masked);

% Force mask variables to logical arrays.
Inc_Mask = logical(Inc_Mask);
Bkg_Mask = logical(Bkg_Mask);

%% ===================== Local pixel-wise accuracy map =====================

% Compute local pixel-wise accuracy error.
% This map is only used for spatial visualization/diagnostics.
Ea_Map = abs(Elastogram_Masked - GT_Mask) ./ GT_Mask * 100;

% Remove Inf/NaN values caused by invalid GT or elastogram pixels.
Ea_Map(~isfinite(Ea_Map)) = nan;

%% ======================= Valid ROI definitions ===========================

% Valid background pixels must belong to the background mask and have finite,
% positive GT values and finite estimated elastogram values.
Valid_Bkg = Bkg_Mask & isfinite(Elastogram_Masked) & isfinite(GT_Mask) & GT_Mask > 0;

% Valid inclusion pixels follow the same logic. In cylindrical simulations,
% this mask may be entirely empty.
Valid_Inc = Inc_Mask & isfinite(Elastogram_Masked) & isfinite(GT_Mask) & GT_Mask > 0;

% Warn if the background ROI is empty, since the background metric is
% required in both inclusion and no-inclusion cases.
if nnz(Valid_Bkg) == 0
    warning('No valid background pixels found for global accuracy estimation.');
end

% Detect whether a valid inclusion ROI exists.
hasInclusion = nnz(Valid_Inc) > 0;

%% ===================== Background accuracy computation ===================

% Compute the representative GT background speed using the robust mode.
[cGT_Bkg, ModeInfo_Bkg] = local_mode_continuous(GT_Mask(Valid_Bkg));

% Extract valid estimated background values.
cHat_Bkg_vals = Elastogram_Masked(Valid_Bkg);

% Compute descriptive statistics of the estimated background speed.
cHat_Bkg_mean = mean(cHat_Bkg_vals, 'omitnan');
cHat_Bkg_std = std(cHat_Bkg_vals, 0, 'omitnan');
cHat_Bkg_median = median(cHat_Bkg_vals, 'omitnan');

% Compute global ROI accuracy error for the background.
Ea_Bkg = abs(cHat_Bkg_mean - cGT_Bkg) / cGT_Bkg * 100;

%% ===================== Inclusion accuracy computation ====================

if hasInclusion

    % Compute the representative GT inclusion speed using the robust mode.
    [cGT_Inc, ModeInfo_Inc] = local_mode_continuous(GT_Mask(Valid_Inc));

    % Extract valid estimated inclusion values.
    cHat_Inc_vals = Elastogram_Masked(Valid_Inc);

    % Compute descriptive statistics of the estimated inclusion speed.
    cHat_Inc_mean = mean(cHat_Inc_vals, 'omitnan');
    cHat_Inc_std = std(cHat_Inc_vals, 0, 'omitnan');
    cHat_Inc_median = median(cHat_Inc_vals, 'omitnan');

    % Compute global ROI accuracy error for the inclusion.
    Ea_Inc = abs(cHat_Inc_mean - cGT_Inc) / cGT_Inc * 100;

    % If both ROIs exist, the global Ea is the mean of the two ROI Ea values.
    Ea_Global = mean([Ea_Bkg, Ea_Inc], 'omitnan');

else

    % No valid inclusion ROI exists, e.g., cylindrical homogeneous medium.
    % Return NaN for inclusion-specific outputs.
    cGT_Inc = nan;
    ModeInfo_Inc = struct('N',0,'ModeCount',0,'ModeFraction',nan,'RoundingTolerance',nan);

    cHat_Inc_vals = [];
    cHat_Inc_mean = nan;
    cHat_Inc_std = nan;
    cHat_Inc_median = nan;

    Ea_Inc = nan;

    % In no-inclusion simulations, the global accuracy is the background
    % accuracy because the background is the only valid tissue class.
    Ea_Global = Ea_Bkg;
end

%% ===================== Local diagnostic statistics =======================

% Extract local pixel-wise Ea values inside each valid ROI.
Local_Bkg = Ea_Map(Valid_Bkg);
Local_Inc = Ea_Map(Valid_Inc);

% Define the local global diagnostic mask.
Valid_Global_Local = (Valid_Bkg | Valid_Inc) & isfinite(Ea_Map);
Local_Global = Ea_Map(Valid_Global_Local);

%% ============================== Stats output =============================

% Store all output statistics in a structured format.
Stats = struct();

% Store method definitions.
Stats.Method = 'Global ROI Accuracy Error Estimation';
Stats.LocalMapDefinition = '|c_hat(x,y) - c_GT(x,y)| / c_GT(x,y) * 100';
Stats.GlobalROIDefinition = '|mean(c_hat_ROI) - mode(c_GT_ROI)| / mode(c_GT_ROI) * 100';

% Store whether the experiment contained a valid inclusion ROI.
Stats.HasInclusion = hasInclusion;

% Background ROI statistics.
Stats.Bkg.Ea = Ea_Bkg;
Stats.Bkg.c_hat_mean = cHat_Bkg_mean;
Stats.Bkg.c_hat_std = cHat_Bkg_std;
Stats.Bkg.c_hat_median = cHat_Bkg_median;
Stats.Bkg.c_GT_mode = cGT_Bkg;
Stats.Bkg.N = nnz(Valid_Bkg);
Stats.Bkg.ModeInfo = ModeInfo_Bkg;

% Inclusion ROI statistics.
Stats.Inc.Ea = Ea_Inc;
Stats.Inc.c_hat_mean = cHat_Inc_mean;
Stats.Inc.c_hat_std = cHat_Inc_std;
Stats.Inc.c_hat_median = cHat_Inc_median;
Stats.Inc.c_GT_mode = cGT_Inc;
Stats.Inc.N = nnz(Valid_Inc);
Stats.Inc.ModeInfo = ModeInfo_Inc;

% Global accuracy statistics.
Stats.Global.Ea = Ea_Global;
Stats.Global.Definition = 'mean([Ea_Bkg, Ea_Inc]) if inclusion exists; otherwise Ea_Bkg';
Stats.Global.N = nnz(Valid_Bkg) + nnz(Valid_Inc);

% Local pixel-wise diagnostic summaries.
Stats.LocalPixelwise.Bkg.mean = mean(Local_Bkg, 'omitnan');
Stats.LocalPixelwise.Bkg.std = std(Local_Bkg, 0, 'omitnan');
Stats.LocalPixelwise.Bkg.median = median(Local_Bkg, 'omitnan');
Stats.LocalPixelwise.Bkg.N = nnz(isfinite(Local_Bkg));

Stats.LocalPixelwise.Inc.mean = mean(Local_Inc, 'omitnan');
Stats.LocalPixelwise.Inc.std = std(Local_Inc, 0, 'omitnan');
Stats.LocalPixelwise.Inc.median = median(Local_Inc, 'omitnan');
Stats.LocalPixelwise.Inc.N = nnz(isfinite(Local_Inc));

Stats.LocalPixelwise.Global.mean = mean(Local_Global, 'omitnan');
Stats.LocalPixelwise.Global.std = std(Local_Global, 0, 'omitnan');
Stats.LocalPixelwise.Global.median = median(Local_Global, 'omitnan');
Stats.LocalPixelwise.Global.N = nnz(isfinite(Local_Global));

% Store the exact valid masks used.
Stats.ValidMasks.Valid_Bkg = Valid_Bkg;
Stats.ValidMasks.Valid_Inc = Valid_Inc;
Stats.ValidMasks.Valid_Global_Local = Valid_Global_Local;

%% ============================ Print summary ==============================

fprintf('\n===== Global ROI Accuracy Error Summary =====\n');

fprintf('Background:\n');
fprintf('  c_hat_mean = %.4f m/s | c_GT_mode = %.4f m/s | Ea = %.3f %% | N = %d pixels\n', ...
    Stats.Bkg.c_hat_mean, Stats.Bkg.c_GT_mode, Stats.Bkg.Ea, Stats.Bkg.N);

if hasInclusion
    fprintf('Inclusion:\n');
    fprintf('  c_hat_mean = %.4f m/s | c_GT_mode = %.4f m/s | Ea = %.3f %% | N = %d pixels\n', ...
        Stats.Inc.c_hat_mean, Stats.Inc.c_GT_mode, Stats.Inc.Ea, Stats.Inc.N);
else
    fprintf('Inclusion:\n');
    fprintf('  No inclusion ROI detected. Ea_Inc = NaN.\n');
end

fprintf('Global:\n');
fprintf('  Ea_Global = %.3f %%\n', Stats.Global.Ea);
fprintf('=============================================\n\n');

end

%% ===================================================================== %%
function [modeValue, ModeInfo] = local_mode_continuous(values)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local_mode_continuous
%
% Robust mode estimator for GT speed values.
%
% The GT medium maps usually contain repeated values corresponding to the
% assigned material speeds. A small rounding tolerance is used to avoid
% floating-point noise from splitting identical GT values into many unique
% values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert input to column vector.
values = values(:);

% Keep only finite values.
values = values(isfinite(values));

% If no values exist, return NaN and empty diagnostics.
if isempty(values)
    modeValue = nan;
    ModeInfo = struct();
    ModeInfo.N = 0;
    ModeInfo.ModeCount = 0;
    ModeInfo.ModeFraction = nan;
    ModeInfo.RoundingTolerance = nan;
    return;
end

% Tolerance used to group nearly identical floating-point values.
roundTol = 1e-6;

% Round values to the chosen tolerance.
valuesRounded = round(values / roundTol) * roundTol;

% Extract unique rounded values.
uniqueVals = unique(valuesRounded);

% Count the occurrence of each unique rounded value.
counts = zeros(numel(uniqueVals),1);

for k = 1:numel(uniqueVals)
    counts(k) = nnz(valuesRounded == uniqueVals(k));
end

% The mode is the value with the largest count.
[modeCount, idxMax] = max(counts);

modeValue = uniqueVals(idxMax);

% Store mode diagnostics.
ModeInfo = struct();
ModeInfo.N = numel(values);
ModeInfo.ModeCount = modeCount;
ModeInfo.ModeFraction = modeCount / numel(values);
ModeInfo.RoundingTolerance = roundTol;
ModeInfo.NumUniqueRoundedValues = numel(uniqueVals);

end