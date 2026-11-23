function [WavePropSim_3DMat, ValidMask_MB, HullMask_MB, InterpOperator_MB] = ...
    Prepare_MB_2DPhaseGradient_v2(SampledFrames_MB, interpMethod, fillOutsideHull, ...
                                  applySmoothing, sigma_px)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 24-03-2026
% Last modified: 21-05-2026
%
% Prepare sparse MB-mode sampled movie for 2D Phase Gradient estimation.
%
% This v2 version is intended for the reduced MB grid workflow:
%
%   SampledFrames_MB size = [nAscans x nAscans x Nt]
%
% It uses MATLAB's scatteredInterpolant directly instead of manually
% constructing barycentric interpolation weights.
%
% Recommended method:
%
%   interpMethod = 'nearest'
%
% because reduced multimeridian MB frames are sparse and radial; nearest
% interpolation is more stable and less prone to angular triangular artifacts
% than barycentric interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    SampledFrames_MB (:,:,:) double
    interpMethod (1,:) char {mustBeMember(interpMethod, {'nearest','linear','natural'})} = 'nearest'
    fillOutsideHull (1,1) logical = false
    applySmoothing (1,1) logical = true
    sigma_px (1,1) double {mustBePositive} = 0.8
end

%% ============================ Basic sizes ===============================

[Nx, Ny, Nt] = size(SampledFrames_MB);

fprintf('\nPreparing MB movie for 2D Phase Gradient...\n');
fprintf('Input sparse movie size: [%d x %d x %d]\n', Nx, Ny, Nt);
fprintf('Interpolation method: %s\n', interpMethod);

%% ======================= Build sampled support mask ======================

ValidMask_MB = any(isfinite(SampledFrames_MB), 3);
nValid = nnz(ValidMask_MB);

fprintf('Number of originally sampled pixels: %d\n', nValid);

if nValid < 3
    error('Not enough valid sampled pixels to perform 2D interpolation.');
end

[Xq, Yq] = ndgrid(1:Nx, 1:Ny);

valid_idx = find(ValidMask_MB);
[xv, yv] = ind2sub([Nx, Ny], valid_idx);

xv = double(xv(:));
yv = double(yv(:));

%% ======================== Build convex hull mask =========================
% HullMask_MB records the region supported by the sampled MB geometry.
% Outside this hull, values are kept as NaN unless fillOutsideHull = true.

F_hull = scatteredInterpolant(xv, yv, ones(size(xv)), 'nearest', 'none');
HullMask_MB = isfinite(F_hull(Xq, Yq));

fprintf('Number of pixels inside interpolation hull: %d\n', nnz(HullMask_MB));

if fillOutsideHull
    extrapMethod = 'nearest';
else
    extrapMethod = 'none';
end

%% ====================== Interpolate temporal frames ======================

WavePropSim_3DMat = nan(Nx, Ny, Nt);

for k = 1:Nt

    frame_sparse = SampledFrames_MB(:,:,k);

    v = frame_sparse(valid_idx);
    valid_v = isfinite(v);

    if nnz(valid_v) < 3
        continue;
    end

    F = scatteredInterpolant( ...
        xv(valid_v), ...
        yv(valid_v), ...
        v(valid_v), ...
        interpMethod, ...
        extrapMethod);

    frame_interp = F(Xq, Yq);

    if ~fillOutsideHull
        frame_interp(~HullMask_MB) = nan;
    end

    WavePropSim_3DMat(:,:,k) = frame_interp;
end

%% ==================== Optional spatial Gaussian filter ==================
% NaN-aware Gaussian filtering.
%
% This smooths spatial discontinuities after nearest-neighbor interpolation
% without allowing NaNs/outside-hull pixels to contaminate valid pixels.

if applySmoothing

    fprintf('Applying NaN-aware spatial Gaussian filter, sigma = %.2f px...\n', sigma_px);

    gaussianKernelSize = 2*ceil(3*sigma_px) + 1;

    for k = 1:Nt

        frame_k = WavePropSim_3DMat(:,:,k);
        validMask_k = isfinite(frame_k);

        if nnz(validMask_k) < 10
            continue;
        end

        frame_zero = frame_k;
        frame_zero(~validMask_k) = 0;

        weight = double(validMask_k);

        frame_blur = imgaussfilt( ...
            frame_zero, ...
            sigma_px, ...
            'FilterSize', gaussianKernelSize, ...
            'Padding', 'replicate');

        weight_blur = imgaussfilt( ...
            weight, ...
            sigma_px, ...
            'FilterSize', gaussianKernelSize, ...
            'Padding', 'replicate');

        frame_filtered = nan(size(frame_k));

        validFiltered = weight_blur > 1e-8;

        frame_filtered(validFiltered) = ...
            frame_blur(validFiltered) ./ weight_blur(validFiltered);

        if ~fillOutsideHull
            frame_filtered(~HullMask_MB) = nan;
        end

        WavePropSim_3DMat(:,:,k) = frame_filtered;
    end
end

%% =========================== Output metadata =============================

InterpOperator_MB = struct();

InterpOperator_MB.method          = interpMethod;
InterpOperator_MB.extrapMethod    = extrapMethod;
InterpOperator_MB.Nx              = Nx;
InterpOperator_MB.Ny              = Ny;
InterpOperator_MB.Nt              = Nt;
InterpOperator_MB.valid_idx       = valid_idx;
InterpOperator_MB.xv              = xv;
InterpOperator_MB.yv              = yv;
InterpOperator_MB.ValidMask       = ValidMask_MB;
InterpOperator_MB.HullMask        = HullMask_MB;
InterpOperator_MB.fillOutsideHull = fillOutsideHull;
InterpOperator_MB.applySmoothing  = applySmoothing;
InterpOperator_MB.sigma_px        = sigma_px;

fprintf('Preparation complete.\n');
fprintf('Output WavePropSim_3DMat size: [%d x %d x %d]\n\n', ...
    size(WavePropSim_3DMat,1), ...
    size(WavePropSim_3DMat,2), ...
    size(WavePropSim_3DMat,3));

%% ============================== Diagnostics =============================

midIdx = round(Nt/2);

figure('Name','MB Preparation for 2D Phase Gradient v2', ...
       'Units','normalized', ...
       'Position',[0.05 0.08 0.85 0.60]);

tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(ValidMask_MB);
axis image;
colorbar;
title('Original Sampled Mask');

nexttile;
imagesc(HullMask_MB);
axis image;
colorbar;
title('Interpolation Hull Mask');

nexttile;
imagesc(SampledFrames_MB(:,:,1));
axis image;
colormap(gca, turbo);
colorbar;
title('Sparse Frame 1');

nexttile;
imagesc(WavePropSim_3DMat(:,:,1));
axis image;
colormap(gca, turbo);
colorbar;
title('Prepared Frame 1');

nexttile;
imagesc(SampledFrames_MB(:,:,midIdx));
axis image;
colormap(gca, turbo);
colorbar;
title(sprintf('Sparse Frame %d', midIdx));

nexttile;
imagesc(WavePropSim_3DMat(:,:,midIdx));
axis image;
colormap(gca, turbo);
colorbar;
title(sprintf('Prepared Frame %d', midIdx));

end