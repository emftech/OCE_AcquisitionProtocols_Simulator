function [WavePropSim_3DMat, ValidMask_SBM, InterpOperator_SBM] = ...
    Prepare_SBM_2DPhaseGradient(SampledFrames_SBM, interpMethod, applySmoothing, sigma_px)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 23-06-2026
% Date of last modification (DD/MM/YY): 23-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare_SBM_2DPhaseGradient
%
% Prepare a Synchronized BM Mode (SBM) sampled movie for 2D Phase Gradient
% and 2D PG-TV phase-speed estimation.
%
% The function guarantees that the spatial movie passed to the 2D estimator
% is defined on a square Cartesian grid.
%
% Main behavior:
%
%   1) If the input is already square and at least 100 x 100:
%
%          [nX x nY] with nX == nY and nX >= 100
%
%      no interpolation is performed.
%
%   2) Otherwise, the movie is interpolated to:
%
%          Nout x Nout x Nt
%
%      where:
%
%          Nout = max([nX, nY, 100])
%
%      This preserves the largest available spatial dimension while also
%      enforcing a minimum 100 x 100 grid for fair comparison with the
%      multimeridian MB-mode protocol.
%
% Optional NaN-aware Gaussian smoothing can be applied after interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    SampledFrames_SBM (:,:,:) double
    interpMethod (1,:) char {mustBeMember(interpMethod, ...
        {'linear','nearest','spline','pchip','cubic'})} = 'linear'
    applySmoothing (1,1) logical = false
    sigma_px (1,1) double {mustBePositive} = 0.8
end

%% ============================ Basic sizes ===============================

[nX, nY, nT] = size(SampledFrames_SBM);

% Enforce both:
%   1) square output grid
%   2) minimum 100 x 100 grid for comparison with MB mode
Nout = max([nX, nY, 100]);

fprintf('\nPreparing SBM movie for 2D Phase Gradient...\n');
fprintf('Input SBM movie size       : [%d x %d x %d]\n', nX, nY, nT);
fprintf('Requested output grid size : [%d x %d x %d]\n', Nout, Nout, nT);

%% ======================= Initialize output info ==========================

InterpOperator_SBM = struct();

InterpOperator_SBM.Method = interpMethod;
InterpOperator_SBM.InputSize = [nX, nY, nT];
InterpOperator_SBM.OutputSize = [Nout, Nout, nT];
InterpOperator_SBM.ApplySmoothing = applySmoothing;
InterpOperator_SBM.SigmaPx = sigma_px;
InterpOperator_SBM.MinimumOutputSize = 100;

%% ====================== Case 1: no interpolation needed ==================

if nX == nY && nX >= 100

    fprintf('SBM spatial grid is already square and >= 100 x 100: [%d x %d].\n', nX, nY);
    fprintf('No spatial interpolation needed.\n');

    WavePropSim_3DMat = SampledFrames_SBM;

    ValidMask_SBM = any(isfinite(WavePropSim_3DMat), 3);

    InterpOperator_SBM.InterpolationNeeded = false;
    InterpOperator_SBM.OriginalXIndex = (1:nX).';
    InterpOperator_SBM.OriginalYIndex = (1:nY).';
    InterpOperator_SBM.OutputXIndex = (1:nX).';
    InterpOperator_SBM.OutputYIndex = (1:nY).';

else

    %% ================= Case 2: interpolate to square grid =================

    if nX == nY
        fprintf('SBM spatial grid is square but smaller than 100 x 100: [%d x %d].\n', nX, nY);
        fprintf('Upsampling spatial frames to minimum fair-comparison grid: [%d x %d].\n', Nout, Nout);
    else
        fprintf('SBM spatial grid is rectangular: [%d x %d].\n', nX, nY);
        fprintf('Interpolating spatial frames to square grid: [%d x %d].\n', Nout, Nout);
    end

    x_original = linspace(1, nX, nX).';
    y_original = linspace(1, nY, nY).';

    x_output = linspace(1, nX, Nout).';
    y_output = linspace(1, nY, Nout).';

    [Xq, Yq] = ndgrid(x_output, y_output);

    WavePropSim_3DMat = nan(Nout, Nout, nT);

    for k = 1:nT

        frame_k = SampledFrames_SBM(:,:,k);

        F = griddedInterpolant( ...
            {x_original, y_original}, ...
            frame_k, ...
            interpMethod, ...
            'none');

        WavePropSim_3DMat(:,:,k) = F(Xq, Yq);
    end

    ValidMask_SBM = any(isfinite(WavePropSim_3DMat), 3);

    InterpOperator_SBM.InterpolationNeeded = true;
    InterpOperator_SBM.OriginalXIndex = x_original;
    InterpOperator_SBM.OriginalYIndex = y_original;
    InterpOperator_SBM.OutputXIndex = x_output;
    InterpOperator_SBM.OutputYIndex = y_output;
end

%% ================= Optional NaN-aware spatial smoothing ==================

if applySmoothing

    fprintf('Applying NaN-aware Gaussian smoothing with sigma = %.2f px...\n', sigma_px);

    for k = 1:nT

        frame_k = WavePropSim_3DMat(:,:,k);

        valid_k = isfinite(frame_k);

        if nnz(valid_k) < 10
            continue;
        end

        frame_zero = frame_k;
        frame_zero(~valid_k) = 0;

        weight = double(valid_k);

        frame_blur  = imgaussfilt(frame_zero, sigma_px);
        weight_blur = imgaussfilt(weight, sigma_px);

        frame_smooth = nan(size(frame_k));

        valid_smooth = weight_blur > 1e-8;

        frame_smooth(valid_smooth) = ...
            frame_blur(valid_smooth) ./ weight_blur(valid_smooth);

        WavePropSim_3DMat(:,:,k) = frame_smooth;
    end

    ValidMask_SBM = any(isfinite(WavePropSim_3DMat), 3);
end

%% ============================== Diagnostics ==============================

fprintf('Preparation complete.\n');
fprintf('Output WavePropSim_3DMat size: [%d x %d x %d]\n', ...
    size(WavePropSim_3DMat,1), ...
    size(WavePropSim_3DMat,2), ...
    size(WavePropSim_3DMat,3));

fprintf('Valid prepared pixels: %d / %d\n', ...
    nnz(ValidMask_SBM), numel(ValidMask_SBM));

figure('Name','SBM Preparation for 2D Phase Gradient', ...
    'Units','normalized', ...
    'Position',[0.08 0.10 0.82 0.55]);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% nexttile;
% imagesc(any(isfinite(SampledFrames_SBM),3));
% axis image;
% title('Original Valid Mask');
% xlabel('Y index');
% ylabel('X index');
% colormap(gca, gray);
% colorbar;

% nexttile;
% imagesc(ValidMask_SBM);
% axis image;
% title('Prepared Valid Mask');
% xlabel('Y index');
% ylabel('X index');
% colormap(gca, gray);
% colorbar;

nexttile;
imagesc(SampledFrames_SBM(:,:,1));
axis image;
title('Original Frame 1');
xlabel('Y index');
ylabel('X index');
colormap(gca, turbo);
clim([-1.6 1.6]);
colorbar;

nexttile;
imagesc(WavePropSim_3DMat(:,:,1));
axis image;
title('Prepared Frame 1');
xlabel('Y index');
ylabel('X index');
colormap(gca, turbo);
clim([-1.6 1.6]);
colorbar;

midIdx = round(nT/2);

nexttile;
imagesc(SampledFrames_SBM(:,:,midIdx));
axis image;
title(sprintf('Original Frame %d', midIdx));
xlabel('Y index');
ylabel('X index');
colormap(gca, turbo);
clim([-1.6 1.6]);
colorbar;

nexttile;
imagesc(WavePropSim_3DMat(:,:,midIdx));
axis image;
title(sprintf('Prepared Frame %d', midIdx));
xlabel('Y index');
ylabel('X index');
colormap(gca, turbo);
clim([-1.6 1.6]);
colorbar;

end