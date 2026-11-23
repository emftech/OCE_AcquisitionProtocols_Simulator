function [Elastogram,R2_map,Elastogram_Masked,Xcrop_mm,Ycrop_mm,...
    Complex_Field_Vz,Inc_Mask_Crop,Bkg_Mask_Crop,Kx_crop,Ky_crop,Kmag_crop] = ...
    PhaseGradient2D_FromSims_v3(WavePropSim_3DMat, Inc_Mask, Bkg_Mask, ...
    kernel_size, Ts_OCT, T_excitation, ...
    horizontalMin, horizontalMax, verticalMin, verticalMax, ...
    applyTemporalWindow, clim_limits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 04-05-2026
% Date of last modification (DD/MM/YY): 04-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PhaseGradient2D_FromSims_v3
%
% Implements the 2D Phase-Gradient Approach for phase-speed estimation from
% simulated wave-based Optical Coherence Elastography data.
%
% The function estimates the local 2D phase-speed field from the complex
% phasor of the simulated wave-propagation video at the excitation
% frequency. For each local spatial window, the wrapped phase is unwrapped
% and fitted to a plane:
%
%       theta(x,y) ≈ kx*x + ky*y + b
%
% where:
%
%       kx = local phase-gradient component along X [rad/m]
%       ky = local phase-gradient component along Y [rad/m]
%
% The local wavenumber magnitude is:
%
%       |k| = sqrt(kx^2 + ky^2)
%
% and the local phase speed is estimated as:
%
%       c = omega / |k|
%
% -------------------------------------------------------------------------
% Mask usage in this version
% -------------------------------------------------------------------------
%
% The function receives:
%
%   1) Inc_Mask:
%      Inclusion ROI mask.
%
%   2) Bkg_Mask:
%      Background ROI mask.
%
% However, Bkg_Mask is NOT used directly as an estimation mask because doing
% so would exclude the inclusion region. Instead, Bkg_Mask and Inc_Mask are
% combined to recover the external circular analysis support:
%
%       AnalysisSupportMask = imfill(Bkg_Mask OR Inc_Mask, 'holes')
%
% This support mask is used only to prevent the phase-gradient estimator
% from running outside the valid circular elastogram domain.
%
% Therefore:
%
%   - outside the circular support:
%       no PG estimation is performed.
%
%   - inside the circular support:
%       PG estimation is allowed, including inside the inclusion.
%
%   - Inc_Mask is still used only for ROI visualization/overlay.
%
% This avoids estimating or displaying the artificial corner regions outside
% the circular support observed in MB/OCE reconstructed elastograms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%   WavePropSim_3DMat:
%       Simulated wave-propagation movie.
%       Expected size:
%
%           [nX x nY x nT]
%
%   Inc_Mask:
%       Inclusion ROI mask.
%       Must have size:
%
%           [nX x nY]
%
%   Bkg_Mask:
%       Background ROI mask.
%       Must have size:
%
%           [nX x nY]
%
%       This mask is used to recover the external circular support, not to
%       exclude the inclusion from the PG estimation.
%
%   kernel_size:
%       Local spatial window size for the phase-gradient estimator.
%       If even, it is automatically incremented by one.
%
%   Ts_OCT:
%       OCT temporal sampling period [s].
%
%   T_excitation:
%       Excitation temporal period [s].
%
%   horizontalMin, horizontalMax:
%       X-axis limits [mm].
%
%   verticalMin, verticalMax:
%       Y-axis limits [mm].
%
%   applyTemporalWindow:
%       true  -> apply Hann window along time before FFT.
%       false -> use raw temporal traces.
%
%   clim_limits:
%       Color limits for elastogram visualization.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%
%   Elastogram:
%       Cropped unmasked phase-speed elastogram [m/s].
%
%   R2_map:
%       Full-size R² map before cropping.
%
%   Elastogram_Masked:
%       Cropped elastogram masked by R² threshold and circular support.
%
%   Xcrop_mm, Ycrop_mm:
%       Cropped spatial axes [mm].
%
%   Complex_Field_Vz:
%       Complex phasor field at excitation frequency.
%
%   Inc_Mask_Crop:
%       Cropped inclusion ROI mask.
%
%   Bkg_Mask_Crop:
%       Cropped background ROI mask.
%
%   Kx_crop, Ky_crop:
%       Cropped phase-gradient components [rad/m].
%
%   Kmag_crop:
%       Cropped wavenumber magnitude [rad/m].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Input defaults =============================

if nargin < 11 || isempty(applyTemporalWindow)
    applyTemporalWindow = false;
end

applyTemporalWindow = logical(applyTemporalWindow);

if nargin < 12 || isempty(clim_limits)
    clim_limits = [4 12];
end

%% ============================ Basic checks ===============================

OCT_fs = 1 / Ts_OCT;

[nX, nY, nT] = size(WavePropSim_3DMat);

if ~isequal(size(Inc_Mask), [nX, nY])
    error(['Inc_Mask must have the same spatial dimensions as ' ...
           'WavePropSim_3DMat, i.e. size(Inc_Mask) = [nX, nY].']);
end

if ~isequal(size(Bkg_Mask), [nX, nY])
    error(['Bkg_Mask must have the same spatial dimensions as ' ...
           'WavePropSim_3DMat, i.e. size(Bkg_Mask) = [nX, nY].']);
end

Inc_Mask = logical(Inc_Mask);
Bkg_Mask = logical(Bkg_Mask);

%% ===================== Analysis support from masks =======================
% Bkg_Mask alone cannot be used as the estimator mask because it intentionally
% excludes the inclusion.
%
% To recover the external circular support, combine the background and
% inclusion masks, then fill internal holes.
%
% This produces the valid computational support where the PG estimator is
% allowed to run.

AnalysisSupportMask = Bkg_Mask | Inc_Mask;

AnalysisSupportMask = imfill(AnalysisSupportMask, 'holes');

% Also reject pixels where the movie is invalid across all time samples.
validMovieSupport = any(isfinite(WavePropSim_3DMat), 3);

AnalysisSupportMask = AnalysisSupportMask & validMovieSupport;

%% ============================ Spatial axes ===============================

Xaxis_mm = linspace(horizontalMin, horizontalMax, nX);
Yaxis_mm = linspace(verticalMin,  verticalMax,  nY);

Xaxis_m  = Xaxis_mm * 1e-3;
Yaxis_m  = Yaxis_mm * 1e-3;

%% ======================== Excitation frequency ===========================

f_exc  = 1 / T_excitation;
omega0 = 2*pi*f_exc;

fprintf('PhaseGradient2D_FromSims_v3()\n');
fprintf('Phase-gradient estimator: f_exc=%.2f Hz, omega0=%.2f rad/s, Ts_OCT=%.3g s (fs=%.2f Hz)\n', ...
    f_exc, omega0, Ts_OCT, OCT_fs);
fprintf('Temporal Hann windowing: %d\n', applyTemporalWindow);

%% ======================= FFT along temporal dimension ====================

if applyTemporalWindow

    hann_win = hanning(nT);

    hann_win1 = zeros(1,1,nT);
    hann_win1(1,1,:) = hann_win;

    hann_win_3dmat = repmat(hann_win1, nX, nY, 1);

    WavePropSim_3DMat_FFTInput = WavePropSim_3DMat .* hann_win_3dmat;

else

    WavePropSim_3DMat_FFTInput = WavePropSim_3DMat;
end

Nfft = 2^10 + 1;

FFT_Vz = fft(WavePropSim_3DMat_FFTInput, Nfft, 3);

freqs = (0:(Nfft-1)) * (OCT_fs/Nfft);

[~, kpos] = min(abs(freqs - f_exc));

f_bin = freqs(kpos);

fprintf('Nearest FFT bin: kpos=%d -> f_bin=%.2f Hz (|f_bin-f_exc|=%.2f Hz)\n', ...
    kpos, f_bin, abs(f_bin - f_exc));

%% ===================== Complex phasor field at f_exc =====================

Complex_Field_Vz = squeeze(FFT_Vz(:,:,kpos));

%% ===================== Sliding-window 2D Phase Gradient ==================

L = kernel_size;

if mod(L,2) == 0
    L = L + 1;
    warning('Window size made odd: L=%d', L);
end

halfL = floor(L/2);

Elastogram_c2D = nan(nX, nY);
Kx_map         = nan(nX, nY);
Ky_map         = nan(nX, nY);
Kmag_map       = nan(nX, nY);
R2_map         = nan(nX, nY);

for ix = (1+halfL):(nX-halfL)

    xw = Xaxis_m((ix-halfL):(ix+halfL));

    for iy = (1+halfL):(nY-halfL)

        support_win = AnalysisSupportMask( ...
            (ix-halfL):(ix+halfL), ...
            (iy-halfL):(iy+halfL));

        % Skip windows that extend outside the valid circular support.
        % This prevents estimation in the artificial corner regions.
        if ~all(support_win(:))
            continue;
        end

        yw = Yaxis_m((iy-halfL):(iy+halfL));

        cwin = Complex_Field_Vz((ix-halfL):(ix+halfL), ...
                                (iy-halfL):(iy+halfL));

        theta = angle(cwin);

        theta_u = unwrap(theta, [], 1);
        theta_u = unwrap(theta_u, [], 2);

        [Yw, Xw] = meshgrid(yw, xw);

        A = [Xw(:), Yw(:), ones(numel(Xw),1)];

        coef = A \ theta_u(:);

        kx = coef(1);
        ky = coef(2);

        theta_hat = A * coef;

        SSres = sum((theta_u(:) - theta_hat).^2);
        SStot = sum((theta_u(:) - mean(theta_u(:))).^2) + eps;

        R2 = 1 - SSres/SStot;

        kmag = hypot(kx, ky);

        if kmag <= eps || ~isfinite(kmag)
            continue;
        end

        c2D = omega0 / kmag;

        Elastogram_c2D(ix, iy) = c2D;
        Kx_map(ix, iy)         = kx;
        Ky_map(ix, iy)         = ky;
        Kmag_map(ix, iy)       = kmag;
        R2_map(ix, iy)         = R2;
    end
end

%% ===================== Enforce support before cropping ===================

Elastogram_c2D(~AnalysisSupportMask) = nan;
R2_map(~AnalysisSupportMask)         = nan;
Kx_map(~AnalysisSupportMask)         = nan;
Ky_map(~AnalysisSupportMask)         = nan;
Kmag_map(~AnalysisSupportMask)       = nan;

%% ===================== Crop to valid kernel region =======================

validX = (1+halfL):(nX-halfL);
validY = (1+halfL):(nY-halfL);

Elastogram = Elastogram_c2D(validX, validY);

Xcrop_mm = Xaxis_mm(validX);
Ycrop_mm = Yaxis_mm(validY);

R2_crop = R2_map(validX, validY);

Inc_Mask_Crop = logical(Inc_Mask(validX, validY));
Bkg_Mask_Crop = logical(Bkg_Mask(validX, validY));

AnalysisSupportMask_Crop = logical(AnalysisSupportMask(validX, validY));

Elastogram(~AnalysisSupportMask_Crop) = nan;

Kx_crop   = Kx_map(validX, validY);
Ky_crop   = Ky_map(validX, validY);
Kmag_crop = Kmag_map(validX, validY);

Kx_crop(~AnalysisSupportMask_Crop)   = nan;
Ky_crop(~AnalysisSupportMask_Crop)   = nan;
Kmag_crop(~AnalysisSupportMask_Crop) = nan;

%% ===================== Display elastogram with inclusion outline =========

figure;
colormap("turbo");
imagesc(Xcrop_mm, Ycrop_mm, Elastogram);
axis image;
xlabel('X [mm]');
ylabel('Y [mm]');
title({sprintf('2D Phase-Gradient Speed Elastogram @ f_{exc}=%.1f Hz', f_exc), ...
    sprintf('Kernel size: L=%dx%d px', L, L)});
colorbar;
clim(clim_limits);
grid off;
hold on;

% Only plot the inclusion outline. The background mask is used internally
% to define the circular computational support.
plot_mask_outline(Xcrop_mm, Ycrop_mm, Inc_Mask_Crop, 'w', 2.5);

hold off;

%% ============================ Fit quality map ============================

figure;
colormap("turbo");
imagesc(Xcrop_mm, Ycrop_mm, R2_crop);
axis image;
xlabel('X [mm]');
ylabel('Y [mm]');
title(sprintf('Plane-fit R^2 (diagnostic), window %dx%d px', L, L));
colorbar;
clim([0 1]);
grid off;
hold on;

plot_mask_outline(Xcrop_mm, Ycrop_mm, Inc_Mask_Crop, 'w', 2.5);

hold off;

%% ============================ Masked elastogram ==========================

R2_thresh = 0.70;

Elastogram_Masked = Elastogram;

Elastogram_Masked(R2_crop < R2_thresh) = nan;
Elastogram_Masked(~AnalysisSupportMask_Crop) = nan;

figure;
colormap("turbo");
imagesc(Xcrop_mm, Ycrop_mm, Elastogram_Masked);
axis image;
xlabel('X [mm]');
ylabel('Y [mm]');
title({'2D Phase-Gradient Masked Elastogram', ...
    sprintf('(R^2 >= %.2f), L=%dx%d px', R2_thresh, L, L)});
colorbar;
clim(clim_limits);
grid off;
hold on;

plot_mask_outline(Xcrop_mm, Ycrop_mm, Inc_Mask_Crop, 'w', 2.5);

hold off;

end

%% ===================================================================== %%
function plot_mask_outline(X_mm, Y_mm, Mask, lineColor, lineWidth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_mask_outline
%
% Plots the contour outline of a binary mask over the current axes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4 || isempty(lineColor)
    lineColor = 'w';
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
    'LineStyle', "-");

end