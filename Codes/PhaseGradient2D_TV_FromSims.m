function [Elastogram_TV,R2_map,Elastogram_Masked_TV,Xcrop_mm,Ycrop_mm,...
    Complex_Field_Vz,Inc_Mask_Crop,Bkg_Mask_Crop,Kx_crop,Ky_crop,Kmag_crop,...
    Kmag_TV_crop,Elastogram_raw,TVInfo] = ...
    PhaseGradient2D_TV_FromSims(WavePropSim_3DMat, Inc_Mask, Bkg_Mask, ...
    kernel_size, Ts_OCT, T_excitation, ...
    horizontalMin, horizontalMax, verticalMin, verticalMax, ...
    applyTemporalWindow, clim_limits, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 04-05-2026
% Date of last modification (DD/MM/YY): 04-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PhaseGradient2D_TV_FromSims
%
% Regularized 2D Phase-Gradient estimator for simulated wave-based OCE data.
%
% This function follows the same indexing, cropping, plotting convention,
% and local plane-fitting workflow of PhaseGradient2D_FromSims_v3(), but
% adds Total Variation regularization on the cropped raw wavenumber
% magnitude map:
%
%       K_TV = argmin_K 0.5 ||K - K_raw||_2^2 + mu TV(K)
%
% Then:
%
%       c_TV = omega0 ./ K_TV
%
% Required external file:
%
%       TVDenoise_ChambollePock2D.m
%
% -------------------------------------------------------------------------
% Mask behavior
% -------------------------------------------------------------------------
%
% Inputs:
%
%   Inc_Mask:
%       Inclusion ROI mask.
%
%   Bkg_Mask:
%       Background ROI mask.
%
% The background mask is NOT used directly as the estimation mask because it
% would exclude the inclusion. Instead, Bkg_Mask and Inc_Mask are combined
% to recover the external circular computational support:
%
%       AnalysisSupportMask = fill_holes(Bkg_Mask OR Inc_Mask)
%
% Therefore:
%
%   - outside the circular support:
%       the 2D PG estimator does not run.
%
%   - inside the circular support:
%       the estimator is allowed to run, including inside the inclusion.
%
%   - Inc_Mask is still used for ROI visualization and later metrics.
%
%   - Bkg_Mask_Crop is returned for downstream background metrics.
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

%% ============================ Optional parameters ========================

p = inputParser;
p.FunctionName = 'PhaseGradient2D_TV_FromSims';

addParameter(p, 'MuTV', 75);
addParameter(p, 'TVIterations', 250);
addParameter(p, 'TVTau', 0.25);
addParameter(p, 'TVSigma', 0.25);
addParameter(p, 'TVTheta', 1.0);
addParameter(p, 'R2Threshold', 0.70);
addParameter(p, 'MinK', 100);
addParameter(p, 'MaxK', inf);
addParameter(p, 'DisplayFigures', true);
addParameter(p, 'Verbose', true);

parse(p, varargin{:});
opts = p.Results;

opts.DisplayFigures = logical(opts.DisplayFigures);
opts.Verbose = logical(opts.Verbose);

%% ============================ Basic checks ===============================

OCT_fs = 1 / Ts_OCT;

[nX, nY, nT] = size(WavePropSim_3DMat);

if ~isequal(size(Inc_Mask), [nX, nY])
    error(['Inc_Mask must have the same spatial dimensions as the input ' ...
           'simulation matrix, i.e. size(Inc_Mask) = [nX, nY].']);
end

if ~isequal(size(Bkg_Mask), [nX, nY])
    error(['Bkg_Mask must have the same spatial dimensions as the input ' ...
           'simulation matrix, i.e. size(Bkg_Mask) = [nX, nY].']);
end

Inc_Mask = logical(Inc_Mask);
Bkg_Mask = logical(Bkg_Mask);

%% ===================== Analysis support from masks =======================
% Combine background and inclusion masks to recover the external circular
% valid domain. The hole left by the inclusion is filled so the estimator is
% still allowed to run inside the inclusion.

AnalysisSupportMask = Bkg_Mask | Inc_Mask;
AnalysisSupportMask = fill_holes_custom(AnalysisSupportMask);

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

if opts.Verbose
    fprintf('PhaseGradient2D_TV_FromSims()\n');
    fprintf('PG-TV estimator: f_exc=%.2f Hz, omega0=%.2f rad/s, Ts_OCT=%.3g s (fs=%.2f Hz)\n', ...
        f_exc, omega0, Ts_OCT, OCT_fs);
    fprintf('Temporal Hann windowing: %d\n', applyTemporalWindow);
end

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

if opts.Verbose
    fprintf('Nearest FFT bin: kpos=%d -> f_bin=%.2f Hz (|f_bin-f_exc|=%.2f Hz)\n', ...
        kpos, f_bin, abs(f_bin - f_exc));
end

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

        % Do not estimate if the full PG kernel extends outside the valid
        % circular support.
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

Elastogram_raw = Elastogram_c2D(validX, validY);

Xcrop_mm = Xaxis_mm(validX);
Ycrop_mm = Yaxis_mm(validY);

R2_crop = R2_map(validX, validY);

Inc_Mask_Crop = logical(Inc_Mask(validX, validY));
Bkg_Mask_Crop = logical(Bkg_Mask(validX, validY));

AnalysisSupportMask_Crop = logical(AnalysisSupportMask(validX, validY));

Elastogram_raw(~AnalysisSupportMask_Crop) = nan;

Kx_crop   = Kx_map(validX, validY);
Ky_crop   = Ky_map(validX, validY);
Kmag_crop = Kmag_map(validX, validY);

Kx_crop(~AnalysisSupportMask_Crop)   = nan;
Ky_crop(~AnalysisSupportMask_Crop)   = nan;
Kmag_crop(~AnalysisSupportMask_Crop) = nan;

%% ================= TV regularization on cropped K magnitude ==============

TVMask = isfinite(Kmag_crop) & ...
         AnalysisSupportMask_Crop & ...
         Kmag_crop >= opts.MinK & ...
         Kmag_crop <= opts.MaxK & ...
         R2_crop >= opts.R2Threshold;

Kmag_for_TV = Kmag_crop;
Kmag_for_TV(~TVMask) = NaN;

Kmag_TV_crop = TVDenoise_ChambollePock2D( ...
    Kmag_for_TV, opts.MuTV, TVMask, ...
    'NIterations', opts.TVIterations, ...
    'Tau', opts.TVTau, ...
    'Sigma', opts.TVSigma, ...
    'Theta', opts.TVTheta);

Kmag_TV_crop(~TVMask) = NaN;
Kmag_TV_crop(Kmag_TV_crop < opts.MinK) = NaN;
Kmag_TV_crop(Kmag_TV_crop > opts.MaxK) = NaN;

Elastogram_TV = omega0 ./ Kmag_TV_crop;
Elastogram_TV(~AnalysisSupportMask_Crop) = nan;

%% ============================ Masked elastogram ==========================

Elastogram_Masked_TV = Elastogram_TV;

Elastogram_Masked_TV(R2_crop < opts.R2Threshold) = nan;
Elastogram_Masked_TV(~AnalysisSupportMask_Crop) = nan;

%% ============================== Info struct ==============================

TVInfo = struct();

TVInfo.Method = 'PG-TV';

TVInfo.f_exc = f_exc;
TVInfo.f_bin = f_bin;
TVInfo.omega0 = omega0;
TVInfo.OCT_fs = OCT_fs;
TVInfo.Ts_OCT = Ts_OCT;
TVInfo.T_excitation = T_excitation;

TVInfo.Nfft = Nfft;
TVInfo.kpos = kpos;

TVInfo.kernel_size = L;
TVInfo.halfL = halfL;

TVInfo.MuTV = opts.MuTV;
TVInfo.TVIterations = opts.TVIterations;
TVInfo.TVTau = opts.TVTau;
TVInfo.TVSigma = opts.TVSigma;
TVInfo.TVTheta = opts.TVTheta;

TVInfo.R2Threshold = opts.R2Threshold;
TVInfo.MinK = opts.MinK;
TVInfo.MaxK = opts.MaxK;

TVInfo.ValidFractionForTV = nnz(TVMask) / numel(TVMask);

TVInfo.validX = validX;
TVInfo.validY = validY;

TVInfo.AnalysisSupportMask = AnalysisSupportMask;
TVInfo.AnalysisSupportMask_Crop = AnalysisSupportMask_Crop;
TVInfo.TVMASK = TVMask;

if opts.Verbose
    fprintf('TV mu = %.4g\n', opts.MuTV);
    fprintf('TV iterations = %d\n', opts.TVIterations);
    fprintf('R2 threshold = %.2f\n', opts.R2Threshold);
    fprintf('Valid TV pixels = %.2f %%\n', 100*TVInfo.ValidFractionForTV);
end

%% ================================ Figures ================================

if opts.DisplayFigures

    figure;
    colormap("turbo");
    imagesc(Xcrop_mm, Ycrop_mm, Elastogram_raw);
    axis image;
    xlabel('X [mm]');
    ylabel('Y [mm]');
    title({sprintf('Raw 2D Phase-Gradient Speed Elastogram @ f_{exc}=%.1f Hz', f_exc), ...
        sprintf('Kernel size: L=%dx%d px', L, L)});
    colorbar;
    clim(clim_limits);
    grid off;
    hold on;
    plot_mask_outline(Xcrop_mm, Ycrop_mm, Inc_Mask_Crop, 'w', 2.5);
    hold off;

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

    figure;
    colormap("turbo");
    imagesc(Xcrop_mm, Ycrop_mm, Elastogram_Masked_TV);
    axis image;
    xlabel('X [mm]');
    ylabel('Y [mm]');
    title({'2D PG-TV Masked Elastogram', ...
        sprintf('(R^2 >= %.2f), L=%dx%d px, \\mu=%.3g', ...
        opts.R2Threshold, L, L, opts.MuTV)});
    colorbar;
    clim(clim_limits);
    grid off;
    hold on;
    plot_mask_outline(Xcrop_mm, Ycrop_mm, Inc_Mask_Crop, 'w', 2.5);
    hold off;

    % figure('Color','w','Name','PG-TV K Magnitude Diagnostics');
    % tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    % 
    % nexttile;
    % imagesc(Xcrop_mm, Ycrop_mm, Kmag_crop);
    % axis image;
    % xlabel('X [mm]');
    % ylabel('Y [mm]');
    % title('Raw K Magnitude [rad/m]');
    % colormap(gca,"turbo");
    % colorbar;
    % grid off;
    % hold on;
    % plot_mask_outline(Xcrop_mm, Ycrop_mm, Inc_Mask_Crop, 'w', 2.0);
    % hold off;
    % 
    % nexttile;
    % imagesc(Xcrop_mm, Ycrop_mm, Kmag_TV_crop);
    % axis image;
    % xlabel('X [mm]');
    % ylabel('Y [mm]');
    % title(sprintf('TV-Regularized K Magnitude [rad/m], \\mu=%.3g', opts.MuTV));
    % colormap(gca,"turbo");
    % colorbar;
    % grid off;
    % hold on;
    % plot_mask_outline(Xcrop_mm, Ycrop_mm, Inc_Mask_Crop, 'w', 2.0);
    % hold off;

end

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

%% ===================================================================== %%
function FilledMask = fill_holes_custom(Mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill_holes_custom
%
% Fills holes in a 2D binary mask without requiring imfill().
%
% Logic:
%   1) Identify background pixels connected to the image boundary.
%   2) Internal background pixels not connected to the boundary are holes.
%   3) Fill those holes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mask = logical(Mask);

[Nx, Ny] = size(Mask);

Background = ~Mask;
Visited = false(Nx, Ny);

queue = zeros(numel(Mask), 2);
head = 1;
tail = 0;

% Add boundary background pixels to the flood-fill queue.
boundary = false(Nx, Ny);
boundary(1,:) = true;
boundary(end,:) = true;
boundary(:,1) = true;
boundary(:,end) = true;

boundary_bg = find(Background & boundary);

for k = 1:numel(boundary_bg)
    [i,j] = ind2sub([Nx,Ny], boundary_bg(k));
    tail = tail + 1;
    queue(tail,:) = [i,j];
    Visited(i,j) = true;
end

% 4-connected flood fill from boundary background.
while head <= tail

    i = queue(head,1);
    j = queue(head,2);
    head = head + 1;

    neighbors = [i-1 j; i+1 j; i j-1; i j+1];

    for n = 1:4

        ii = neighbors(n,1);
        jj = neighbors(n,2);

        if ii < 1 || ii > Nx || jj < 1 || jj > Ny
            continue;
        end

        if Background(ii,jj) && ~Visited(ii,jj)
            Visited(ii,jj) = true;
            tail = tail + 1;
            queue(tail,:) = [ii,jj];
        end
    end
end

BoundaryConnectedBackground = Visited;
InternalHoles = Background & ~BoundaryConnectedBackground;

FilledMask = Mask | InternalHoles;

end