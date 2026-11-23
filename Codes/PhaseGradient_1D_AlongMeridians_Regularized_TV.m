function [Elastogram_TV,R2_map,Elastogram_Masked_TV,Xgrid_mm,Ygrid_mm,...
    SparseElastogram,Inc_Mask_Grid,KVectorField] = ...
    PhaseGradient_1D_AlongMeridians_Regularized_TV(Meridian_Frames,Inc_Mask,...
    kernel_size,Ts_OCT,T_excitation,...
    horizontalMin,horizontalMax,verticalMin,verticalMax,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:
%   Estéfano Machado Ferrer
%
% Function:
%   PhaseGradient_1D_AlongMeridians_Regularized_TV
%
% Description:
%   This function implements a regularized MB-mode 1D Phase Gradient
%   estimator for wave-based OCE simulations.
%
%   The estimator follows the acquisition-consistent 1D Phase Gradient
%   workflow:
%
%       1) Estimate local phase gradients along each MB meridian.
%       2) Convert each local 1D phase-gradient estimate into a local
%          wave-number magnitude.
%       3) Project the sparse meridian estimates into the reduced Cartesian
%          en-face MB grid.
%       4) Interpolate the sparse wave-number map into a dense Cartesian
%          wave-number map.
%       5) Apply 2D Total Variation regularization to the interpolated
%          wave-number map.
%       6) Convert the TV-regularized wave-number map into phase speed:
%
%              c_TV(x,y) = omega0 / K_TV(x,y)
%
%   The key point is that the TV regularization is applied to the
%   interpolated wave-number map K(x,y), not directly to the phase-speed
%   elastogram c(x,y). This matches the PG-TV inverse problem:
%
%       K_TV = argmin_K 0.5 ||K - K_raw||_2^2 + mu TV(K)
%
%   where:
%
%       K_raw : interpolated wave-number magnitude map from 1D PG
%       K_TV  : TV-regularized wave-number magnitude map
%       mu    : TV regularization parameter
%
%   The final elastogram is then computed from K_TV.
%
% Required external function:
%   TVDenoise_ChambollePock2D.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%   Meridian_Frames:
%       3D matrix containing MB-mode meridian space-time data.
%       Expected size:
%
%           [nAscans x nT x nBscans]
%
%       where:
%           nAscans = spatial samples per meridian
%           nT      = temporal samples
%           nBscans = number of meridians
%
%   Inc_Mask:
%       Binary inclusion mask. If its size does not match the reduced
%       output grid [nAscans x nAscans], it is resized using nearest-neighbor
%       interpolation.
%
%   kernel_size:
%       Local 1D fitting window size. If even, it is automatically increased
%       by one to guarantee an odd window.
%
%   Ts_OCT:
%       OCT sampling period [s].
%
%   T_excitation:
%       Excitation period [s].
%
%   horizontalMin, horizontalMax:
%       X-axis limits of the simulated field [mm].
%
%   verticalMin, verticalMax:
%       Y-axis limits of the simulated field [mm].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional name-value parameters:
%
%   'R2_thresh':
%       Final R² threshold used for the masked elastogram.
%       Default: 0.70
%
%   'DisplayResults':
%       If true, display sparse, raw interpolated, regularized, R², and
%       masked elastograms.
%       Default: true
%
%   'ApplyTemporalHann':
%       If true, apply Hann windowing before temporal FFT.
%       Default: true
%
%   'InterpMethod':
%       Scattered interpolation method.
%       Typical values: 'natural', 'linear', 'nearest'
%       Default: 'natural'
%
%   'MuTV':
%       TV regularization parameter applied to the interpolated K-map.
%       Since K is expressed in rad/m, values such as 50–500 are reasonable.
%       Default: 250
%
%   'TVIterations':
%       Number of Chambolle-Pock iterations.
%       Default: 500
%
%   'TVTau', 'TVSigma', 'TVTheta':
%       Chambolle-Pock numerical parameters.
%       Defaults: 0.25, 0.25, 1.0
%
%   'TV_UseR2Mask':
%       If true, TV is applied only where R² >= TV_R2Threshold.
%       If false, TV is applied over the valid interpolated K-map support.
%       Default: false
%
%   'TV_R2Threshold':
%       R² threshold used only if TV_UseR2Mask is true.
%       Default: 0.70
%
%   'MinK', 'MaxK':
%       Valid wave-number range [rad/m] for TV regularization.
%       Default: 500 and 8000
%
%   'CLim':
%       Color limits for displayed speed maps [m/s].
%       Default: [4 12]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%
%   Elastogram_TV:
%       TV-regularized interpolated phase-speed elastogram [m/s].
%
%   R2_map:
%       Interpolated R² map obtained from 1D PG local linear fits.
%
%   Elastogram_Masked_TV:
%       TV-regularized elastogram masked using R2_thresh.
%
%   Xgrid_mm, Ygrid_mm:
%       Reduced Cartesian output grid axes [mm].
%
%   SparseElastogram:
%       Sparse Cartesian elastogram before interpolation.
%
%   Inc_Mask_Grid:
%       Inclusion mask resized/aligned to the reduced MB grid.
%
%   KVectorField:
%       Struct containing raw and regularized K maps, k-vector maps,
%       meridian estimates, masks, and TV metadata.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options ==============================

P = inputParser;

addParameter(P,'R2_thresh',0.70);
addParameter(P,'DisplayResults',true);
addParameter(P,'ApplyTemporalHann',true);
addParameter(P,'InterpMethod','natural');

addParameter(P,'MuTV',250);
addParameter(P,'TVIterations',500);
addParameter(P,'TVTau',0.25);
addParameter(P,'TVSigma',0.25);
addParameter(P,'TVTheta',1.0);
addParameter(P,'TV_UseR2Mask',false);
addParameter(P,'TV_R2Threshold',0.70);
addParameter(P,'MinK',500);
addParameter(P,'MaxK',8000);
addParameter(P,'CLim',[4 12]);

parse(P,varargin{:});
Opt = P.Results;

%% ===================== Preliminary calculations ==========================

OCT_fs = 1 / Ts_OCT;
f_exc  = 1 / T_excitation;
omega0 = 2*pi*f_exc;

[nAscans,nT,nBscans] = size(Meridian_Frames);

fprintf('\nPhaseGradient_1D_AlongMeridians_Regularized_TV()\n');
fprintf('Input Meridian_Frames size: [%d x %d x %d]\n', nAscans, nT, nBscans);
fprintf('f_exc = %.2f Hz, omega0 = %.2f rad/s, Ts_OCT = %.3g s, fs = %.2f Hz\n', ...
    f_exc, omega0, Ts_OCT, OCT_fs);

L = kernel_size;

if mod(L,2) == 0
    L = L + 1;
    warning('kernel_size was even; using L = %d instead.', L);
end

halfL = floor(L/2);

%% ===================== Reduced output Cartesian grid =====================

Nx_out = nAscans;
Ny_out = nAscans;

Xgrid_mm = linspace(horizontalMin,horizontalMax,Nx_out);
Ygrid_mm = linspace(verticalMin,verticalMax,Ny_out);

Xgrid_m = Xgrid_mm * 1e-3;
Ygrid_m = Ygrid_mm * 1e-3;

fprintf('Output reduced MB elastogram grid: [%d x %d]\n', Nx_out, Ny_out);

x_center_mm = (horizontalMin + horizontalMax) / 2;
y_center_mm = (verticalMin   + verticalMax)   / 2;

Rx_mm = (horizontalMax - horizontalMin) / 2;
Ry_mm = (verticalMax   - verticalMin)   / 2;

%% ===================== Inclusion mask alignment ==========================

if nargin >= 2 && ~isempty(Inc_Mask)

    if ~isequal(size(Inc_Mask),[Nx_out,Ny_out])
        Inc_Mask_Grid = imresize(logical(Inc_Mask),[Nx_out,Ny_out],'nearest');
    else
        Inc_Mask_Grid = logical(Inc_Mask);
    end

else
    Inc_Mask_Grid = false(Nx_out,Ny_out);
end

%% ===================== Rebuild exact MB acquisition geometry =============

r = linspace(-1,1,nAscans).';

theta = linspace(0,pi,nBscans+1);
theta = theta(1:end-1);

rr = repmat(r,1,nBscans);
tt = repmat(theta,nAscans,1);

% Match multimeridian acquisition convention:
% every second meridian is acquired in reverse radial order.
rr(:,2:2:end) = rr(end:-1:1,2:2:end);

X_mer_mm = x_center_mm + (rr .* cos(tt)) * Rx_mm;
Y_mer_mm = y_center_mm + (rr .* sin(tt)) * Ry_mm;

Rref_mm = mean([Rx_mm,Ry_mm]);
s_mm = rr * Rref_mm;
s_m  = s_mm * 1e-3;

fprintf('Reconstructed MB acquisition geometry using reduced MB grid.\n');

%% ===================== Clean temporal traces =============================

Meridian_Clean = nan(size(Meridian_Frames));

for b = 1:nBscans
    for a = 1:nAscans

        trace = squeeze(Meridian_Frames(a,:,b)).';

        if all(isnan(trace))
            continue;
        end

        trace = fill_nans_1d(trace);
        Meridian_Clean(a,:,b) = trace;
    end
end

%% ===================== FFT at excitation frequency =======================

if Opt.ApplyTemporalHann

    hann_win = hanning(size(Meridian_Clean,2)).';

    hann_win_3dmat = repmat(reshape(hann_win,1,[],1), ...
        size(Meridian_Clean,1),1,size(Meridian_Clean,3));

    Meridian_Clean = Meridian_Clean .* hann_win_3dmat;
end

Nfft = 2^nextpow2(max(nT,1024));

FFT_Mer = fft(Meridian_Clean,Nfft,2);

freqs = (0:(Nfft-1)) * (OCT_fs / Nfft);

[~,kpos] = min(abs(freqs - f_exc));
f_bin = freqs(kpos);

fprintf('Nearest FFT bin: kpos = %d -> f_bin = %.2f Hz, error = %.2f Hz\n', ...
    kpos, f_bin, abs(f_bin - f_exc));

Complex_Field = squeeze(FFT_Mer(:,kpos,:));

%% ===================== 1D Phase Gradient per meridian ====================

Speed_Meridian = nan(nAscans,nBscans);
R2_Meridian    = nan(nAscans,nBscans);
K_Meridian     = nan(nAscans,nBscans);

for b = 1:nBscans

    cplx_line = Complex_Field(:,b);

    if all(isnan(cplx_line))
        continue;
    end

    theta_line   = angle(cplx_line);
    theta_line_u = unwrap(theta_line);

    amp_line = abs(cplx_line);
    amp_thresh = 0.05 * max(amp_line(:) + eps);

    s_line = s_m(:,b);

    for i = (1+halfL):(nAscans-halfL)

        idx = (i-halfL):(i+halfL);

        s_win     = s_line(idx);
        theta_win = theta_line_u(idx);
        amp_win   = amp_line(idx);

        if any(~isfinite(theta_win)) || any(~isfinite(s_win))
            continue;
        end

        if mean(amp_win,'omitnan') < amp_thresh
            continue;
        end

        A = [s_win(:), ones(numel(s_win),1)];
        coef = A \ theta_win(:);

        k_est = coef(1);

        theta_hat = A * coef;

        SSres = sum((theta_win(:) - theta_hat).^2);
        SStot = sum((theta_win(:) - mean(theta_win(:))).^2) + eps;

        R2 = 1 - SSres/SStot;

        if abs(k_est) < eps
            continue;
        end

        c_est = omega0 / abs(k_est);

        Speed_Meridian(i,b) = c_est;
        R2_Meridian(i,b)    = R2;
        K_Meridian(i,b)     = k_est;
    end
end

%% ===================== Build 2D k-vector field ===========================

Kx_Meridian   = nan(nAscans,nBscans);
Ky_Meridian   = nan(nAscans,nBscans);
Kmag_Meridian = abs(K_Meridian);

for b = 1:nBscans

    dx_ds = Rx_mm * cos(theta(b));
    dy_ds = Ry_mm * sin(theta(b));

    dir_norm = hypot(dx_ds,dy_ds) + eps;

    ux = dx_ds / dir_norm;
    uy = dy_ds / dir_norm;

    Kx_tmp = K_Meridian(:,b) * ux;
    Ky_tmp = K_Meridian(:,b) * uy;

    rx = X_mer_mm(:,b) - x_center_mm;
    ry = Y_mer_mm(:,b) - y_center_mm;

    inward = (Kx_tmp .* rx + Ky_tmp .* ry) < 0;

    Kx_tmp(inward) = -Kx_tmp(inward);
    Ky_tmp(inward) = -Ky_tmp(inward);

    Kx_Meridian(:,b) = Kx_tmp;
    Ky_Meridian(:,b) = Ky_tmp;
end

%% ===================== Project estimates to en face ======================

X_scatter_mm = X_mer_mm(:);
Y_scatter_mm = Y_mer_mm(:);

C_scatter    = Speed_Meridian(:);
R_scatter    = R2_Meridian(:);

Kx_scatter   = Kx_Meridian(:);
Ky_scatter   = Ky_Meridian(:);
Kmag_scatter = Kmag_Meridian(:);

valid_scatter = isfinite(X_scatter_mm) & isfinite(Y_scatter_mm) & ...
                isfinite(C_scatter)    & isfinite(R_scatter)    & ...
                isfinite(Kx_scatter)   & isfinite(Ky_scatter)   & ...
                isfinite(Kmag_scatter) & ...
                (C_scatter > 0);

Xv_mm = X_scatter_mm(valid_scatter);
Yv_mm = Y_scatter_mm(valid_scatter);

Cv    = C_scatter(valid_scatter);
Rv    = R_scatter(valid_scatter);

Kxv   = Kx_scatter(valid_scatter);
Kyv   = Ky_scatter(valid_scatter);
Kmagv = Kmag_scatter(valid_scatter);

fprintf('Valid meridian estimates projected to en face: %d\n',numel(Cv));

if isempty(Cv)
    error('No valid 1D phase-gradient estimates were obtained.');
end

Xv_m = Xv_mm * 1e-3;
Yv_m = Yv_mm * 1e-3;

%% ===================== Sparse Cartesian elastogram =======================

SparseElastogram = nan(Nx_out,Ny_out);
SparseR2         = nan(Nx_out,Ny_out);
SparseKmag       = nan(Nx_out,Ny_out);

SumC   = zeros(Nx_out,Ny_out);
CountC = zeros(Nx_out,Ny_out);

SumR   = zeros(Nx_out,Ny_out);
CountR = zeros(Nx_out,Ny_out);

SumK   = zeros(Nx_out,Ny_out);
CountK = zeros(Nx_out,Ny_out);

for k = 1:numel(Cv)

    [~,ix] = min(abs(Xgrid_mm - Xv_mm(k)));
    [~,iy] = min(abs(Ygrid_mm - Yv_mm(k)));

    SumC(ix,iy)   = SumC(ix,iy) + Cv(k);
    CountC(ix,iy) = CountC(ix,iy) + 1;

    SumR(ix,iy)   = SumR(ix,iy) + Rv(k);
    CountR(ix,iy) = CountR(ix,iy) + 1;

    SumK(ix,iy)   = SumK(ix,iy) + Kmagv(k);
    CountK(ix,iy) = CountK(ix,iy) + 1;
end

valid_sparse_c = CountC > 0;
valid_sparse_r = CountR > 0;
valid_sparse_k = CountK > 0;

SparseElastogram(valid_sparse_c) = SumC(valid_sparse_c) ./ CountC(valid_sparse_c);
SparseR2(valid_sparse_r)         = SumR(valid_sparse_r) ./ CountR(valid_sparse_r);
SparseKmag(valid_sparse_k)       = SumK(valid_sparse_k) ./ CountK(valid_sparse_k);

%% ===================== Interpolate AFTER 1D estimation ===================

[Xq_m,Yq_m] = ndgrid(Xgrid_m,Ygrid_m);

F_c    = scatteredInterpolant(Xv_m,Yv_m,Cv,    Opt.InterpMethod,'none');
F_r2   = scatteredInterpolant(Xv_m,Yv_m,Rv,    Opt.InterpMethod,'none');
F_kx   = scatteredInterpolant(Xv_m,Yv_m,Kxv,   Opt.InterpMethod,'none');
F_ky   = scatteredInterpolant(Xv_m,Yv_m,Kyv,   Opt.InterpMethod,'none');
F_kmag = scatteredInterpolant(Xv_m,Yv_m,Kmagv, Opt.InterpMethod,'none');

Elastogram_raw = F_c(Xq_m,Yq_m);
R2_map         = F_r2(Xq_m,Yq_m);

Kx_map      = F_kx(Xq_m,Yq_m);
Ky_map      = F_ky(Xq_m,Yq_m);
Kmag_rawmap = F_kmag(Xq_m,Yq_m);

[Xq_mm,Yq_mm] = ndgrid(Xgrid_mm,Ygrid_mm);

rx_map = Xq_mm - x_center_mm;
ry_map = Yq_mm - y_center_mm;

inward_map = (Kx_map .* rx_map + Ky_map .* ry_map) < 0;

Kx_map(inward_map) = -Kx_map(inward_map);
Ky_map(inward_map) = -Ky_map(inward_map);

% K-vector magnitude from interpolated Kx/Ky is kept as a diagnostic.
Kmag_vecmap = hypot(Kx_map,Ky_map);

%% ===================== Elliptical support mask ===========================

Rnorm = ((Xq_m - x_center_mm*1e-3)/(Rx_mm*1e-3)).^2 + ...
        ((Yq_m - y_center_mm*1e-3)/(Ry_mm*1e-3)).^2;

inside_mask = Rnorm <= 1;

Elastogram_raw(~inside_mask) = nan;
R2_map(~inside_mask)         = nan;

Kx_map(~inside_mask)         = nan;
Ky_map(~inside_mask)         = nan;
Kmag_rawmap(~inside_mask)    = nan;
Kmag_vecmap(~inside_mask)    = nan;

%% ===================== 2D TV on interpolated K-map =======================

TVMask = inside_mask & isfinite(Kmag_rawmap) & ...
         Kmag_rawmap >= Opt.MinK & ...
         Kmag_rawmap <= Opt.MaxK;

if Opt.TV_UseR2Mask
    TVMask = TVMask & isfinite(R2_map) & (R2_map >= Opt.TV_R2Threshold);
end

fprintf('Applying 2D TV to interpolated 1D-PG K-map: mu = %.3g, iterations = %d\n', ...
    Opt.MuTV, Opt.TVIterations);
fprintf('Valid TV pixels = %.2f %%\n',100*nnz(TVMask)/numel(TVMask));

Kmag_TV_map = TVDenoise_ChambollePock2D( ...
    Kmag_rawmap, Opt.MuTV, TVMask, ...
    'NIterations', Opt.TVIterations, ...
    'Tau', Opt.TVTau, ...
    'Sigma', Opt.TVSigma, ...
    'Theta', Opt.TVTheta);

Kmag_TV_map(~inside_mask) = nan;
Kmag_TV_map(Kmag_TV_map < Opt.MinK) = nan;
Kmag_TV_map(Kmag_TV_map > Opt.MaxK) = nan;

%% ===================== Convert regularized K-map to speed ================

Elastogram_TV = omega0 ./ Kmag_TV_map;
Elastogram_TV(~isfinite(Kmag_TV_map)) = nan;
Elastogram_TV(Kmag_TV_map <= eps) = nan;

%% ===================== R2 masking after TV ===============================

R2_thresh = Opt.R2_thresh;

Elastogram_Masked_TV = Elastogram_TV;
Elastogram_Masked_TV(R2_map < R2_thresh) = nan;

%% ===================== Store k-vector field and TV info ==================

KVectorField = struct();

KVectorField.Xgrid_mm = Xgrid_mm;
KVectorField.Ygrid_mm = Ygrid_mm;

KVectorField.Kx_map = Kx_map;
KVectorField.Ky_map = Ky_map;

KVectorField.Kmag_rawmap = Kmag_rawmap;
KVectorField.Kmag_vecmap = Kmag_vecmap;
KVectorField.Kmag_TV_map = Kmag_TV_map;

KVectorField.Kx_Meridian = Kx_Meridian;
KVectorField.Ky_Meridian = Ky_Meridian;
KVectorField.Kmag_Meridian = Kmag_Meridian;
KVectorField.K_Meridian = K_Meridian;
KVectorField.Speed_Meridian = Speed_Meridian;

KVectorField.X_mer_mm = X_mer_mm;
KVectorField.Y_mer_mm = Y_mer_mm;

KVectorField.R2_map = R2_map;
KVectorField.R2_Meridian = R2_Meridian;

KVectorField.SparseKmag = SparseKmag;
KVectorField.SparseR2 = SparseR2;

KVectorField.Elastogram_raw = Elastogram_raw;
KVectorField.Elastogram_TV = Elastogram_TV;
KVectorField.Elastogram_Masked_TV = Elastogram_Masked_TV;

KVectorField.TVMask = TVMask;
KVectorField.inside_mask = inside_mask;
KVectorField.R2_thresh = R2_thresh;

KVectorField.TVInfo = struct();
KVectorField.TVInfo.Method = '1D-PG followed by interpolation of K-map and 2D-TV on K';
KVectorField.TVInfo.MuTV = Opt.MuTV;
KVectorField.TVInfo.TVIterations = Opt.TVIterations;
KVectorField.TVInfo.TVTau = Opt.TVTau;
KVectorField.TVInfo.TVSigma = Opt.TVSigma;
KVectorField.TVInfo.TVTheta = Opt.TVTheta;
KVectorField.TVInfo.TV_UseR2Mask = Opt.TV_UseR2Mask;
KVectorField.TVInfo.TV_R2Threshold = Opt.TV_R2Threshold;
KVectorField.TVInfo.ValidTVFraction = nnz(TVMask)/numel(TVMask);
KVectorField.TVInfo.MinK = Opt.MinK;
KVectorField.TVInfo.MaxK = Opt.MaxK;
KVectorField.TVInfo.f_exc = f_exc;
KVectorField.TVInfo.omega0 = omega0;
KVectorField.TVInfo.Nfft = Nfft;
KVectorField.TVInfo.f_bin = f_bin;

%% ===================== Display results ==================================

if Opt.DisplayResults

    figure;
    colormap("turbo");
    imagesc(Ygrid_mm,Xgrid_mm,SparseElastogram);
    axis image;
    xlabel('Y [mm]');
    ylabel('X [mm]');
    title('Sparse En Face Elastogram Before Interpolation');
    colorbar;
    clim(Opt.CLim);
    hold on;
    plot_mask_outline(Ygrid_mm,Xgrid_mm,Inc_Mask_Grid,'w',2.5);
    hold off;

    figure;
    colormap("turbo");
    imagesc(Ygrid_mm,Xgrid_mm,Elastogram_raw);
    axis image;
    xlabel('Y [mm]');
    ylabel('X [mm]');
    title('Interpolated MB Elastogram from 1D Phase Gradient');
    colorbar;
    clim(Opt.CLim);
    % grid on;
    hold on;
    plot_mask_outline(Ygrid_mm,Xgrid_mm,Inc_Mask_Grid,'w',2.5);
    hold off;
    % 
    % figure;
    % colormap("turbo");
    % imagesc(Ygrid_mm,Xgrid_mm,Kmag_rawmap);
    % axis image;
    % xlabel('Y [mm]');
    % ylabel('X [mm]');
    % title('Interpolated K Magnitude Map from 1D Phase Gradient [rad/m]');
    % colorbar;
    % grid on;
    % hold on;
    % plot_mask_outline(Ygrid_mm,Xgrid_mm,Inc_Mask_Grid,'w',2.5);
    % hold off;
    % 
    % figure;
    % colormap("turbo");
    % imagesc(Ygrid_mm,Xgrid_mm,Kmag_TV_map);
    % axis image;
    % xlabel('Y [mm]');
    % ylabel('X [mm]');
    % title(sprintf('TV-Regularized K Magnitude Map from 1D PG, \\mu = %.3g',Opt.MuTV));
    % colorbar;
    % grid on;
    % hold on;
    % plot_mask_outline(Ygrid_mm,Xgrid_mm,Inc_Mask_Grid,'w',2.5);
    % hold off;

    % figure;
    % colormap("turbo");
    % imagesc(Ygrid_mm,Xgrid_mm,Elastogram_TV);
    % axis image;
    % xlabel('Y [mm]');
    % ylabel('X [mm]');
    % title(sprintf('Regularized Interpolated MB Elastogram from 1D PG-TV, \\mu = %.3g',Opt.MuTV));
    % colorbar;
    % clim(Opt.CLim);
    % grid on;
    % hold on;
    % plot_mask_outline(Ygrid_mm,Xgrid_mm,Inc_Mask_Grid,'w',2.5);
    % hold off;
    % 
    % figure;
    % colormap("turbo");
    % imagesc(Ygrid_mm,Xgrid_mm,R2_map);
    % axis image;
    % xlabel('Y [mm]');
    % ylabel('X [mm]');
    % title('Interpolated R^2 Map from 1D Phase Gradient');
    % colorbar;
    % clim([0 1]);
    % grid on;
    % hold on;
    % plot_mask_outline(Ygrid_mm,Xgrid_mm,Inc_Mask_Grid,'w',2.5);
    % hold off;

    figure;
    colormap("turbo");
    imagesc(Ygrid_mm,Xgrid_mm,Elastogram_Masked_TV);
    axis image;
    xlabel('Y [mm]');
    ylabel('X [mm]');
    title(sprintf('Masked Regularized MB Elastogram from 1D PG-TV, R^2 >= %.2f',R2_thresh));
    colorbar;
    clim(Opt.CLim);
    grid off;
    hold on;
    plot_mask_outline(Ygrid_mm,Xgrid_mm,Inc_Mask_Grid,'w',2.5);
    hold off;
end

end

%% ===================================================================== %%
function y = fill_nans_1d(x)
% Fill NaN values in a 1D trace using linear interpolation.
% If only one finite value exists, the trace is filled with that value.

x = x(:);
idx = (1:numel(x)).';

valid = isfinite(x);

if ~any(valid)
    y = x;
    return;
end

if sum(valid) == 1
    y = repmat(x(valid),size(x));
    return;
end

y = x;
y(~valid) = interp1(idx(valid),x(valid),idx(~valid),'linear','extrap');

end

%% ===================================================================== %%
function plot_mask_outline(X_axis,Y_axis,Mask,lineColor,lineWidth)
% Plot binary mask outline over an existing axes.

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

contour(X_axis,Y_axis,double(Mask),[0.5 0.5], ...
    'LineColor',lineColor,'LineWidth',lineWidth);

end