function [Elastogram,R2_map,Elastogram_Masked,Xgrid_mm,Ygrid_mm,SparseElastogram,Inc_Mask_Grid,KVectorField] = ...
    PhaseGradient_1D_AlongMeridians_MBMode_v2_1(Meridian_Frames,Inc_Mask,...
    kernel_size,Ts_OCT,T_excitation,...
    horizontalMin,horizontalMax,verticalMin,verticalMax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 14-04-2026
% Date of last modification (DD/MM/YY): 28-04-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Function's Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function 'PhaseGradient_1D_AlongMeridians_MBMode_v2_1' implements the
% MB-mode estimation framework used with real data:
%
%   1) Estimate phase speed independently along each meridian using a
%      1D phase-gradient estimator on Meridian_Frames.
%   2) Project the estimated values into an en face spatial frame using the
%      EXACT geometry defined in AcqVectors_MultimeridianMB().
%   3) Interpolate only AFTER estimation.
%   4) Plot the sparse, interpolated, and masked elastograms with the
%      outline of the inclusion mask.
%   5) Plot the 2D k-vector field induced by the 1D meridian estimates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ===================== Preliminary calculations =========================

OCT_fs = 1 / Ts_OCT;
f_exc  = 1 / T_excitation;
omega0 = 2*pi*f_exc;

[nAscans, nT, nBscans] = size(Meridian_Frames);

fprintf('PhaseGradient_1D_AlongMeridians_MBMode_v2_1()\n');
fprintf('Input Meridian_Frames size: [%d x %d x %d]\n', nAscans, nT, nBscans);
fprintf('f_exc = %.2f Hz, omega0 = %.2f rad/s, Ts_OCT = %.3g s (fs = %.2f Hz)\n', ...
    f_exc, omega0, Ts_OCT, OCT_fs);

L = kernel_size;
if mod(L,2) == 0
    L = L + 1;
    warning('kernel_size was even; using L = %d instead.', L);
end
halfL = floor(L/2);

%% ===================== Output Cartesian grid ============================

Nx_out = nAscans;
Ny_out = nAscans;

Xgrid_mm = linspace(horizontalMin, horizontalMax, Nx_out);
Ygrid_mm = linspace(verticalMin,   verticalMax,   Ny_out);

Xgrid_m = Xgrid_mm * 1e-3;
Ygrid_m = Ygrid_mm * 1e-3;

x_center_mm = (horizontalMin + horizontalMax) / 2;
y_center_mm = (verticalMin   + verticalMax)   / 2;

Rx_mm = (horizontalMax - horizontalMin) / 2;
Ry_mm = (verticalMax   - verticalMin)   / 2;

fprintf('Output grid: [%d x %d]\n', Nx_out, Ny_out);

%% ===================== Inclusion mask alignment =========================

if nargin >= 2 && ~isempty(Inc_Mask)
    if ~isequal(size(Inc_Mask), [Nx_out, Ny_out])
        Inc_Mask_Grid = imresize(logical(Inc_Mask), [Nx_out, Ny_out], 'nearest');
    else
        Inc_Mask_Grid = logical(Inc_Mask);
    end
else
    Inc_Mask_Grid = false(Nx_out, Ny_out);
end

%% ===================== Rebuild MB acquisition geometry ==================

r = linspace(-1, 1, nAscans).';
theta = linspace(0, pi, nBscans+1);
theta = theta(1:end-1);

rr = repmat(r, 1, nBscans);
tt = repmat(theta, nAscans, 1);

rr(:,2:2:end) = rr(end:-1:1, 2:2:end);

X_mer_mm = x_center_mm + (rr .* cos(tt)) * Rx_mm;
Y_mer_mm = y_center_mm + (rr .* sin(tt)) * Ry_mm;

Rref_mm = mean([Rx_mm, Ry_mm]);
s_mm = rr * Rref_mm;
s_m  = s_mm * 1e-3;

fprintf('Reconstructed exact MB acquisition geometry from vector logic.\n');

%% ===================== Clean traces and FFT =============================

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

% Temporal Hann windowing along the MB dwell-time dimension.
hann_win = hanning(size(Meridian_Clean,2)).';
hann_win_3dmat = repmat(reshape(hann_win,1,[],1), ...
    size(Meridian_Clean,1), 1, size(Meridian_Clean,3));

Meridian_Clean_win = Meridian_Clean .* hann_win_3dmat;

Nfft = 2^nextpow2(max(nT, 1024));
FFT_Mer = fft(Meridian_Clean_win, Nfft, 2);

freqs = (0:(Nfft-1)) * (OCT_fs / Nfft);
[~, kpos] = min(abs(freqs - f_exc));
f_bin = freqs(kpos);

fprintf('Nearest FFT bin: kpos = %d -> f_bin = %.2f Hz (|f_bin-f_exc| = %.2f Hz)\n', ...
    kpos, f_bin, abs(f_bin - f_exc));

Complex_Field = squeeze(FFT_Mer(:,kpos,:));
Amp_Field     = abs(Complex_Field);

%% ===================== 1D Phase Gradient per meridian ===================

Speed_Meridian = nan(nAscans, nBscans);
R2_Meridian    = nan(nAscans, nBscans);
K_Meridian     = nan(nAscans, nBscans);

for b = 1:nBscans

    cplx_line = Complex_Field(:,b);

    if all(isnan(cplx_line))
        continue;
    end

    theta_line = angle(cplx_line);
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

        if mean(amp_win, 'omitnan') < amp_thresh
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

%% ===================== Build 2D k-vector field ==========================
% The 1D meridian-estimated scalar k values are projected into Cartesian
% components and then interpolated into an en face kx, ky map. The vectors
% plotted later are taken from these interpolated kx, ky maps, not directly
% from the meridian lines.

Kx_Meridian = nan(nAscans, nBscans);
Ky_Meridian = nan(nAscans, nBscans);
Kmag_Meridian = abs(K_Meridian);

for b = 1:nBscans

    dx_ds = Rx_mm * cos(theta(b));
    dy_ds = Ry_mm * sin(theta(b));

    dir_norm = hypot(dx_ds, dy_ds) + eps;

    ux = dx_ds / dir_norm;
    uy = dy_ds / dir_norm;

    Kx_tmp = K_Meridian(:,b) * ux;
    Ky_tmp = K_Meridian(:,b) * uy;

    % Force outward direction from source/center
    rx = X_mer_mm(:,b) - x_center_mm;
    ry = Y_mer_mm(:,b) - y_center_mm;

    inward = (Kx_tmp .* rx + Ky_tmp .* ry) < 0;

    Kx_tmp(inward) = -Kx_tmp(inward);
    Ky_tmp(inward) = -Ky_tmp(inward);

    Kx_Meridian(:,b) = Kx_tmp;
    Ky_Meridian(:,b) = Ky_tmp;
end

%% ===================== Project estimates to en face =====================

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

fprintf('Valid meridian estimates projected to en face: %d\n', numel(Cv));

if isempty(Cv)
    error('No valid 1D phase-gradient estimates were obtained.');
end

Xv_m = Xv_mm * 1e-3;
Yv_m = Yv_mm * 1e-3;

%% ===================== Sparse Cartesian elastogram ======================

SparseElastogram = nan(Nx_out, Ny_out);
SparseR2         = nan(Nx_out, Ny_out);

SumC   = zeros(Nx_out, Ny_out);
CountC = zeros(Nx_out, Ny_out);

SumR   = zeros(Nx_out, Ny_out);
CountR = zeros(Nx_out, Ny_out);

for k = 1:numel(Cv)

    [~, ix] = min(abs(Xgrid_mm - Xv_mm(k)));
    [~, iy] = min(abs(Ygrid_mm - Yv_mm(k)));

    SumC(ix,iy)   = SumC(ix,iy) + Cv(k);
    CountC(ix,iy) = CountC(ix,iy) + 1;

    SumR(ix,iy)   = SumR(ix,iy) + Rv(k);
    CountR(ix,iy) = CountR(ix,iy) + 1;
end

valid_sparse_c = CountC > 0;
valid_sparse_r = CountR > 0;

SparseElastogram(valid_sparse_c) = SumC(valid_sparse_c) ./ CountC(valid_sparse_c);
SparseR2(valid_sparse_r)         = SumR(valid_sparse_r) ./ CountR(valid_sparse_r);

%% ===================== Interpolate AFTER estimation =====================

[Xq_m, Yq_m] = ndgrid(Xgrid_m, Ygrid_m);

F_c  = scatteredInterpolant(Xv_m, Yv_m, Cv,  'natural', 'none');
F_r2 = scatteredInterpolant(Xv_m, Yv_m, Rv,  'natural', 'none');
F_kx = scatteredInterpolant(Xv_m, Yv_m, Kxv, 'natural', 'none');
F_ky = scatteredInterpolant(Xv_m, Yv_m, Kyv, 'natural', 'none');

Elastogram = F_c(Xq_m, Yq_m);
R2_map     = F_r2(Xq_m, Yq_m);

Kx_map = F_kx(Xq_m, Yq_m);
Ky_map = F_ky(Xq_m, Yq_m);

% Force interpolated vectors outward from the divergent source
[Xq_mm, Yq_mm] = ndgrid(Xgrid_mm, Ygrid_mm);

rx_map = Xq_mm - x_center_mm;
ry_map = Yq_mm - y_center_mm;

inward_map = (Kx_map .* rx_map + Ky_map .* ry_map) < 0;

Kx_map(inward_map) = -Kx_map(inward_map);
Ky_map(inward_map) = -Ky_map(inward_map);

Kmag_map = hypot(Kx_map, Ky_map);

%% ===================== Circular support mask ============================

Rnorm = ((Xq_m - x_center_mm*1e-3)/(Rx_mm*1e-3)).^2 + ...
        ((Yq_m - y_center_mm*1e-3)/(Ry_mm*1e-3)).^2;

inside_mask = Rnorm <= 1;

Elastogram(~inside_mask) = nan;
R2_map(~inside_mask)     = nan;

Kx_map(~inside_mask)     = nan;
Ky_map(~inside_mask)     = nan;
Kmag_map(~inside_mask)   = nan;

%% ===================== R2 masking =======================================

R2_thresh = 0.70;

Elastogram_Masked = Elastogram;
Elastogram_Masked(R2_map < R2_thresh) = nan;

%% ===================== Store k-vector field output ======================

KVectorField = struct();
KVectorField.Xgrid_mm = Xgrid_mm;
KVectorField.Ygrid_mm = Ygrid_mm;
KVectorField.Kx_map = Kx_map;
KVectorField.Ky_map = Ky_map;
KVectorField.Kmag_map = Kmag_map;
KVectorField.Kx_Meridian = Kx_Meridian;
KVectorField.Ky_Meridian = Ky_Meridian;
KVectorField.Kmag_Meridian = Kmag_Meridian;
KVectorField.X_mer_mm = X_mer_mm;
KVectorField.Y_mer_mm = Y_mer_mm;
KVectorField.R2_map = R2_map;
KVectorField.R2_Meridian = R2_Meridian;

%% ===================== Display results ==================================

figure;
colormap("turbo");
imagesc(Ygrid_mm, Xgrid_mm, SparseElastogram);
axis image;
xlabel('Y [mm]');
ylabel('X [mm]');
title('Sparse En Face Elastogram (Before Interpolation)');
colorbar;
clim([4 12]);
hold on;
plot_mask_outline(Ygrid_mm, Xgrid_mm, Inc_Mask_Grid, 'w', 2.5);
hold off;

figure;
colormap("turbo");
imagesc(Ygrid_mm, Xgrid_mm, Elastogram);
axis image;
xlabel('Y [mm]');
ylabel('X [mm]');
title(sprintf('Interpolated MB Elastogram from 1D Phase Gradient (L=%d)', L));
colorbar;
clim([4 12]);
grid on;
hold on;
plot_mask_outline(Ygrid_mm, Xgrid_mm, Inc_Mask_Grid, 'w', 2.5);
hold off;

figure;
colormap("turbo");
imagesc(Ygrid_mm, Xgrid_mm, R2_map);
axis image;
xlabel('Y [mm]');
ylabel('X [mm]');
title(sprintf('Interpolated R^2 Map from 1D Phase Gradient (L=%d)', L));
colorbar;
clim([0 1]);
grid on;
hold on;
plot_mask_outline(Ygrid_mm, Xgrid_mm, Inc_Mask_Grid, 'w', 2.5);
hold off;

figure;
colormap("turbo");
imagesc(Ygrid_mm, Xgrid_mm, Elastogram_Masked);
axis image;
xlabel('Y [mm]');
ylabel('X [mm]');
title(sprintf('Masked MB Elastogram from 1D Phase Gradient (R^2 >= %.2f)', R2_thresh));
colorbar;
clim([4 12]);
grid off;
hold on;
plot_mask_outline(Ygrid_mm, Xgrid_mm, Inc_Mask_Grid, 'w', 2.5);
hold off;

%% ===================== Display k-vector field ===========================

plot_k_vector_field_MB( ...
    X_mer_mm, Y_mer_mm, ...
    Xgrid_mm, Ygrid_mm, ...
    Kx_map, Ky_map, Kmag_map, ...
    Inc_Mask_Grid, inside_mask, R2_map, R2_thresh);

end

%% ===================================================================== %%
function y = fill_nans_1d(x)
% Fill NaNs in a 1D vector using linear interpolation. If NaNs appear at
% the ends, nearest extrapolation is used.

x = x(:);
idx = (1:numel(x)).';

valid = isfinite(x);

if ~any(valid)
    y = x;
    return;
end

if sum(valid) == 1
    y = repmat(x(valid), size(x));
    return;
end

y = x;
y(~valid) = interp1(idx(valid), x(valid), idx(~valid), 'linear', 'extrap');
end

%% ===================================================================== %%
function plot_k_vector_field_MB(X_mer_mm, Y_mer_mm, ...
    Xgrid_mm, Ygrid_mm, Kx_map, Ky_map, Kmag_map, ...
    Inc_Mask_Grid, inside_mask, R2_map, R2_thresh)

valid = isfinite(Kx_map) & isfinite(Ky_map) & isfinite(Kmag_map) & ...
        inside_mask & isfinite(R2_map) & (R2_map >= R2_thresh);

if ~any(valid(:))
    warning('No valid k-vectors available for plotting.');
    return;
end

[Xg_mm, Yg_mm] = ndgrid(Xgrid_mm, Ygrid_mm);

% Subsample the interpolated k-map to avoid crowding
stride = max(1, round(numel(Xgrid_mm)/24));

sample_mask = false(size(valid));
sample_mask(1:stride:end, 1:stride:end) = true;

plot_mask = valid & sample_mask;

if nnz(plot_mask) < 5
    plot_mask = valid;
end

Xv = Xg_mm(plot_mask);
Yv = Yg_mm(plot_mask);
Kxv = Kx_map(plot_mask);
Kyv = Ky_map(plot_mask);
Kmagv = Kmag_map(plot_mask);

Knorm = hypot(Kxv, Kyv) + eps;
Ux = Kxv ./ Knorm;
Uy = Kyv ./ Knorm;

figure('Name','MB 1D Phase-Gradient k-Vector Field', ...
    'Units','normalized','Position',[0.08 0.10 0.70 0.72]);

hold on;

% Plot multimeridian acquisition pattern in the background
for b = 1:size(X_mer_mm,2)
    plot(X_mer_mm(:,b), Y_mer_mm(:,b), '-', ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    plot(X_mer_mm(:,b), Y_mer_mm(:,b), '.', ...
        'Color', [0.65 0.65 0.65], 'MarkerSize', 5);
end

% Inclusion outline
contour(Xgrid_mm, Ygrid_mm, double(Inc_Mask_Grid), [0.5 0.5], ...
    'LineColor', 'k', 'LineWidth', 2.0);

% Color-coded normalized vectors from interpolated Kx_map, Ky_map
cmap = turbo(256);
kmin = min(Kmagv);
kmax = max(Kmagv);

if abs(kmax - kmin) < eps
    color_idx = ones(size(Kmagv));
else
    color_idx = round(1 + 255*(Kmagv - kmin)/(kmax - kmin));
end
color_idx = min(max(color_idx,1),256);

vec_scale_mm = 0.35;

for i = 1:numel(Kmagv)
    col = cmap(color_idx(i),:);

    quiver(Xv(i), Yv(i), ...
        vec_scale_mm*Ux(i), vec_scale_mm*Uy(i), ...
        0, 'Color', col, 'LineWidth', 1.5, 'MaxHeadSize', 1.8);
end

axis image;
set(gca, 'YDir', 'reverse');
grid on;
xlabel('X [mm]');
ylabel('Y [mm]');
title({'2D k-Vector Field from Interpolated MB 1D Phase Gradient', ...
       'Vectors from K_x and K_y maps; color indicates |k| [rad/m]'});

colormap(gca, turbo);
cb = colorbar;
clim([kmin kmax]);
ylabel(cb, '|k| [rad/m]');

hold off;

end

%% ===================================================================== %%
function plot_mask_outline(X_axis, Y_axis, Mask, lineColor, lineWidth)
% Plot the outline of a binary mask over an existing axes.
%
% This function is intentionally consistent with the plotting convention:
%
%   imagesc(Ygrid_mm, Xgrid_mm, Elastogram)
%
% where the matrix is stored as [X x Y]. Therefore, the mask is also passed
% without transposition.

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

contour(X_axis, Y_axis, double(Mask), [0.5 0.5], ...
    'LineColor', lineColor, 'LineWidth', lineWidth);
end