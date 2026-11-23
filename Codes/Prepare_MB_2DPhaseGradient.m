function [WavePropSim_3DMat,ValidMask_MB,HullMask_MB,InterpOperator_MB] = ...
    Prepare_MB_2DPhaseGradient(SampledFrames_MB,interpMethod,fillOutsideHull,...
                                applySmoothing,sigma_px)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Description:
% Prepare the sparse MB-mode sampled movie for 2D Phase Gradient
% estimation by performing a physically consistent spatial interpolation
% over time.
% Date of creation (DD/MM/YY): 24-03-2026
% Last modified (DD/MM/YY): 24-03-2026
%
% Inputs:
%   - SampledFrames_MB (3D matrix, double):
%       Sparse MB sampled movie of size [Nx x Ny x Nt], typically with NaNs
%       in non-sampled pixels.
%
%   - interpMethod (char):
%       Interpolation method:
%           'linearDT' -> fixed Delaunay triangulation + barycentric weights
%           'nearest'  -> nearest-neighbor interpolation
%
%   - fillOutsideHull (logical):
%       false -> keep pixels outside convex hull as NaN
%       true  -> fill outside hull using nearest-neighbor extrapolation
%
%   - applySmoothing (logical):
%       If true, apply light Gaussian smoothing after interpolation.
%
%   - sigma_px (scalar, double):
%       Standard deviation of Gaussian smoothing in pixels.
%
% Outputs:
%   - WavePropSim_3DMat (3D matrix, double):
%       Dense/interpolated movie ready to be used directly as input to
%       PhaseGradient2D_FromSims().
%
%   - ValidMask_MB (2D logical matrix):
%       Mask of originally sampled pixels.
%
%   - HullMask_MB (2D logical matrix):
%       Mask of pixels inside the interpolation convex hull.
%
%   - InterpOperator_MB (struct):
%       Fixed interpolation operator reused across all temporal frames.
%
% Notes:
% - The interpolation operator is built once and reused for all time
%   frames, which is physically consistent because the MB sampling geometry
%   is fixed over time.
% - This is preferable to re-triangulating each frame independently.
% - If applySmoothing = true, smoothing is performed in a NaN-aware way.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    SampledFrames_MB (:,:,:) double
    interpMethod (1,:) char {mustBeMember(interpMethod,{'linearDT','nearest'})} = 'linearDT'
    fillOutsideHull (1,1) logical = false
    applySmoothing (1,1) logical = false
    sigma_px (1,1) double {mustBePositive} = 0.8
end

%% ============================ Basic sizes ===============================
[Nx,Ny,Nt] = size(SampledFrames_MB);

fprintf('Preparing MB movie for 2D Phase Gradient...\n');
fprintf('Input sparse movie size: [%d x %d x %d]\n', Nx, Ny, Nt);

%% ======================= Build support mask =============================
% Since the MB acquisition geometry is fixed in time, we define the valid
% support from all temporal frames.

ValidMask_MB = any(~isnan(SampledFrames_MB), 3);
nValid = nnz(ValidMask_MB);

fprintf('Number of originally sampled pixels: %d\n', nValid);

if nValid < 3
    error('Not enough valid sampled pixels to perform 2D interpolation.');
end

%% ================= Build fixed interpolation operator ===================
InterpOperator_MB = localBuildInterpOperator(ValidMask_MB, interpMethod, fillOutsideHull);
HullMask_MB = InterpOperator_MB.HullMask;

fprintf('Number of pixels inside interpolation hull: %d\n', nnz(HullMask_MB));

%% ==================== Interpolate all temporal frames ===================
WavePropSim_3DMat = nan(Nx,Ny,Nt);

for k = 1:Nt
    frame_sparse = SampledFrames_MB(:,:,k);
    frame_interp = localApplyInterpOperator(frame_sparse, InterpOperator_MB);

    % Keep outside-hull region as NaN unless user explicitly requested fill
    if ~fillOutsideHull
        frame_interp(~HullMask_MB) = nan;
    end

    WavePropSim_3DMat(:,:,k) = frame_interp;
end

%% ================= Optional NaN-aware spatial smoothing =================
if applySmoothing
    fprintf('Applying NaN-aware Gaussian smoothing (sigma = %.2f px)...\n', sigma_px);

    for k = 1:Nt
        frame_k = WavePropSim_3DMat(:,:,k);
        mask_k = ~isnan(frame_k);

        if nnz(mask_k) < 10
            continue;
        end

        % NaN-aware smoothing by normalized convolution
        frame_zero = frame_k;
        frame_zero(~mask_k) = 0;

        weight = double(mask_k);

        frame_blur = imgaussfilt(frame_zero, sigma_px);
        weight_blur = imgaussfilt(weight, sigma_px);

        frame_smooth = nan(size(frame_k));
        valid = weight_blur > 1e-8;
        frame_smooth(valid) = frame_blur(valid) ./ weight_blur(valid);

        % Unless explicitly extrapolated, keep outside hull as NaN
        if ~fillOutsideHull
            frame_smooth(~HullMask_MB) = nan;
        end

        WavePropSim_3DMat(:,:,k) = frame_smooth;
    end
end

fprintf('Preparation complete.\n');
fprintf('Output WavePropSim_3DMat size: [%d x %d x %d]\n', ...
    size(WavePropSim_3DMat,1), size(WavePropSim_3DMat,2), size(WavePropSim_3DMat,3));

%% ============================== Diagnostics =============================
figure('Name','Multimeridian MB Preparation for 2D Phase Gradient','Units','normalized','Position',[0.05 0.08 0.85 0.6]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(ValidMask_MB);
axis image;
title('Original Sampled Mask');

nexttile;
imagesc(HullMask_MB);
axis image;
title('Interpolation Hull Mask');

nexttile;
imagesc(SampledFrames_MB(:,:,1));
axis image; colorbar;
title('Sparse Frame 1');

nexttile;
imagesc(WavePropSim_3DMat(:,:,1));
axis image; colorbar;
title('Prepared Frame 1');

midIdx = round(Nt/2);
nexttile;
imagesc(SampledFrames_MB(:,:,midIdx));
axis image; colorbar;
title(sprintf('Sparse Frame %d', midIdx));

nexttile;
imagesc(WavePropSim_3DMat(:,:,midIdx));
axis image; colorbar; colormap turbo;
title(sprintf('Prepared Frame %d', midIdx));

end


function InterpOperator_MB = localBuildInterpOperator(ValidMask_MB, interpMethod, fillOutsideHull)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% localBuildInterpOperator
% Build a fixed spatial interpolation operator for MB sparse geometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nx,Ny] = size(ValidMask_MB);

[Xq,Yq] = ndgrid(1:Nx,1:Ny);

valid_idx = find(ValidMask_MB);
[xv,yv] = ind2sub([Nx,Ny], valid_idx);

P_valid = [xv(:), yv(:)];
P_query = [Xq(:), Yq(:)];

switch interpMethod
    case 'linearDT'
        DT = delaunayTriangulation(P_valid);

        tri_id = pointLocation(DT, P_query);
        insideHull = ~isnan(tri_id);
        HullMask = reshape(insideHull, Nx, Ny);

        query_inside = P_query(insideHull,:);
        tri_inside   = tri_id(insideHull);

        bary = cartesianToBarycentric(DT, tri_inside, query_inside);

        T = DT.ConnectivityList;
        vertex_ids = T(tri_inside,:);

        if fillOutsideHull
            nearest_ids = knnsearch(P_valid, P_query(~insideHull,:));
        else
            nearest_ids = [];
        end

        InterpOperator_MB = struct();
        InterpOperator_MB.method = interpMethod;
        InterpOperator_MB.Nx = Nx;
        InterpOperator_MB.Ny = Ny;
        InterpOperator_MB.valid_idx = valid_idx;
        InterpOperator_MB.HullMask = HullMask;
        InterpOperator_MB.insideHull = insideHull;
        InterpOperator_MB.vertex_ids = vertex_ids;
        InterpOperator_MB.bary = bary;
        InterpOperator_MB.fillOutsideHull = fillOutsideHull;
        InterpOperator_MB.nearest_ids = nearest_ids;

    case 'nearest'
        nearest_ids_all = knnsearch(P_valid, P_query);

        InterpOperator_MB = struct();
        InterpOperator_MB.method = interpMethod;
        InterpOperator_MB.Nx = Nx;
        InterpOperator_MB.Ny = Ny;
        InterpOperator_MB.valid_idx = valid_idx;
        InterpOperator_MB.HullMask = true(Nx,Ny);
        InterpOperator_MB.nearest_ids_all = nearest_ids_all;
        InterpOperator_MB.fillOutsideHull = true;
end

end


function frame_interp = localApplyInterpOperator(frame_sparse, InterpOperator_MB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% localApplyInterpOperator
% Apply the precomputed interpolation operator to one sparse frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx = InterpOperator_MB.Nx;
Ny = InterpOperator_MB.Ny;

frame_interp = nan(Nx,Ny);

v_valid = frame_sparse(InterpOperator_MB.valid_idx);

switch InterpOperator_MB.method
    case 'linearDT'
        insideHull = InterpOperator_MB.insideHull;
        vertex_ids = InterpOperator_MB.vertex_ids;
        bary       = InterpOperator_MB.bary;

        vq_inside = ...
            bary(:,1).*v_valid(vertex_ids(:,1)) + ...
            bary(:,2).*v_valid(vertex_ids(:,2)) + ...
            bary(:,3).*v_valid(vertex_ids(:,3));

        frame_interp(insideHull) = vq_inside;

        if InterpOperator_MB.fillOutsideHull
            outsideHull = ~insideHull;
            nearest_ids = InterpOperator_MB.nearest_ids;
            frame_interp(outsideHull) = v_valid(nearest_ids);
        end

    case 'nearest'
        nearest_ids_all = InterpOperator_MB.nearest_ids_all;
        frame_interp(:) = v_valid(nearest_ids_all);
end

end