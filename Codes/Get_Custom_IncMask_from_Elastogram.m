function [Custom_IncMask, Custom_BackgroundMask, EdgeMask, DebugInfo] = ...
    Get_Custom_IncMask_from_Elastogram( ...
    Elastogram_Masked, Inc_Mask_Ideal, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get_Custom_IncMask_from_Elastogram
%
% Builds:
%
%   1) A custom reconstructed inclusion ROI mask
%   2) A custom background mask
%
% from a masked elastogram.
%
% -------------------------------------------------------------------------
% Key idea
% -------------------------------------------------------------------------
%
% The reconstructed inclusion geometry does not necessarily coincide with
% the ideal/theoretical inclusion mask because of:
%
%   - reflections,
%   - angular undersampling,
%   - interpolation artifacts,
%   - wavefront distortion,
%   - windowing artifacts.
%
% Therefore:
%
%   - the ideal mask is used ONLY as a spatial prior,
%   - while the final ROI is inferred from the elastogram itself.
%
% -------------------------------------------------------------------------
% Background mask logic
% -------------------------------------------------------------------------
%
% The custom background mask is obtained as:
%
%       Background = NOT(Inclusion)
%
% but leaving a transition band around the reconstructed ROI so that:
%
%   - edge-transition pixels,
%   - interpolation boundaries,
%   - partial-volume regions,
%   - wavefront discontinuities
%
% are NOT included in the background statistics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse inputs ===============================

P = inputParser;

addParameter(P, 'SearchDilateRadius_px', 35);
addParameter(P, 'BoundaryBandRadius_px', 28);

addParameter(P, 'CannySigma', 1.6);
addParameter(P, 'CannyThreshold', []);

addParameter(P, 'DisplayCLim', [4 12]);

addParameter(P, 'CloseRadius_px', 10);
addParameter(P, 'OpenRadius_px', 3);

addParameter(P, 'MinArea_px', 80);

addParameter(P, 'LowSpeedQuantile', 0.45);

addParameter(P, 'RequireOverlapWithIdealPrior', true);

% NEW:
% Radius used to exclude transition pixels from the background mask.
addParameter(P, 'BackgroundGapRadius_px', 12);

addParameter(P, 'DisplayDiagnostics', true);

parse(P, varargin{:});
Opt = P.Results;

Elastogram_Masked = double(Elastogram_Masked);
Inc_Mask_Ideal    = logical(Inc_Mask_Ideal);

if ~isequal(size(Elastogram_Masked), size(Inc_Mask_Ideal))
    error('Elastogram_Masked and Inc_Mask_Ideal must have the same size.');
end

[Nx, Ny] = size(Elastogram_Masked);

%% ====================== Spatial prior from ideal mask ====================

searchMask = imdilate(Inc_Mask_Ideal, ...
    strel('disk', Opt.SearchDilateRadius_px));

searchMask = searchMask & isfinite(Elastogram_Masked);

idealBoundary = bwperim(Inc_Mask_Ideal);

distToIdealBoundary = bwdist(idealBoundary);

boundaryBand = distToIdealBoundary <= Opt.BoundaryBandRadius_px;

boundaryBand = boundaryBand & ...
    imdilate(Inc_Mask_Ideal, ...
    strel('disk', Opt.SearchDilateRadius_px));

boundaryBand = boundaryBand & isfinite(Elastogram_Masked);

%% ====================== CLim-based normalization =========================

climLow  = Opt.DisplayCLim(1);
climHigh = Opt.DisplayCLim(2);

E = Elastogram_Masked;

E_clip = E;

E_clip(E_clip < climLow)  = climLow;
E_clip(E_clip > climHigh) = climHigh;

E_norm = (E_clip - climLow) ./ (climHigh - climLow + eps);

E_norm(~isfinite(E_norm)) = 0;

E_norm_search = E_norm;
E_norm_search(~searchMask) = 0;

%% ============================= Edge detection ============================

if isempty(Opt.CannyThreshold)

    EdgeMask_raw = edge( ...
        E_norm_search, ...
        'canny', ...
        [], ...
        Opt.CannySigma);

else

    EdgeMask_raw = edge( ...
        E_norm_search, ...
        'canny', ...
        Opt.CannyThreshold, ...
        Opt.CannySigma);
end

EdgeMask = EdgeMask_raw & boundaryBand;

EdgeMask = bwareaopen(EdgeMask, 5);

EdgeMask_closed = imclose(EdgeMask, strel('disk', 2));

EdgeMask_closed = bwmorph(EdgeMask_closed, 'bridge');
EdgeMask_closed = bwmorph(EdgeMask_closed, 'clean');

%% ====================== Low-speed reconstructed ROI ======================

validLocalVals = E_norm(searchMask);

q = Opt.LowSpeedQuantile;
q = min(max(q, 0.05), 0.95);

lowSpeedThreshold = quantile(validLocalVals, q);

LowSpeedMask = false(Nx, Ny);

LowSpeedMask(searchMask) = ...
    E_norm(searchMask) <= lowSpeedThreshold;

priorForRegion = imdilate( ...
    Inc_Mask_Ideal, ...
    strel('disk', Opt.SearchDilateRadius_px));

LowSpeedMask = LowSpeedMask & priorForRegion;

%% ====================== Morphological ROI reconstruction =================

closeSE = strel('disk', Opt.CloseRadius_px);
openSE  = strel('disk', Opt.OpenRadius_px);

edgeSupport = imdilate( ...
    EdgeMask_closed, ...
    strel('disk', max(2, round(Opt.CloseRadius_px/3))));

candidateMask = LowSpeedMask | edgeSupport;

candidateMask = candidateMask & priorForRegion;

candidateMask = imclose(candidateMask, closeSE);

candidateMask = imfill(candidateMask, 'holes');

candidateMask = imopen(candidateMask, openSE);

candidateMask = bwareaopen(candidateMask, Opt.MinArea_px);

%% ====================== Select best connected component ==================

CC = bwconncomp(candidateMask);

if CC.NumObjects == 0

    warning('No valid custom ROI detected. Returning ideal mask.');

    Custom_IncMask = Inc_Mask_Ideal;

else

    stats = regionprops( ...
        CC, ...
        'Area', ...
        'PixelIdxList', ...
        'Centroid');

    idealStats = regionprops( ...
        Inc_Mask_Ideal, ...
        'Area', ...
        'Centroid');

    if isempty(idealStats)

        idealArea = nnz(Inc_Mask_Ideal);
        idealCentroidXY = [Ny/2, Nx/2];

    else

        idealArea = idealStats(1).Area;
        idealCentroidXY = idealStats(1).Centroid;
    end

    scores = inf(CC.NumObjects,1);

    for c = 1:CC.NumObjects

        comp = false(Nx, Ny);
        comp(stats(c).PixelIdxList) = true;

        overlapIdeal = ...
            nnz(comp & Inc_Mask_Ideal) / ...
            max(nnz(comp | Inc_Mask_Ideal), 1);

        overlapPrior = ...
            nnz(comp & priorForRegion) / ...
            max(nnz(comp), 1);

        areaPenalty = ...
            abs(stats(c).Area - idealArea) / ...
            max(idealArea, 1);

        centroidDist = norm( ...
            stats(c).Centroid - idealCentroidXY);

        centroidPenalty = centroidDist / max([Nx, Ny]);

        edgeContact = ...
            nnz(bwperim(comp) & ...
            imdilate(EdgeMask_closed, strel('disk',2))) / ...
            max(nnz(bwperim(comp)),1);

        if Opt.RequireOverlapWithIdealPrior && overlapPrior < 0.80
            priorPenalty = 10;
        else
            priorPenalty = 0;
        end

        scores(c) = ...
            -2.5*overlapIdeal ...
            -1.5*edgeContact ...
            +0.9*areaPenalty ...
            +0.7*centroidPenalty ...
            +priorPenalty;
    end

    [~, bestIdx] = min(scores);

    Custom_IncMask = false(Nx, Ny);

    Custom_IncMask(stats(bestIdx).PixelIdxList) = true;

    Custom_IncMask = imclose(Custom_IncMask, closeSE);

    Custom_IncMask = imfill(Custom_IncMask, 'holes');

    Custom_IncMask = imopen(Custom_IncMask, openSE);

    Custom_IncMask = bwareaopen( ...
        Custom_IncMask, ...
        Opt.MinArea_px);

    Custom_IncMask = Custom_IncMask & priorForRegion;
end

%% ====================== Build custom background mask =====================
% The background is defined as:
%
%   Background = NOT(Custom_IncMask)
%
% but excluding:
%   1) a transition band around the reconstructed inclusion ROI,
%   2) the lower-right / 4th-quadrant region of the elastogram,
%      where strong reflection artifacts often appear.

transitionBand = imdilate( ...
    Custom_IncMask, ...
    strel('disk', Opt.BackgroundGapRadius_px));

% -------------------------------------------------------------
% Exclude 4th quadrant / lower-right image region from background
% -------------------------------------------------------------
% In matrix/image coordinates:
%   rows increase downward,
%   columns increase to the right.
%
% Therefore, the lower-right quadrant is:
%
%   row > centerRow
%   col > centerCol
%
% This matches the artifact-heavy region seen near the reconstructed
% inclusion in the lower-right part of the elastogram.

[Nx, Ny] = size(Elastogram_Masked);

centerRow = round(Nx/2);
centerCol = round(Ny/2);

[RR, CC] = ndgrid(1:Nx, 1:Ny);

FourthQuadrantMask = (RR > centerRow) & (CC > centerCol);

Custom_BackgroundMask = ...
    ~transitionBand & ...
    ~FourthQuadrantMask & ...
    isfinite(Elastogram_Masked);

%% ============================= Debug info ================================

DebugInfo = struct();

DebugInfo.Options = Opt;

DebugInfo.SearchMask = searchMask;

DebugInfo.IdealBoundary = idealBoundary;

DebugInfo.BoundaryBand = boundaryBand;

DebugInfo.E_clip = E_clip;

DebugInfo.E_norm = E_norm;

DebugInfo.EdgeMask_raw = EdgeMask_raw;

DebugInfo.EdgeMask = EdgeMask;

DebugInfo.EdgeMask_closed = EdgeMask_closed;

DebugInfo.LowSpeedMask = LowSpeedMask;

DebugInfo.CandidateMask = candidateMask;

DebugInfo.transitionBand = transitionBand;

DebugInfo.FourthQuadrantMask = FourthQuadrantMask;

%% ============================= Diagnostics ===============================

if Opt.DisplayDiagnostics

    figure( ...
        'Name', 'Custom Inclusion ROI Detection', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.08 0.88 0.72]);

    tiledlayout(2,3, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    %% 1) Input masked elastogram
    nexttile;

    imagesc(Elastogram_Masked);
    axis image;
    colormap(gca, turbo);
    clim(Opt.DisplayCLim);
    colorbar;

    title('Input Masked Elastogram');

    hold on;
    contour(Inc_Mask_Ideal, [0.5 0.5], ...
        'w', 'LineWidth', 2);
    hold off;

    %% 2) CLim-normalized elastogram
    nexttile;

    imagesc(E_norm_search);
    axis image;
    colormap(gca, turbo);
    colorbar;

    title('CLim-Normalized Elastogram');

    %% 3) Candidate ROI
    nexttile;

    imagesc(candidateMask);
    axis image;
    colormap(gca, turbo);
    colorbar;

    title('Candidate ROI');

    %% 4) Final custom ROI overlay
    nexttile;

    imagesc(Elastogram_Masked);
    axis image;
    colormap(gca, turbo);
    clim(Opt.DisplayCLim);
    colorbar;

    title('Final Custom ROI Overlay');

    hold on;
    contour(Inc_Mask_Ideal, [0.5 0.5], ...
        'w--', 'LineWidth', 1.8);

    contour(Custom_IncMask, [0.5 0.5], ...
        'Color', [0.45 0 0.85], ...
        'LineWidth', 3.0);
    hold off;

    %% 5) Custom background mask
    nexttile;

    imagesc(Custom_BackgroundMask);
    axis image;
    colormap(gca, turbo);
    colorbar;

    title('Custom Background Mask');

end

end