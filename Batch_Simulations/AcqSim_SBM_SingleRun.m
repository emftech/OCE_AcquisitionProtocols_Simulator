function Results = AcqSim_SBM_SingleRun(cfg, W, nBscans, run_tag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcqSim_SBM_SingleRun
%
% Function version of Sampling_Acq_Sim_SBM_v4_Cylindrical_Ideal.mlx
% (sections "Definition of Acquisition Vectors" to "Store Metrics"),
% applied to a precomputed causal harmonic video W (AcqSim_WaveProp,
% protocol "SBM"). Same framework calls and arguments as the live script:
%
%   AcqVectors_SBM -> SamplingSim_fromAcqVecs -> Recovery_fromSampledFrames_SBM
%   -> Prepare_SBM_2DPhaseGradient -> PhaseGradient2D_FromSims_v3 /
%   PhaseGradient2D_TV_FromSims (Ts_SBM = phaseStepSamples*Ts_OCT)
%
% Protocol parameters follow the SBM log (19-10-2026): N_A = 37,
% N_F = P = T/Ts, N_P = 4 -> phaseStepSamples = N_phi = N_P*P/N_F = 4,
% N_turn = smallest k*P - N_A >= 20, startDelay = P - ExcPulseWidth so the
% raster starts at phase 0 of the excitation ("causal").
%
% Note: the thesis SBM chain has no 1D estimator (1D PG is defined along MB
% meridians), so SBM is evaluated with 2D PG and 2D PG-TV, as in the thesis.
%
% Storage (thesis-like, no saccades):
%   <results>/SBM/<NA>_Ascans_<NB>_BScans_<NT>_TurnSamples/Harmonic_Ideal_Cylindrical_<T>us/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4, run_tag = ""; end
t_run = tic;

T_excitation = W.T_excitation;
T_us   = round(T_excitation*1e6);
Ts_OCT = cfg.Ts_OCT;
P      = round(T_excitation/Ts_OCT);
hMin = cfg.horizontalMin; hMax = cfg.horizontalMax;
vMin = cfg.verticalMin;   vMax = cfg.verticalMax;
GT_SpeedMap = W.GT_SpeedMap; Bkg_Mask = W.Bkg_Mask; Inc_Mask = W.Inc_Mask; Ann_Mask = W.Ann_Mask;

%% SBM protocol parameters (log rules)
nAscans = cfg.SBM.nAscans;
nFrames = P;                                   % N_F = P
phaseStepSamples = cfg.SBM.N_P * P / nFrames;  % N_phi = 4
nTurnSamples = ceil((nAscans + cfg.SBM.MinTurn)/P)*P - nAscans;
ExcPulseWidth = cfg.SBM.ExcPulseWidth;
AcqPulseWidth = cfg.SBM.AcqPulseWidth;
startDelay_samples = P - ExcPulseWidth;        % thesis: 47 = 57 - 10

%% Output folder
exp_type = "Harmonic_Ideal_Cylindrical";
prop_type_medium = 'Cylindrical';
pattern_path = sprintf("%d_Ascans_%d_BScans_%d_TurnSamples", nAscans, nBscans, nTurnSamples);
run_path = fullfile(cfg.results_path, "SBM", pattern_path, sprintf("%s_%dus", exp_type, T_us));
if ~isfolder(run_path), mkdir(run_path); end
tagFile = sprintf("RasterSBM_%dB_%dA_%dT", nBscans, nAscans, nTurnSamples);

%% Acquisition vectors (SBM)
[Exc_Trigger, Acq_Trigger, GalvoX, GalvoY, Exc_aux] = AcqVectors_SBM(hMin, hMax, vMin, ...
    vMax, nAscans, nBscans, nFrames, ...
    startDelay_samples, nTurnSamples, AcqPulseWidth, ...
    ExcPulseWidth, T_excitation, Ts_OCT, phaseStepSamples);
close all force;

%% Sampling (same output assignment as the thesis live script)
[SampledFrames_vis, SampledFrames_raw] = SamplingSim_fromAcqVecs(W.Vid, W.x_grid, ...
    W.y_grid, W.t_grid, GalvoX, GalvoY, Acq_Trigger, Exc_aux, nFrames, nAscans, nBscans, ...
    nTurnSamples, T_excitation, Ts_OCT);

%% Recovery (Doppler / accumulated-phase compensation)
[Acq_Phase_Matrix, Recovered_Frames] = Recovery_fromSampledFrames_SBM( ...
    nTurnSamples, T_excitation, Ts_OCT, nFrames, SampledFrames_raw, phaseStepSamples);
Recovered_Frames = flipud(Recovered_Frames);

%% Interpolation to a square grid
[WavePropSim_3DMat, ValidMask_SBM, InterpOperator_SBM] = ...
    Prepare_SBM_2DPhaseGradient(Recovered_Frames, 'spline', false, true);
close all force;
Inc_Mask_ReducedGrid = imresize(Inc_Mask, size(WavePropSim_3DMat, 1:2));
Bkg_Mask_ReducedGrid = imresize(Bkg_Mask, size(WavePropSim_3DMat, 1:2));
Ts_SBM = Ts_OCT*phaseStepSamples;

%% Near-field blanking BEFORE the estimators (experimental Cfg.NeedleRadius_mm rule)
% Prepare_SBM_2DPhaseGradient interpolates with griddedInterpolant 'spline',
% which turns the whole frame NaN if a hole is present, so the disk is
% removed on the interpolated movie (rows = y, cols = x). Difference with
% the experimental SBM chain: the spline still uses the near-field samples
% to build the frame outside the disk. Windows of the PG kernel that touch
% the hole are skipped by the estimators (unwrap/fit -> NaN).
R_nf = cfg.NeedleRadius_mm; Src = cfg.Source_mm;
xi = linspace(hMin, hMax, size(WavePropSim_3DMat,2));
yi = linspace(vMin, vMax, size(WavePropSim_3DMat,1));
[Xi, Yi] = meshgrid(xi, yi);
NearI = hypot(Xi - Src(1), Yi - Src(2)) < R_nf;
WavePropSim_3DMat_Unblanked = WavePropSim_3DMat;            % display only (full-field video)
WavePropSim_3DMat(repmat(NearI, 1, 1, size(WavePropSim_3DMat,3))) = NaN;

%% 2D PG
[Elastogram_SBM_2DPG, R2_map_SBM_2DPG, Elastogram_Masked_SBM_2DPG, ...
 Xcrop_mm_SBM_2DPG, Ycrop_mm_SBM_2DPG, Complex_Field_Vz_SBM_2DPG, ...
 Inc_Mask_Crop_SBM_2DPG, Bkg_Mask_Crop_SBM_2DPG, Kx_SBM_2DPG, Ky_SBM_2DPG, Kmag_SBM_2DPG] = ...
    PhaseGradient2D_FromSims_v3( ...
        WavePropSim_3DMat, Inc_Mask_ReducedGrid, Bkg_Mask_ReducedGrid, ...
        cfg.kernel_size, Ts_SBM, T_excitation, hMin, hMax, vMin, vMax, ...
        cfg.applyTemporalWindow, cfg.clim);
close all force;

%% 2D PG-TV
[Elastogram_SBM_2DPG_TV, R2_map_SBM_2DPG_TV, Elastogram_Masked_SBM_2DPG_TV, ...
 Xcrop_mm_SBM_2DPG_TV, Ycrop_mm_SBM_2DPG_TV, Complex_Field_Vz_SBM_2DPG_TV, ...
 Inc_Mask_Crop_SBM_2DPG_TV, Bkg_Mask_Crop_SBM_2DPG_TV, Kx_SBM_2DPG_TV, Ky_SBM_2DPG_TV, Kmag_SBM_2DPG_TV] = ...
    PhaseGradient2D_TV_FromSims( ...
        WavePropSim_3DMat, Inc_Mask_ReducedGrid, Bkg_Mask_ReducedGrid, ...
        cfg.kernel_size, Ts_SBM, T_excitation, hMin, hMax, vMin, vMax, ...
        cfg.applyTemporalWindow, cfg.clim, ...
        'MuTV', cfg.MuTV, 'TVIterations', cfg.TVIterations, 'TVTau', cfg.TVTau, ...
        'TVSigma', cfg.TVSigma, 'TVTheta', cfg.TVTheta, 'MinK', cfg.MinK, 'MaxK', cfg.MaxK, ...
        'R2Threshold', cfg.R2Threshold, 'DisplayFigures', false, 'Verbose', false);
close all force;

%% Remove the same disk from every map (NearSource), as in the experiments
Elastogram_SBM_2DPG           = blankDisk(Elastogram_SBM_2DPG,           Xcrop_mm_SBM_2DPG,    Ycrop_mm_SBM_2DPG,    Src, R_nf);
Elastogram_Masked_SBM_2DPG    = blankDisk(Elastogram_Masked_SBM_2DPG,    Xcrop_mm_SBM_2DPG,    Ycrop_mm_SBM_2DPG,    Src, R_nf);
Elastogram_SBM_2DPG_TV        = blankDisk(Elastogram_SBM_2DPG_TV,        Xcrop_mm_SBM_2DPG_TV, Ycrop_mm_SBM_2DPG_TV, Src, R_nf);
Elastogram_Masked_SBM_2DPG_TV = blankDisk(Elastogram_Masked_SBM_2DPG_TV, Xcrop_mm_SBM_2DPG_TV, Ycrop_mm_SBM_2DPG_TV, Src, R_nf);

%% Metrics
Est = struct();
Est.SBM_2DPG    = struct('E', Elastogram_Masked_SBM_2DPG,    'x', Xcrop_mm_SBM_2DPG,    'y', Ycrop_mm_SBM_2DPG);
Est.SBM_2DPG_TV = struct('E', Elastogram_Masked_SBM_2DPG_TV, 'x', Xcrop_mm_SBM_2DPG_TV, 'y', Ycrop_mm_SBM_2DPG_TV);
Metrics = AcqSim_ComputeMetrics(Est, GT_SpeedMap, Inc_Mask, Bkg_Mask, Ann_Mask, W.Far_Mask);

%% Store (thesis Results struct + batch extras)
Results = struct();
Results.Metadata = struct('exp_type', exp_type, 'protocol_type', "SBM", 'prop_type_medium', prop_type_medium, ...
    'nAscans', nAscans, 'nBscans', nBscans, 'nFrames', nFrames, 'nTurnSamples', nTurnSamples, ...
    'kernel_size', cfg.kernel_size, 'applyTemporalWindow', cfg.applyTemporalWindow, ...
    'clim_limits_2DPG', cfg.clim, 'clim_limits_PGTV', cfg.clim, ...
    'date_saved', string(datetime("now")), 'run_tag', run_tag);
Results.Parameters = struct('c_medium', cfg.c_ref, 'Ts_OCT', Ts_OCT, 'T_excitation', T_excitation, ...
    'f_excitation', 1/T_excitation, 'P', P, 'N_P', cfg.SBM.N_P, 'phaseStepSamples', phaseStepSamples, ...
    'Ts_SBM', Ts_SBM, 'startDelay_samples', startDelay_samples, 'ExcPulseWidth', ExcPulseWidth, ...
    'AcqPulseWidth', AcqPulseWidth, 'alpha', cfg.alpha, 'N_excitations', cfg.SBM.N_excitations, ...
    'A_Peak', cfg.A_Peak, 'horizontalMin', hMin, 'horizontalMax', hMax, 'verticalMin', vMin, ...
    'verticalMax', vMax, 'nPoints', cfg.nPoints, 'prop_type', W.prop_type, 'MotionCase', "ideal", ...
    'AcqFramework', "SBM", 'SourceFWHM_mm', cfg.SourceFWHM_mm, 'NeedleExclusion_mm', cfg.NeedleExclusion_mm, ...
    'CausalClip', W.clipInfo, 'TotalAcqTime_ms', numel(GalvoX)*Ts_OCT*1e3, ...
    'NearFieldBlanking', struct('Radius_mm', R_nf, 'Source_mm', Src, ...
        'Rule', "interpolated (100x100, spline) movie set to NaN where distance to the tip < R before 2D PG / 2D PG-TV; same disk set to NaN on every map", ...
        'DifferenceToExperiment', "the spline interpolation of the recovered frames still uses the near-field samples (blanking is applied after it)"));
Results.Acquisition = struct('GalvoX', GalvoX, 'GalvoY', GalvoY, 'Acq_Trigger', Acq_Trigger, ...
    'Exc_Trigger', Exc_Trigger, 'Exc_aux', Exc_aux);
Results.SampledData = struct('SampledFrames_vis', SampledFrames_vis, 'SampledFrames_raw', SampledFrames_raw, ...
    'Acq_Phase_Matrix', Acq_Phase_Matrix, 'Recovered_Frames', Recovered_Frames, ...
    'WavePropSim_3DMat', WavePropSim_3DMat, 'ValidMask_SBM', ValidMask_SBM, 'InterpOperator_SBM', InterpOperator_SBM);
Results.GroundTruth = struct('GT_SpeedMap', GT_SpeedMap, 'Inc_Mask', Inc_Mask, 'Bkg_Mask', Bkg_Mask, ...
    'Ann_Mask', Ann_Mask, 'Inc_Mask_ReducedGrid', Inc_Mask_ReducedGrid, 'Bkg_Mask_ReducedGrid', Bkg_Mask_ReducedGrid);
Results.SBM_2DPG = struct('Elastogram', Elastogram_SBM_2DPG, 'Elastogram_Masked', Elastogram_Masked_SBM_2DPG, ...
    'R2_map', R2_map_SBM_2DPG, 'Xaxis_mm', Xcrop_mm_SBM_2DPG, 'Yaxis_mm', Ycrop_mm_SBM_2DPG, ...
    'Complex_Field_Vz', Complex_Field_Vz_SBM_2DPG, 'Inc_Mask_Crop', Inc_Mask_Crop_SBM_2DPG, ...
    'Bkg_Mask_Crop', Bkg_Mask_Crop_SBM_2DPG, 'Kx', Kx_SBM_2DPG, 'Ky', Ky_SBM_2DPG, 'Kmag', Kmag_SBM_2DPG);
Results.SBM_2DPG_TV = struct('Elastogram', Elastogram_SBM_2DPG_TV, 'Elastogram_Masked', Elastogram_Masked_SBM_2DPG_TV, ...
    'R2_map', R2_map_SBM_2DPG_TV, 'Xaxis_mm', Xcrop_mm_SBM_2DPG_TV, 'Yaxis_mm', Ycrop_mm_SBM_2DPG_TV, ...
    'Complex_Field_Vz', Complex_Field_Vz_SBM_2DPG_TV, 'Inc_Mask_Crop', Inc_Mask_Crop_SBM_2DPG_TV, ...
    'Bkg_Mask_Crop', Bkg_Mask_Crop_SBM_2DPG_TV, 'Kx', Kx_SBM_2DPG_TV, 'Ky', Ky_SBM_2DPG_TV, 'Kmag', Kmag_SBM_2DPG_TV);
fn = fieldnames(Metrics);
for k = 1:numel(fn)
    Results.(fn{k}).AccuracyError  = Metrics.(fn{k}).AccuracyError;
    Results.(fn{k}).PrecisionError = Metrics.(fn{k}).PrecisionError;
    Results.(fn{k}).Annulus        = Metrics.(fn{k}).Annulus;
    Results.(fn{k}).FarField       = Metrics.(fn{k}).FarField;
end
Results.OG_WaveProp  = W.OG;
Results.SummaryTable = AcqSim_SummaryTable(Metrics);

results_filename = sprintf("SBM_SimResults_%d_Ascans_%d_BScans_%d_TurnSamples_%s_%s_%dus.mat", ...
    nAscans, nBscans, nTurnSamples, exp_type, prop_type_medium, T_us);
Results.Metadata.results_filename = results_filename;
Results.Metadata.results_fullpath = fullfile(run_path, results_filename);
save(fullfile(run_path, results_filename), "Results", "-v7.3");

%% Elastogram figures (thesis naming)
names = {'SBM_2DPG','2DPG'; 'SBM_2DPG_TV','2DPGTV'};
for k = 1:size(names,1)
    e = Est.(names{k,1});
    f = sprintf("Elastogram_%s_%s_Masked_R2gt0p7_%dusExc_%s", names{k,2}, tagFile, T_us, exp_type);
    AcqSim_SaveElastogram(e.E, e.x, e.y, cfg.clim, ...
        {sprintf('SBM %d A x %d B (%d turn) - %s - %d \\mus', nAscans, nBscans, nTurnSamples, names{k,2}, T_us), ...
         sprintf('mean %.2f m/s (annulus r \\geq %.1f mm)', Metrics.(names{k,1}).Annulus.mean, cfg.NeedleExclusion_mm)}, ...
        fullfile(run_path, f), Src);
end

%% Videos / figures shown in the live script (NaN copies for display only: gray disk)
xs = linspace(hMin, hMax, nAscans);
y_desc = linspace(vMax, vMin, nBscans);   % acquisition order: first B-scan at y = +FOV/2
y_asc  = linspace(vMin, vMax, nBscans);
[Xd, Yd] = meshgrid(xs, y_desc); NearDesc = hypot(Xd - Src(1), Yd - Src(2)) < R_nf;
[Xa, Ya] = meshgrid(xs, y_asc);  NearAsc  = hypot(Xa - Src(1), Ya - Src(2)) < R_nf;
vopt = struct('FrameStep', 1, 'FrameRate', 10, 'dt_s', Ts_SBM);
tagV = sprintf("%s_%dus_%s", tagFile, T_us, exp_type);
sfx = sprintf(' (r < %.2f mm blanked)', R_nf);
AcqSim_WriteVideo(SampledFrames_vis, xs, y_desc, ...
    sprintf('Corrected-orientation sampled frames - SBM - %d \\mus%s', T_us, sfx), fullfile(run_path, "SampledVisFrames_" + tagV), vmask(vopt, NearDesc));
AcqSim_WriteVideo(SampledFrames_raw, xs, y_desc, ...
    sprintf('Raw sampled frames (no serpentine compensation) - SBM - %d \\mus%s', T_us, sfx), fullfile(run_path, "SampledRawFrames_" + tagV), vmask(vopt, NearDesc));
Recovered_Blanked = nanDisk(Recovered_Frames, NearAsc);
AcqSim_WriteVideo(Recovered_Frames, xs, y_asc, ...
    sprintf('Recovered frames (acquisition-phase compensated) - %d \\mus%s', T_us, sfx), fullfile(run_path, "RecoveredFrames_" + tagV), vmask(vopt, NearAsc));
AcqSim_WriteVideo(WavePropSim_3DMat_Unblanked, xi, yi, ...
    sprintf('Interpolated frames for 2D PG - SBM - %d \\mus%s', T_us, sfx), fullfile(run_path, "InterpFrames_" + tagV), vmask(vopt, NearI));
f = figure('Visible','off','Color','w','Position',[100 100 1500 420]);
tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
ax = nexttile; imagesc(ax, xs, y_desc, Acq_Phase_Matrix); set(ax,'YDir','normal'); axis(ax,'image');
colormap(ax, hsv); clim(ax, [-pi pi]); colorbar(ax);
title(ax, 'Accumulated acquisition phase', 'FontWeight', 'bold'); xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');
ax = nexttile; iy = ceil(nBscans/2); ST = squeeze(Recovered_Blanked(iy,:,:)).'; lim = prctile(abs(ST(isfinite(ST))),97);
im = imagesc(ax, xs, (0:nFrames-1)*Ts_SBM*1e3, ST); set(im, 'AlphaData', isfinite(ST)); set(ax, 'Color', [0.6 0.6 0.6]);
colormap(ax, MMB_RedBlue(255)); clim(ax,[-lim lim]); colorbar(ax);
xlabel(ax,'x [mm]'); ylabel(ax,'equivalent time [ms]');
title(ax, {'S-T map of recovered frames at y = 0', sprintf('near field r < %.2f mm in gray', R_nf)}, 'FontWeight', 'bold');
ax = nexttile; STf = squeeze(Recovered_Frames(iy,:,:)).';
im = imagesc(ax, xs, (0:nFrames-1)*Ts_SBM*1e3, STf); set(im, 'AlphaData', isfinite(STf)); set(ax, 'Color', [0.6 0.6 0.6]);
colormap(ax, MMB_RedBlue(255)); clim(ax,[-lim lim]); colorbar(ax);
xlabel(ax,'x [mm]'); ylabel(ax,'equivalent time [ms]');
title(ax, {'Full field (same colour limits:', 'near field saturated)'}, 'FontWeight', 'bold');
exportgraphics(f, fullfile(run_path, "AcqPhase_and_ST_" + tagV + ".png"), 'Resolution', 150);
savefig(f, fullfile(run_path, "AcqPhase_and_ST_" + tagV + ".fig")); close(f);

fprintf('[SBM] %d A x %d B x %d T, N_F = %d @ %d us done in %.1f s -> %s\n', ...
    nAscans, nBscans, nTurnSamples, nFrames, T_us, toc(t_run), run_path);
end

%% ------------------------------------------------------------------------
function E = blankDisk(E, x_mm, y_mm, Src, R)
% NaN the disk |r - Src| < R on a map whose axes are x_mm / y_mm (works for
% [Ny x Nx] and [Nx x Ny] maps; the disk is symmetric for a centred source).
if isempty(E), return; end
if isequal(size(E), [numel(y_mm) numel(x_mm)])
    [X, Y] = meshgrid(x_mm, y_mm);
else
    [X, Y] = ndgrid(x_mm, y_mm);
end
E(hypot(X - Src(1), Y - Src(2)) < R) = NaN;
end

function V = nanDisk(V, Near)
% Display copy of a [Ny x Nx x Nt] movie with the near-field disk set to NaN.
V = double(V);
V(repmat(Near, 1, 1, size(V,3))) = NaN;
end

function o = vmask(o, M)
% Video options with the near-field mask (cropped + _Full videos).
o.NaNMask = M;
end
