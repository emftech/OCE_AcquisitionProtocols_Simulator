function W = AcqSim_WaveProp(cfg, T_excitation, protocol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcqSim_WaveProp
%
% Headless wrapper around Wave_Prop_Simulation_v3 (thesis framework) for
% one excitation period. The k-Wave video is computed ONCE per period and
% protocol and then reused by every acquisition pattern of that period.
%
%   protocol = "MB"  -> 'K_Wave_Quasi_Harmonic_CustomMedium' (3-cycle QH
%                       bursts, ensemble = M_repetitions + nSyncSamples),
%                       exactly as Sampling_Acq_Sim_MB_v4_Ideal_Cylindrical_QH.
%   protocol = "SBM" -> 'K_Wave_Harmonic_CustomMedium' (continuous sine).
%                       "Causal" excitation: the source runs before the
%                       acquisition, and the video is clipped at an absolute
%                       k-Wave sample that is (i) after the wave filled the
%                       cfg.FOV_mm domain and (ii) an integer multiple of
%                       P = T_exc/Ts, so video t = 0 is phase 0 of the
%                       steady-state excitation (scan starts at phase 0
%                       because startDelay = P - ExcPulseWidth).
%
% Output struct W: Vid, Vid_PG (clip used by the 2D PG on the full field),
% x_grid, y_grid, t_grid, GT_SpeedMap, Bkg_Mask, Inc_Mask, Ann_Mask and
% the 2D PG of the original (fully sampled) wave field as reference.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

protocol = upper(string(protocol));

%% Headless homogeneous medium (no inclusions -> cylindrical propagation)
MediumConfig = struct();
MediumConfig.background_speed = cfg.c_ref;
MediumConfig.density          = cfg.rho;
MediumConfig.apply_smoothing  = false;
MediumConfig.sigma_smooth_mm  = 0;
MediumConfig.inclusions = struct('type',{},'center_x_mm',{},'center_y_mm',{}, ...
    'size_x_mm',{},'size_y_mm',{},'angle_deg',{},'speed',{});

sim_opts = struct();
sim_opts.visualize             = false;
sim_opts.useGUI                = false;
sim_opts.saveConfig            = false;   % never touch the user's saved GUI medium
sim_opts.MediumConfig          = MediumConfig;
sim_opts.showSourceDiagnostics = false;
sim_opts.SourceFWHM_mm         = cfg.SourceFWHM_mm;

switch protocol
    case "MB"
        prop_type = 'K_Wave_Quasi_Harmonic_CustomMedium';
        sim_opts.M_repetitions = cfg.MB.M_repetitions;
        sim_opts.nSyncSamples  = cfg.MB.nSyncSamples;
        N_pulses      = cfg.MB.N_pulses;
        N_excitations = cfg.MB.N_excitations;
    case "SBM"
        prop_type = 'K_Wave_Harmonic_CustomMedium';
        sim_opts.M_repetitions = cfg.SBM.M_repetitions;
        sim_opts.nSyncSamples  = cfg.SBM.nSyncSamples;
        N_pulses      = 1;
        N_excitations = cfg.SBM.N_excitations;
    otherwise
        error('protocol must be "MB" or "SBM".');
end

[Vid, x_grid, y_grid, t_grid, GT_SpeedMap] = Wave_Prop_Simulation_v3( ...
    cfg.c_ref, cfg.Ts_OCT, T_excitation, cfg.alpha, N_pulses, cfg.A_Peak, ...
    cfg.horizontalMin, cfg.horizontalMax, cfg.verticalMin, cfg.verticalMax, ...
    cfg.nPoints, N_excitations, prop_type, sim_opts);

%% Protocol-specific clipping (thesis conventions)
P = round(T_excitation/cfg.Ts_OCT);
switch protocol
    case "MB"
        % Thesis QH script: 2D PG on the full field uses the second ensemble.
        Vid_PG = Vid(:,:,round(size(Vid,3)/2)+1:end);
        clipInfo = struct('type','QH second half for 2D PG on full field');
    case "SBM"
        % Absolute k-Wave sample of video index 1 (same T_transient as the
        % child solver: 0.185 ms).
        n_offset = round(0.185e-3/cfg.Ts_OCT);
        n_clip   = ceil(cfg.SBM.SteadyState_s/T_excitation) * P;    % multiple of P
        idx0     = n_clip - n_offset + 1;
        Vid      = Vid(:,:,idx0:end);
        t_grid   = t_grid(1:size(Vid,3));
        Vid_PG   = Vid;
        clipInfo = struct('type','harmonic causal clip', 'AbsoluteSample', n_clip, ...
            'AbsoluteTime_s', n_clip*cfg.Ts_OCT, 'VideoIndex', idx0, 'P', P);
end

%% Masks: thesis background disk + experiment-style annulus (needle blanking)
[Nx, Ny] = size(GT_SpeedMap);
[xIdx, yIdx] = ndgrid(1:Nx, 1:Ny);
cx = (Nx + 1)/2; cy = (Ny + 1)/2;
R  = min(Nx, Ny)/2 - cfg.MarginPx;
Bkg_Mask = ((xIdx - cx).^2 + (yIdx - cy).^2) <= R^2;
Inc_Mask = false(size(GT_SpeedMap));
[Xm, Ym] = ndgrid(x_grid*1e3, y_grid*1e3);
Ann_Mask = Bkg_Mask & (hypot(Xm, Ym) >= cfg.NeedleExclusion_mm);
% Far field: outside the blanked near field (source support) plus the PG-kernel
% half-width margin (edge-free ROI). Previously: outside the compact support of the k-Wave source
% (G > exp(-4) <=> r < 2*sqrt(2)*sigma_s), where the field is no longer forced.
sigma_s_mm = cfg.SourceFWHM_mm/(2*sqrt(2*log(2)));
FarRadius_mm = max(cfg.NeedleExclusion_mm, 2*sqrt(2)*sigma_s_mm) + cfg.EdgeMargin_mm;   % blank radius + kernel half-width
Far_Mask = Bkg_Mask & (hypot(Xm, Ym) >= FarRadius_mm);

%% 2D PG on the original (fully sampled) wave field: reference estimate
[Elast_OG, R2_OG, ElastM_OG, Xc_OG, Yc_OG, ~, ~, ~, Kx_OG, Ky_OG] = ...
    PhaseGradient2D_FromSims_v3(Vid_PG, Inc_Mask, Bkg_Mask, ...
        cfg.kernel_size, cfg.Ts_OCT, T_excitation, ...
        cfg.horizontalMin, cfg.horizontalMax, cfg.verticalMin, cfg.verticalMax, ...
        cfg.applyTemporalWindow, cfg.clim);
close all force;

%% Pack
W = struct();
W.protocol    = protocol;
W.prop_type   = prop_type;
W.T_excitation = T_excitation;
W.P           = P;
W.Vid         = Vid;
W.Vid_PG      = Vid_PG;
W.x_grid      = x_grid;
W.y_grid      = y_grid;
W.t_grid      = t_grid;
W.GT_SpeedMap = GT_SpeedMap;
W.Bkg_Mask    = Bkg_Mask;
W.Inc_Mask    = Inc_Mask;
W.Ann_Mask    = Ann_Mask;
W.Far_Mask    = Far_Mask;
W.FarRadius_mm = FarRadius_mm;
W.clipInfo    = clipInfo;
W.sim_opts    = rmfield(sim_opts, 'MediumConfig');
W.MediumConfig = MediumConfig;
W.OG = struct('Elastogram', Elast_OG, 'R2_map', R2_OG, 'Elastogram_Masked', ElastM_OG, ...
    'Xaxis_mm', Xc_OG, 'Yaxis_mm', Yc_OG, 'Kx', Kx_OG, 'Ky', Ky_OG);

%% Exports of the live-script displays: wave video, S-T map, GT, OG elastogram
T_us = round(T_excitation*1e6);
switch protocol
    case "MB", tp = "Quasi_Harmonic";  tag = sprintf("%dusExc_QH_Cylindrical", T_us);
    case "SBM", tp = "Harmonic_Causal"; tag = sprintf("%dusExc_Harmonic_Causal_Cylindrical", T_us);
end
wp_path = fullfile(cfg.results_path, "Wave_Propagation", tp, sprintf("%dus", T_us));
if ~isfolder(wp_path), mkdir(wp_path); end
W.wp_path = wp_path;
x_mm = x_grid*1e3; y_mm = y_grid*1e3;
% Near-field display blanking (display copies only; W.Vid used for sampling is untouched)
R_nf = cfg.NeedleRadius_mm; Src = cfg.Source_mm;
NearXY = hypot(Xm - Src(1), Ym - Src(2)) < R_nf;                % [Nx x Ny], Vid orientation
AcqSim_WriteVideo(Vid, x_mm, y_mm, sprintf('k-Wave wave propagation - %d \\mus - FWHM %.2f mm (r < %.2f mm gray)', ...
    T_us, cfg.SourceFWHM_mm, R_nf), ...
    fullfile(wp_path, "WavePropVid_" + tag), ...
    struct('Transpose', true, 'FrameStep', 2, 'FrameRate', 30, 'dt_s', cfg.Ts_OCT, 'NaNMask', NearXY));
f = figure('Visible','off','Color','w','Position',[100 100 1250 520]);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
iy0 = round(numel(y_grid)/2);
STfull = squeeze(double(Vid(:, iy0, :))).';
ST = STfull;
ST(:, abs(x_mm - Src(1)) < R_nf & abs(y_mm(iy0) - Src(2)) < R_nf) = NaN;   % band |x| < R on the y = 0 line
lim = prctile(abs(ST(isfinite(ST))), 97);                                % far-field limits for both panels
panels = {ST, STfull};
ttls = {{sprintf('Space-time map at y = 0 mm - %d \\mus', T_us), sprintf('near field |x| < %.2f mm in gray', R_nf)}, ...
        {'Full field, same colour limits', '(near field saturated)'}};
for q = 1:2
    ax = nexttile(tl);
    im = imagesc(ax, x_mm, t_grid*1e3, panels{q}); set(im, 'AlphaData', isfinite(panels{q})); set(ax, 'Color', [0.6 0.6 0.6]);
    colormap(ax, MMB_RedBlue(255)); clim(ax, [-lim lim]);
    xlabel(ax, 'x [mm]','FontWeight','bold'); ylabel(ax, 't [ms]','FontWeight','bold');
    title(ax, ttls{q}, 'FontWeight','bold'); set(ax,'FontWeight','bold','FontSize',11);
end
cb = colorbar(ax); cb.Layout.Tile = 'east';
title(tl, char(clipInfo.type), 'FontWeight', 'bold');
savefig(f, fullfile(wp_path, "ST_Map_" + tag + ".fig")); exportgraphics(f, fullfile(wp_path, "ST_Map_" + tag + ".png"), 'Resolution', 150); close(f);
AcqSim_SaveElastogram(GT_SpeedMap.', x_mm, y_mm, cfg.clim, {'Ground truth phase speed', sprintf('c = %.2f m/s', cfg.c_ref)}, ...
    fullfile(wp_path, "GT_Map_" + tag));
AcqSim_SaveElastogram(ElastM_OG, Xc_OG, Yc_OG, cfg.clim, {sprintf('2D PG on the full wave field - %d \\mus', T_us), ...
    sprintf('far field (r \\geq %.2f mm) median %.2f m/s', FarRadius_mm, median(ElastM_OG(imresize(Far_Mask, size(ElastM_OG), 'nearest')), 'omitnan'))}, ...
    fullfile(wp_path, "2D_PG_Masked_Elastogram_WaveProp_" + tag));   % unblanked: near-field justification
% Same map with the experimental near-field blanking (W.OG stays unblanked
% for the near-field radius analysis)
[Xog, Yog] = ndgrid(Xc_OG, Yc_OG);
if ~isequal(size(Xog), size(ElastM_OG)), [Xog, Yog] = meshgrid(Xc_OG, Yc_OG); end
ElastM_OG_Blanked = ElastM_OG; ElastM_OG_Blanked(hypot(Xog - Src(1), Yog - Src(2)) < R_nf) = NaN;
AcqSim_SaveElastogram(ElastM_OG_Blanked, Xc_OG, Yc_OG, cfg.clim, {sprintf('2D PG on the full wave field - %d \\mus', T_us), ...
    sprintf('near field r < %.2f mm blanked; median %.2f m/s', R_nf, median(ElastM_OG_Blanked(:), 'omitnan'))}, ...
    fullfile(wp_path, "2D_PG_Masked_Elastogram_WaveProp_Blanked_" + tag), Src);
W.OG.Elastogram_Masked_Blanked = ElastM_OG_Blanked;
W.OG.FarField_median = median(ElastM_OG(imresize(Far_Mask, size(ElastM_OG), 'nearest')), 'omitnan');

end
