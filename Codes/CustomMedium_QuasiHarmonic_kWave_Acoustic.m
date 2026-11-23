function [SimVideo, x_grid, y_grid, t_grid, GT_Mask] = CustomMedium_QuasiHarmonic_kWave_Acoustic( ...
                               c_medium, Ts_OCT, ...
                               T_excitation, alpha, N_pulses, A_Peak, ...
                               horizontalMin, horizontalMax, verticalMin, ...
                               verticalMax, nPoints, N_excitations, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 04-05-2026
% Date of last modification (DD/MM/YY): 04-05-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Function's Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CustomMedium_QuasiHarmonic_kWave_Acoustic
%
% This function generates quasi-harmonic wave propagation simulations using
% the k-Wave acoustic solver in a user-defined custom medium.
%
% The function supports both homogeneous and heterogeneous media. Therefore,
% it can be used to simulate:
%
%   1) Cylindrical-like propagation in a homogeneous medium
%      by defining only the background medium.
%
%   2) Non-cylindrical propagation in a heterogeneous medium
%      by adding inclusions with user-defined sound speeds.
%
% -----------------------------------------------------------------------
% QUASI-HARMONIC TEMPORAL EXCITATION
% -----------------------------------------------------------------------
%
% The excitation is not modeled as an infinite-duration harmonic signal.
% Instead, it is modeled as a finite-duration sinusoidal burst repeated in
% time, which better represents experimental wave-based OCE excitation.
%
% Each excitation ensemble contains:
%
%   - A sinusoidal carrier with period T_excitation
%   - A finite burst of N_pulses periods
%   - A silent interval before the next excitation ensemble
%
% The ensemble duration is matched to the MB-mode temporal sampling at one
% spatial location:
%
%       T_ensemble = (M_repetitions + nSyncSamples) * Ts_OCT
%
% where:
%
%   - M_repetitions:
%       Number of temporal A-lines acquired at one spatial position.
%
%   - nSyncSamples:
%       Number of synchronization/transition samples between consecutive
%       spatial positions.
%
% The full source signal is then built by repeating this ensemble
% N_excitations times.
%
% -----------------------------------------------------------------------
% CUSTOM MEDIUM BUILDER
% -----------------------------------------------------------------------
%
% The medium is defined using a custom medium builder that supports both GUI
% and headless operation.
%
% The GUI allows the user to:
%
%   - Define the background sound speed and density.
%   - Add inclusions of different geometries:
%       * circular
%       * elliptical
%       * rectangular
%   - Assign an independent sound speed to each inclusion.
%   - Preview the final sound-speed map using the "turbo" colormap.
%   - Apply optional Gaussian spatial smoothing to reduce sharp transition
%     artifacts at inclusion boundaries.
%   - Tune the smoothing strength using sigma_smooth_mm.
%   - Keep the last-used configuration in memory for repeated simulations.
%
% The headless mode allows the user to bypass the GUI and provide a
% predefined MediumConfig structure directly. This is essential for batch
% simulations, statistical repetitions, and reproducible parameter sweeps.
%
% Inclusions are painted sequentially, meaning that later inclusions
% overwrite earlier ones. This allows concentric or nested geometries such
% as keratoconus-like regions with multiple stiffness levels.
%
% The final medium used by k-Wave is stored as:
%
%       medium.sound_speed = sound_speed_map;
%       medium.density     = rho;
%
% and the same sound_speed_map is returned as:
%
%       GT_Mask
%
% -----------------------------------------------------------------------
% ATTENUATION MODEL
% -----------------------------------------------------------------------
%
% The attenuation parameter alpha is provided in [dB/m].
%
% In this implementation, attenuation is intended to reduce the amplitude of
% the propagating wave without changing the wavelength or phase velocity.
% Therefore, the k-Wave attenuation model is configured as:
%
%       medium.alpha_power = 0;
%       medium.alpha_mode  = 'no_dispersion';
%
% This makes the attenuation approximately frequency-independent and avoids
% dispersive phase effects that would bias phase-speed estimation.
%
% -----------------------------------------------------------------------
% RELATION TO MB-MODE SIMULATION FRAMEWORK
% -----------------------------------------------------------------------
%
% The generated wavefield is designed to be sampled later by the
% multimeridian MB-mode acquisition-vector sampling functions.
%
% This function therefore provides the physical wave propagation video used
% to assess phase-speed estimation error under:
%
%   - different numbers of meridians
%   - homogeneous or heterogeneous media
%   - ideal acquisition or saccadic motion artifacts
%   - quasi-harmonic excitation conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%   - c_medium (scalar, double):
%       Background sound speed [m/s].
%
%   - Ts_OCT (scalar, double):
%       OCT sampling period [s].
%
%   - T_excitation (scalar, double):
%       Period of the sinusoidal excitation waveform [s].
%
%   - alpha (scalar, double):
%       Attenuation coefficient [dB/m].
%
%   - N_pulses (scalar, integer):
%       Number of sinusoidal periods inside each excitation burst.
%
%   - A_Peak (scalar, double):
%       Peak pressure amplitude of the excitation [Pa].
%
%   - horizontalMin, horizontalMax (scalar, double):
%       Horizontal simulation bounds [mm].
%
%   - verticalMin, verticalMax (scalar, double):
%       Vertical simulation bounds [mm].
%
%   - nPoints (scalar, integer):
%       Number of grid points per spatial dimension.
%
%   - N_excitations (scalar, integer):
%       Number of quasi-harmonic excitation ensembles to simulate.
%
%   - varargin{1} (optional, struct):
%       Options structure for GUI/headless control.
%
%       Supported fields:
%
%       * useGUI (logical):
%           true  -> open medium-builder GUI
%           false -> run in headless mode
%
%       * MediumConfig (struct):
%           Predefined medium configuration used in headless mode.
%
%       * useLastConfig (logical):
%           Use the most recent saved medium configuration if available.
%
%       * saveConfig (logical):
%           Save the selected/final configuration in memory and preferences.
%
%       * forceDefaultConfig (logical):
%           Force the default keratoconus-like medium configuration.
%
%       * M_repetitions (scalar, integer):
%           Number of MB-mode temporal samples acquired at one spatial
%           position.
%
%       * nSyncSamples (scalar, integer):
%           Number of synchronization samples between consecutive spatial
%           positions.
%
%       * showSourceDiagnostics (logical):
%           If true, displays source footprint, temporal waveform, and
%           source matrix diagnostics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%
%   - SimVideo (3D matrix, double):
%       Simulated acoustic pressure field arranged as [Nx x Ny x Nt].
%
%   - x_grid (vector, double):
%       Horizontal spatial axis [m].
%
%   - y_grid (vector, double):
%       Vertical spatial axis [m].
%
%   - t_grid (vector, double):
%       Temporal axis of the cropped useful simulation video [s].
%
%   - GT_Mask (2D matrix, double):
%       Final sound-speed map used in the k-Wave simulation [m/s].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Persistent memory for GUI configuration
persistent LAST_MEDIUM_CONFIG

%% Parse optional options structure
opts = parse_optional_inputs(varargin{:});

%% Error checks
if horizontalMin >= horizontalMax || verticalMin >= verticalMax
    error(['Invalid spatial bounds: ensure horizontalMin < horizontalMax ' ...
           'and verticalMin < verticalMax.']);
end

if N_pulses <= 0 || floor(N_pulses) ~= N_pulses
    error('N_pulses must be a strictly positive integer.');
end

if N_excitations <= 0 || floor(N_excitations) ~= N_excitations
    error('N_excitations must be a strictly positive integer.');
end

if alpha < 0 || ~isfinite(alpha)
    error('alpha must be a finite non-negative scalar in [dB/m].');
end

if opts.M_repetitions <= 0 || floor(opts.M_repetitions) ~= opts.M_repetitions
    error('opts.M_repetitions must be a strictly positive integer.');
end

if opts.nSyncSamples < 0 || floor(opts.nSyncSamples) ~= opts.nSyncSamples
    error('opts.nSyncSamples must be a non-negative integer.');
end

ratio_Texc_Ts = T_excitation / Ts_OCT;
if abs(ratio_Texc_Ts - round(ratio_Texc_Ts)) > 1e-12
    error('T_excitation must be an integer multiple of Ts_OCT.');
end

%% Initialize the computational grid
Nx = nPoints;
Ny = nPoints;

Lx_m = (horizontalMax - horizontalMin) * 1e-3;
Ly_m = (verticalMax - verticalMin) * 1e-3;

dx = Lx_m / nPoints;
dy = Ly_m / nPoints;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Physical axes [m], consistent with user-defined rectangular domain
x_grid = linspace(horizontalMin*1e-3, horizontalMax*1e-3, Nx);
y_grid = linspace(verticalMin*1e-3, verticalMax*1e-3, Ny);

%% Time grid
dt = Ts_OCT;

% One ensemble corresponds to one MB-mode temporal record at a single
% spatial position:
samples_per_ensemble = opts.M_repetitions + opts.nSyncSamples;
T_ensemble = samples_per_ensemble * Ts_OCT;

% Excitation burst duration inside each ensemble
T_burst = N_pulses * T_excitation;
burst_samples = round(T_burst / dt);

if burst_samples > samples_per_ensemble
    error(['The excitation burst does not fit inside one MB-mode ensemble. ' ...
           'Increase M_repetitions+nSyncSamples or reduce N_pulses.']);
end

% Time required for the excitation to appear in the medium
T_transient = 0.185e-3;% + T_excitation;

% Total simulation duration
T_end = T_transient + N_excitations*T_ensemble - Ts_OCT;
Nt = round(T_end / dt) + 1;

kgrid.setTime(Nt, dt);

%% Medium creation
default_cfg = build_default_medium_config( ...
    c_medium, horizontalMin, horizontalMax, verticalMin, verticalMax);

if opts.useGUI
    % Select starting config for GUI
    if opts.forceDefaultConfig
        cfg = default_cfg;
    else
        if opts.useLastConfig
            if ~isempty(LAST_MEDIUM_CONFIG)
                cfg = LAST_MEDIUM_CONFIG;
            else
                cfg = recover_saved_config(default_cfg);
            end
        else
            cfg = default_cfg;
        end
    end

    % Adapt config to current domain
    cfg = adapt_config_to_current_domain(cfg, default_cfg, ...
        horizontalMin, horizontalMax, verticalMin, verticalMax, c_medium);

    [cfg, sound_speed_map, rho] = launch_medium_builder_gui( ...
        cfg, x_grid, y_grid, dx, dy, ...
        horizontalMin, horizontalMax, verticalMin, verticalMax);

    if opts.saveConfig
        LAST_MEDIUM_CONFIG = cfg;
        save_config_to_preferences(cfg);
    end

else
    % Headless mode
    if isempty(opts.MediumConfig)
        if opts.forceDefaultConfig
            cfg = default_cfg;
        else
            if opts.useLastConfig
                if ~isempty(LAST_MEDIUM_CONFIG)
                    cfg = LAST_MEDIUM_CONFIG;
                else
                    cfg = recover_saved_config(default_cfg);
                end
            else
                cfg = default_cfg;
            end
        end
    else
        cfg = opts.MediumConfig;
    end

    cfg = adapt_config_to_current_domain(cfg, default_cfg, ...
        horizontalMin, horizontalMax, verticalMin, verticalMax, c_medium);

    rho = cfg.density;
    [~, sound_speed_map] = build_medium_from_config(cfg, x_grid, y_grid, dx, dy);

    if opts.saveConfig
        LAST_MEDIUM_CONFIG = cfg;
        save_config_to_preferences(cfg);
    end
end

% Ground truth map actually used in the simulation
GT_Mask = sound_speed_map;

% Assign medium
medium.sound_speed = sound_speed_map;
medium.density     = rho;

% Assign acoustic attenuation if requested.
% alpha is given in [dB/m].
%
% Goal:
%   - attenuate propagating amplitude over distance
%   - avoid frequency-dependent attenuation
%   - avoid dispersion / wavelength alteration
%
% k-Wave expects alpha_coeff in [dB / (MHz^y cm)].
% With alpha_power = 0, this becomes simply [dB/cm].

if alpha > 0
    medium.alpha_power = 0;              % frequency-independent attenuation
    medium.alpha_coeff = alpha / 100;    % [dB/m] -> [dB/cm]
    medium.alpha_mode  = 'no_dispersion';

    fprintf('\n===== k-Wave Attenuation Enabled =====\n');
    fprintf('alpha input        = %.4f dB/m\n', alpha);
    fprintf('alpha_coeff        = %.6f dB/cm\n', medium.alpha_coeff);
    fprintf('alpha_power        = %.2f\n', medium.alpha_power);
    fprintf('alpha_mode         = %s\n', medium.alpha_mode);

    R_domain_m = 0.5 * sqrt(Lx_m^2 + Ly_m^2);
    fprintf('domain radius      = %.3f mm\n', R_domain_m*1e3);
    fprintf('attenuation radius = %.4f dB\n', alpha * R_domain_m);
    fprintf('======================================\n\n');
else
    fprintf('k-Wave attenuation disabled: alpha = 0 dB/m\n');
end

%% Excitation source and sensor definition

source_freq = 1 / T_excitation;   % [Hz]
P0 = A_Peak;                      % [Pa]

% Build one MB-mode quasi-harmonic ensemble:
% - its total duration matches one spatial-position temporal record
% - the excitation is active only during the first N_pulses periods
t_ensemble_local = (0:samples_per_ensemble-1) * dt;
carrier_square_OOK = (t_ensemble_local < N_pulses*T_excitation);
% carrier_square_OOK = ones(1,length(t_ensemble_local));
ensemble_sig = (P0 * sin(2*pi*source_freq*t_ensemble_local)) .* carrier_square_OOK;

% Build quasi-harmonic multi-ensemble source signal
time_sig = zeros(1, Nt);

for kExc = 1:N_excitations
    idx_start = 1 + (kExc - 1) * samples_per_ensemble;
    idx_end   = idx_start + samples_per_ensemble - 1;

    if idx_start > Nt
        break;
    end

    idx_end = min(idx_end, Nt);
    n_local = idx_end - idx_start + 1;

    time_sig(idx_start:idx_end) = ensemble_sig(1:n_local);
end

% Keep within grid temporal bandwidth
time_sig = filterTimeSeries(kgrid, medium, time_sig);

%% Define spatial waveform (Gaussian footprint)

[X_src, Y_src] = ndgrid(x_grid, y_grid);

% Excitation center at (0,0) if contained inside the domain, otherwise use
% the geometric center of the user-defined rectangle.
if (horizontalMin <= 0) && (horizontalMax >= 0)
    x0 = 0;
else
    x0 = mean([horizontalMin, horizontalMax]) * 1e-3;
end

if (verticalMin <= 0) && (verticalMax >= 0)
    y0 = 0;
else
    y0 = mean([verticalMin, verticalMax]) * 1e-3;
end

FWHM_xy  = opts.SourceFWHM_mm * 1e-3;    % [m] (default 0.3 mm, thesis value)
sigma_xy = FWHM_xy/(2*sqrt(2*log(2))); % [m]

Gxy = exp(-((X_src - x0).^2 + (Y_src - y0).^2) / (2*sigma_xy^2));
Gxy = Gxy ./ (max(Gxy(:)) + eps);

% Compact support
support = Gxy > exp(-4);

%% Build k-Wave acoustic source

source = struct();
source.p_mask = false(Nx, Ny);
source.p_mask(support) = true;

weights = single(Gxy(support));
source.p = single(weights) * single(time_sig(:)).';

%% Source Diagnostics Figure
% Optional visualization of:
%   1) Spatial Gaussian source footprint
%   2) Temporal excitation waveform
%   3) Full source matrix source.p = Gxy(support) * time_sig(t)

if opts.showSourceDiagnostics

    % Reconstruct full-grid source footprint for visualization
    Gxy_vis = zeros(Nx, Ny);
    Gxy_vis(support) = Gxy(support);

    % Full source matrix: [Nsrc x Nt]
    SourceMatrix_vis = double(source.p);

    % Source spatial coordinates for compact display
    x_src_mm = x_grid * 1e3;
    y_src_mm = y_grid * 1e3;

    % Temporal axis
    t_ms = kgrid.t_array * 1e3;

    figure('Name', 'k-Wave Source Diagnostics', ...
           'Units', 'normalized', ...
           'Position', [0.05 0.18 0.90 0.48]);

    % --------------------------------------------------------------- %
    % 1) Spatial Gaussian Footprint
    % --------------------------------------------------------------- %
    subplot(1,3,1);
    imagesc(x_src_mm, y_src_mm, Gxy_vis);
    set(gca, 'YDir', 'reverse');
    axis image;
    colormap(gca, turbo);
    colorbar;
    xlabel('x (mm)');
    ylabel('y (mm)');
    title('Spatial Gaussian Footprint G_{xy}');

    % --------------------------------------------------------------- %
    % 2) Temporal waveform
    % --------------------------------------------------------------- %
    subplot(1,3,2);
    plot(t_ms, double(time_sig), 'r', 'LineWidth', 1.5);
    grid on;
    xlabel('time (ms)');
    ylabel('amplitude');
    title('Temporal Waveform time\_sig');
    xlim([t_ms(1), t_ms(end)]);

    % --------------------------------------------------------------- %
    % 3) Full source matrix
    % --------------------------------------------------------------- %
    subplot(1,3,3);
    imagesc(t_ms, 1:size(SourceMatrix_vis,1), SourceMatrix_vis);
    axis tight;
    colormap(gca, turbo);
    colorbar;
    xlabel('time (ms)');
    ylabel('source point index');
    title('Full Source Matrix s_z(Pa) [N_{src} x N_t]');

end

%% Sensor mask

sensor = struct();
sensor.mask   = true(Nx, Ny);
sensor.record = {'p', 'u'};

%% Run the 2D acoustic solver

input_args = {'PMLInside', true, 'PMLSize', [20 20], 'DataCast', 'single'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

%% Reorder to full (x,y,t) cubes

p_xy_t  = kwave_reshape_fullgrid(kgrid, sensor.mask, sensor_data.p);
ux_xy_t = kwave_reshape_fullgrid(kgrid, sensor.mask, sensor_data.ux);
uy_xy_t = kwave_reshape_fullgrid(kgrid, sensor.mask, sensor_data.uy);

%% Crop video

% Start the useful movie when the excitation is expected to appear in the
% medium, so that the resulting video is directly compatible with the
% MB-mode sampling function.
idx_start = round(T_transient / Ts_OCT) + 1;
idx_start = min(max(idx_start, 1), size(p_xy_t, 3));

% Keep exactly the duration corresponding to the requested number of
% excitation ensembles.
T_visible = N_excitations * T_ensemble;
idx_end = idx_start - 1 + round(T_visible / Ts_OCT);
idx_end = min(size(p_xy_t,3), max(idx_start, idx_end));

p_xy_t_h  = p_xy_t(:,:,idx_start:idx_end); %#ok<NASGU>
ux_xy_t_h = ux_xy_t(:,:,idx_start:idx_end); %#ok<NASGU>
uy_xy_t_h = uy_xy_t(:,:,idx_start:idx_end); %#ok<NASGU>

t_h = kgrid.t_array(idx_start:idx_end);
t_h = t_h - t_h(1);

%% Choose output video
SimVideo = p_xy_t(:,:,idx_start:idx_end);

%% Export time axis
t_grid = t_h;
%% Choose output video
SimVideo = p_xy_t(:,:,idx_start:idx_end);

%% Export time axis
t_grid = t_h;

end

%% ===================================================================== %%
function opts = parse_optional_inputs(varargin)

opts = struct();
opts.useGUI = true;
opts.MediumConfig = [];
opts.useLastConfig = true;
opts.saveConfig = true;
opts.forceDefaultConfig = false;
opts.M_repetitions = 200;
opts.nSyncSamples = 50;
opts.showSourceDiagnostics = false;
opts.SourceFWHM_mm = 0.3; % Gaussian source FWHM [mm] (added 20-10-2026 for needle-matched sims)

if nargin == 0 || isempty(varargin)
    return;
end

candidate = varargin{1};

if ~isstruct(candidate)
    error('Optional input must be a struct named opts.');
end

fn = fieldnames(candidate);
for k = 1:numel(fn)
    opts.(fn{k}) = candidate.(fn{k});
end

end

%% ===================================================================== %%
function cfg = build_default_medium_config(c_medium, horizontalMin, horizontalMax, verticalMin, verticalMax)

% Domain center in mm
cx0 = mean([horizontalMin, horizontalMax]);

% Keratoconic region shifted by -2 mm in Y direction
cy_kc = -2.0;   % [mm]

% Background properties
rho_default = 1000;   % [kg/m^3]
E_0 = 3 * rho_default * c_medium^2;

% Inclusion sound speeds from previous keratoconic case
Cs_Background = c_medium;
Cs_outer_ring = sqrt((0.8 * E_0) / (3 * rho_default));
Cs_inner_ring = sqrt((0.6 * E_0) / (3 * rho_default));
Cs_inner_disc = sqrt((0.4 * E_0) / (3 * rho_default));

cfg = struct();
cfg.background_speed = Cs_Background;
cfg.density = rho_default;
cfg.apply_smoothing = true;
cfg.sigma_smooth_mm = 0.30;

cfg.domain.horizontalMin = horizontalMin;
cfg.domain.horizontalMax = horizontalMax;
cfg.domain.verticalMin   = verticalMin;
cfg.domain.verticalMax   = verticalMax;

% Default three-inclusion keratoconic medium:
%   - outer ring: radius 2.4 mm (diameter 4.8 mm)
%   - inner ring: radius 1.5 mm (diameter 3.0 mm)
%   - inner disc: radius 1.0 mm (diameter 2.0 mm)
%
% Inclusions are painted in order, so later ones overwrite earlier ones.
cfg.inclusions = struct( ...
    'type',        {'circle',       'circle',       'circle'}, ...
    'center_x_mm', {cx0,            cx0,            cx0}, ...
    'center_y_mm', {cy_kc,          cy_kc,          cy_kc}, ...
    'size_x_mm',   {2.4,            1.5,            1.0}, ...
    'size_y_mm',   {2.4,            1.5,            1.0}, ...
    'angle_deg',   {0,              0,              0}, ...
    'speed',       {Cs_outer_ring,  Cs_inner_ring,  Cs_inner_disc});
end

%% ===================================================================== %%
function cfg = recover_saved_config(default_cfg)

cfg = default_cfg;

try
    if ispref('OCE_MB_Sim', 'LastMediumConfig')
        cfg = getpref('OCE_MB_Sim', 'LastMediumConfig');
    end
catch
    cfg = default_cfg;
end

end

%% ===================================================================== %%
function save_config_to_preferences(cfg)

try
    setpref('OCE_MB_Sim', 'LastMediumConfig', cfg);
catch
    % If preferences cannot be written, continue silently.
end

end

%% ===================================================================== %%
function cfg = adapt_config_to_current_domain(cfg, default_cfg, ...
    horizontalMin, horizontalMax, verticalMin, verticalMax, c_medium)

if ~isstruct(cfg)
    cfg = default_cfg;
    return;
end

required_fields = {'background_speed','density','apply_smoothing','sigma_smooth_mm','inclusions'};
for i = 1:numel(required_fields)
    if ~isfield(cfg, required_fields{i})
        cfg = default_cfg;
        return;
    end
end

cfg.domain.horizontalMin = horizontalMin;
cfg.domain.horizontalMax = horizontalMax;
cfg.domain.verticalMin   = verticalMin;
cfg.domain.verticalMax   = verticalMax;

if ~isfield(cfg, 'background_speed') || isempty(cfg.background_speed) || ~isfinite(cfg.background_speed)
    cfg.background_speed = c_medium;
end

if ~isfield(cfg, 'density') || isempty(cfg.density) || ~isfinite(cfg.density)
    cfg.density = 1000;
end

if isempty(cfg.inclusions)
    cfg.inclusions = default_cfg.inclusions([]);
end

for k = 1:numel(cfg.inclusions)
    if ~isfield(cfg.inclusions(k), 'type');        cfg.inclusions(k).type = 'circle'; end
    if ~isfield(cfg.inclusions(k), 'center_x_mm'); cfg.inclusions(k).center_x_mm = mean([horizontalMin, horizontalMax]); end
    if ~isfield(cfg.inclusions(k), 'center_y_mm'); cfg.inclusions(k).center_y_mm = mean([verticalMin, verticalMax]); end
    if ~isfield(cfg.inclusions(k), 'size_x_mm');   cfg.inclusions(k).size_x_mm = 1.0; end
    if ~isfield(cfg.inclusions(k), 'size_y_mm');   cfg.inclusions(k).size_y_mm = 1.0; end
    if ~isfield(cfg.inclusions(k), 'angle_deg');   cfg.inclusions(k).angle_deg = 0; end
    if ~isfield(cfg.inclusions(k), 'speed');       cfg.inclusions(k).speed = c_medium; end

    cfg.inclusions(k).center_x_mm = min(max(cfg.inclusions(k).center_x_mm, horizontalMin), horizontalMax);
    cfg.inclusions(k).center_y_mm = min(max(cfg.inclusions(k).center_y_mm, verticalMin), verticalMax);
    cfg.inclusions(k).size_x_mm = max(cfg.inclusions(k).size_x_mm, eps);
    cfg.inclusions(k).size_y_mm = max(cfg.inclusions(k).size_y_mm, eps);
end

end

%% ===================================================================== %%
function [cfg, sound_speed_map, rho] = launch_medium_builder_gui( ...
    cfg, x_grid, y_grid, dx, dy, ...
    horizontalMin, horizontalMax, verticalMin, verticalMax)

rho = cfg.density;

currentSelection = [];
sound_speed_map = [];
ui_finished = false;
ui_cancelled = false;

fig = figure( ...
    'Name', 'Medium Builder - Non Cylindrical QH k-Wave v3', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Color', [0.94 0.94 0.94], ...
    'Units', 'normalized', ...
    'Position', [0.08 0.08 0.84 0.80], ...
    'CloseRequestFcn', @on_close_request);

% Axes for preview
ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.42 0.12 0.55 0.80]);

% --------------------------------------------------------------------- %
% Background / global controls
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.93 0.12 0.03], 'String', 'Background c [m/s]', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

edit_bg = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.93 0.10 0.035], ...
    'String', num2str(cfg.background_speed), ...
    'BackgroundColor', 'w', 'Callback', @refresh_preview_callback);

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.885 0.12 0.03], 'String', 'Density [kg/m^3]', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

edit_rho = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.885 0.10 0.035], ...
    'String', num2str(cfg.density), ...
    'BackgroundColor', 'w', 'Callback', @refresh_preview_callback);

check_smooth = uicontrol(fig, 'Style', 'checkbox', 'Units', 'normalized', ...
    'Position', [0.02 0.835 0.18 0.035], 'String', 'Apply smoothing', ...
    'Value', logical(cfg.apply_smoothing), ...
    'Callback', @refresh_preview_callback);

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.795 0.12 0.03], 'String', '\sigma smooth [mm]', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

edit_sigma = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.795 0.10 0.035], ...
    'String', num2str(cfg.sigma_smooth_mm), ...
    'BackgroundColor', 'w', 'Callback', @refresh_preview_callback);

% --------------------------------------------------------------------- %
% Inclusion editor
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.735 0.18 0.03], 'String', 'Inclusion Type', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

popup_type = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', ...
    'Position', [0.15 0.735 0.16 0.04], ...
    'String', {'circle', 'ellipse', 'rectangle'}, ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.69 0.12 0.03], 'String', 'Center y [mm]', ...
    'HorizontalAlignment', 'left');

edit_cx = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.69 0.10 0.035], ...
    'String', num2str(mean([horizontalMin horizontalMax])), ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.65 0.12 0.03], 'String', 'Center x [mm]', ...
    'HorizontalAlignment', 'left');

edit_cy = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.65 0.10 0.035], ...
    'String', num2str(mean([verticalMin verticalMax])), ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.61 0.18 0.03], 'String', 'Size y [mm]', ...
    'HorizontalAlignment', 'left');

edit_sx = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.61 0.10 0.035], ...
    'String', '1.0', ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.57 0.18 0.03], 'String', 'Size x / radius [mm]', ...
    'HorizontalAlignment', 'left');

edit_sy = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.57 0.10 0.035], ...
    'String', '1.0', ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.53 0.12 0.03], 'String', 'Angle [deg]', ...
    'HorizontalAlignment', 'left');

edit_ang = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.53 0.10 0.035], ...
    'String', '0', ...
    'BackgroundColor', 'w');

uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.49 0.12 0.03], 'String', 'Speed [m/s]', ...
    'HorizontalAlignment', 'left');

edit_spd = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.15 0.49 0.10 0.035], ...
    'String', num2str(0.8*cfg.background_speed), ...
    'BackgroundColor', 'w');

btn_add = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.02 0.43 0.09 0.045], 'String', 'Add', ...
    'Callback', @on_add_inclusion);

btn_update = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.12 0.43 0.09 0.045], 'String', 'Update', ...
    'Callback', @on_update_inclusion);

btn_remove = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.22 0.43 0.09 0.045], 'String', 'Remove', ...
    'Callback', @on_remove_inclusion);

list_inclusions = uicontrol(fig, 'Style', 'listbox', 'Units', 'normalized', ...
    'Position', [0.02 0.18 0.29 0.22], ...
    'BackgroundColor', 'w', ...
    'Callback', @on_select_inclusion);

btn_preview = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.02 0.12 0.14 0.045], 'String', 'Refresh Preview', ...
    'Callback', @refresh_preview_callback);

btn_done = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.17 0.12 0.14 0.045], 'String', 'Use Configuration', ...
    'FontWeight', 'bold', 'Callback', @on_done);

btn_reset = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.02 0.06 0.14 0.045], 'String', 'Reset Defaults', ...
    'Callback', @on_reset_defaults);

btn_cancel = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.17 0.06 0.14 0.045], 'String', 'Cancel', ...
    'Callback', @on_cancel);

% Initialize list and preview
refresh_inclusion_list();
refresh_preview();

uiwait(fig);

if ui_cancelled
    if isvalid(fig)
        delete(fig);
    end
    error('Medium builder GUI cancelled by user.');
end

if ~ui_finished
    if isvalid(fig)
        delete(fig);
    end
    error('Medium builder GUI closed before confirming the configuration.');
end

if isvalid(fig)
    delete(fig);
end

rho = cfg.density;
[~, sound_speed_map] = build_medium_from_config(cfg, x_grid, y_grid, dx, dy);

% ========================= nested callbacks ============================ %
    function on_add_inclusion(~, ~)
        incl = read_inclusion_editor();
        cfg.inclusions = [cfg.inclusions, incl];
        currentSelection = numel(cfg.inclusions);
        refresh_inclusion_list();
        list_inclusions.Value = currentSelection;
        refresh_preview();
    end

    function on_update_inclusion(~, ~)
        if isempty(currentSelection) || currentSelection < 1 || currentSelection > numel(cfg.inclusions)
            warndlg('Select an inclusion first.', 'Update Inclusion');
            return;
        end
        cfg.inclusions(currentSelection) = read_inclusion_editor();
        refresh_inclusion_list();
        list_inclusions.Value = currentSelection;
        refresh_preview();
    end

    function on_remove_inclusion(~, ~)
        if isempty(currentSelection) || currentSelection < 1 || currentSelection > numel(cfg.inclusions)
            warndlg('Select an inclusion first.', 'Remove Inclusion');
            return;
        end
        cfg.inclusions(currentSelection) = [];
        if isempty(cfg.inclusions)
            currentSelection = [];
        else
            currentSelection = max(1, min(currentSelection, numel(cfg.inclusions)));
        end
        refresh_inclusion_list();
        if ~isempty(currentSelection)
            list_inclusions.Value = currentSelection;
        end
        refresh_preview();
    end

    function on_select_inclusion(~, ~)
        if isempty(cfg.inclusions)
            currentSelection = [];
            return;
        end
        val = list_inclusions.Value;
        if isempty(val) || val < 1 || val > numel(cfg.inclusions)
            return;
        end
        currentSelection = val;
        load_inclusion_to_editor(cfg.inclusions(val));
    end

    function on_done(~, ~)
        if ~validate_global_inputs()
            return;
        end

        % Capture the latest global controls before closing
        cfg.background_speed = str2double(edit_bg.String);
        cfg.density = str2double(edit_rho.String);
        cfg.apply_smoothing = logical(check_smooth.Value);
        cfg.sigma_smooth_mm = str2double(edit_sigma.String);

        ui_finished = true;
        uiresume(fig);
    end

    function on_cancel(~, ~)
        ui_cancelled = true;
        uiresume(fig);
    end

    function on_close_request(~, ~)
        ui_cancelled = true;
        uiresume(fig);
    end

    function on_reset_defaults(~, ~)
        cfg = build_default_medium_config( ...
            str2double_safe(edit_bg.String, cfg.background_speed), ...
            horizontalMin, horizontalMax, verticalMin, verticalMax);

        edit_bg.String = num2str(cfg.background_speed);
        edit_rho.String = num2str(cfg.density);
        check_smooth.Value = cfg.apply_smoothing;
        edit_sigma.String = num2str(cfg.sigma_smooth_mm);

        currentSelection = [];
        refresh_inclusion_list();
        refresh_preview();
    end

    function refresh_preview_callback(~, ~)
        refresh_preview();
    end

    function refresh_preview()
        if ~validate_global_inputs(false)
            cla(ax);
            return;
        end

        cfg.background_speed = str2double(edit_bg.String);
        cfg.density = str2double(edit_rho.String);
        cfg.apply_smoothing = logical(check_smooth.Value);
        cfg.sigma_smooth_mm = str2double(edit_sigma.String);

        [~, map_preview] = build_medium_from_config(cfg, x_grid, y_grid, dx, dy);

        imagesc(ax, x_grid*1e3, y_grid*1e3, map_preview);
        set(ax, 'YDir', 'reverse');
        axis(ax, 'image');
        colormap(ax, turbo);
        colorbar(ax);
        xlabel(ax, 'x [mm]');
        ylabel(ax, 'y [mm]');
        title(ax, 'Medium Preview - Sound Speed Map [m/s]');
        drawnow;
    end

    function refresh_inclusion_list()
        if isempty(cfg.inclusions)
            list_inclusions.String = {'<no inclusions>'};
            list_inclusions.Value = 1;
            return;
        end

        strs = cell(1, numel(cfg.inclusions));
        for kk = 1:numel(cfg.inclusions)
            inc = cfg.inclusions(kk);
            strs{kk} = sprintf('%d) %s | c=%.3g | (%.2f, %.2f) mm', ...
                kk, inc.type, inc.speed, inc.center_y_mm, inc.center_x_mm);
        end
        list_inclusions.String = strs;
        if isempty(currentSelection)
            list_inclusions.Value = 1;
        else
            list_inclusions.Value = min(max(currentSelection,1),numel(cfg.inclusions));
        end
    end

    function ok = validate_global_inputs(showDialog)
        if nargin < 1
            showDialog = true;
        end

        ok = true;

        bg = str2double(edit_bg.String);
        rr = str2double(edit_rho.String);
        sg = str2double(edit_sigma.String);

        if ~isfinite(bg) || bg <= 0
            ok = false;
            if showDialog
                errordlg('Background speed must be a positive scalar.', 'Invalid Input');
            end
            return;
        end

        if ~isfinite(rr) || rr <= 0
            ok = false;
            if showDialog
                errordlg('Density must be a positive scalar.', 'Invalid Input');
            end
            return;
        end

        if ~isfinite(sg) || sg < 0
            ok = false;
            if showDialog
                errordlg('Smoothing sigma must be a non-negative scalar.', 'Invalid Input');
            end
            return;
        end
    end

    function incl = read_inclusion_editor()
        typeList = popup_type.String;
        typeVal  = popup_type.Value;

        incl = struct();
        incl.type        = typeList{typeVal};
        incl.center_x_mm = str2double_safe(edit_cx.String, mean([horizontalMin horizontalMax]));
        incl.center_y_mm = str2double_safe(edit_cy.String, mean([verticalMin verticalMax]));
        incl.size_x_mm   = abs(str2double_safe(edit_sx.String, 1.0));
        incl.size_y_mm   = abs(str2double_safe(edit_sy.String, 1.0));
        incl.angle_deg   = str2double_safe(edit_ang.String, 0);
        incl.speed       = str2double_safe(edit_spd.String, cfg.background_speed);

        if strcmpi(incl.type, 'circle')
            incl.size_y_mm = incl.size_x_mm;
        end
    end

    function load_inclusion_to_editor(incl)
        switch lower(incl.type)
            case 'circle'
                popup_type.Value = 1;
            case 'ellipse'
                popup_type.Value = 2;
            otherwise
                popup_type.Value = 3;
        end

        edit_cx.String  = num2str(incl.center_x_mm);
        edit_cy.String  = num2str(incl.center_y_mm);
        edit_sx.String  = num2str(incl.size_x_mm);
        edit_sy.String  = num2str(incl.size_y_mm);
        edit_ang.String = num2str(incl.angle_deg);
        edit_spd.String = num2str(incl.speed);
    end

end

%% ===================================================================== %%
function [E_map_smooth, sound_speed_map] = build_medium_from_config(cfg, x_grid, y_grid, dx, dy)

rho = cfg.density;
c_bg = cfg.background_speed;

E0 = 3 * rho * c_bg^2;
E_map = E0 * ones(numel(x_grid), numel(y_grid));

[X, Y] = ndgrid(x_grid*1e3, y_grid*1e3); % use mm coordinates

% Paint inclusions in order; later inclusions overwrite previous ones
for k = 1:numel(cfg.inclusions)
    inc = cfg.inclusions(k);
    mask = build_inclusion_mask(X, Y, inc);

    c_region = inc.speed;
    E_region = 3 * rho * c_region^2;
    E_map(mask) = E_region;
end

if cfg.apply_smoothing && cfg.sigma_smooth_mm > 0
    sigma_smooth_x = cfg.sigma_smooth_mm / (dx*1e3);
    sigma_smooth_y = cfg.sigma_smooth_mm / (dy*1e3);

    kernel_half_size_x = max(1, ceil(3*sigma_smooth_x));
    kernel_half_size_y = max(1, ceil(3*sigma_smooth_y));

    kernel_size_x = 2*kernel_half_size_x + 1;
    kernel_size_y = 2*kernel_half_size_y + 1;

    G_kernel = gaussian2D_kernel(kernel_size_x, kernel_size_y, sigma_smooth_x, sigma_smooth_y);
    E_map_smooth = conv2_symmetric(E_map, G_kernel);
else
    E_map_smooth = E_map;
end

sound_speed_map = sqrt(E_map_smooth ./ (3*rho));

end

%% ===================================================================== %%
function mask = build_inclusion_mask(X_mm, Y_mm, inc)

cx = inc.center_x_mm;
cy = inc.center_y_mm;
sx = max(inc.size_x_mm, eps);
sy = max(inc.size_y_mm, eps);
th = deg2rad(inc.angle_deg);

% Shift and rotate coordinates
Xc = X_mm - cx;
Yc = Y_mm - cy;

Xr =  cos(th)*Xc + sin(th)*Yc;
Yr = -sin(th)*Xc + cos(th)*Yc;

switch lower(inc.type)
    case 'circle'
        r = sx;
        mask = (Xr.^2 + Yr.^2) <= r^2;

    case 'ellipse'
        ax = sx;
        by = sy;
        mask = (Xr./ax).^2 + (Yr./by).^2 <= 1;

    case 'rectangle'
        hx = sx/2;
        hy = sy/2;
        mask = abs(Xr) <= hx & abs(Yr) <= hy;

    otherwise
        error('Unsupported inclusion type: %s', inc.type);
end

end

%% ===================================================================== %%
function val = str2double_safe(str, fallback)
val = str2double(str);
if ~isfinite(val)
    val = fallback;
end
end

%% ===================================================================== %%
function vol_xy_t = kwave_reshape_fullgrid(kgrid, mask, data)
% data : [M x Nt] time series (M = number of active mask points)
% mask : [Nx x Ny] logical full-grid mask
% out  : [Nx x Ny x Nt]

    Nx = kgrid.Nx;
    Ny = kgrid.Ny;
    Nt = size(data,2);

    lin_idx = find(mask(:));
    tmp = zeros(Nx*Ny, Nt, 'like', data);
    tmp(lin_idx, :) = data;

    vol_xy_t = reshape(tmp, Nx, Ny, Nt);
end

%% ===================================================================== %%
function G = gaussian2D_kernel(nx, ny, sigma_x, sigma_y)
% Build a normalized 2D Gaussian kernel

    cx = (nx + 1)/2;
    cy = (ny + 1)/2;

    [X, Y] = ndgrid(1:nx, 1:ny);

    G = exp( -((X - cx).^2)/(2*sigma_x^2) - ((Y - cy).^2)/(2*sigma_y^2) );
    G = G ./ sum(G(:));
end

%% ===================================================================== %%
function out = conv2_symmetric(in, kernel)
% 2D convolution with symmetric padding to avoid edge artifacts.
% This avoids requiring the Image Processing Toolbox.

    [kx, ky] = size(kernel);
    pad_x = floor(kx/2);
    pad_y = floor(ky/2);

    in_pad = padarray_symmetric_local(in, [pad_x, pad_y]);
    out = conv2(in_pad, kernel, 'valid');
end

%% ===================================================================== %%
function A_pad = padarray_symmetric_local(A, pad_sz)
% Symmetric padding without relying on padarray()

    px = pad_sz(1);
    py = pad_sz(2);

    % Pad in x direction
    top    = A(px:-1:1, :);
    bottom = A(end:-1:end-px+1, :);
    A1 = [top; A; bottom];

    % Pad in y direction
    left  = A1(:, py:-1:1);
    right = A1(:, end:-1:end-py+1);
    A_pad = [left, A1, right];
end