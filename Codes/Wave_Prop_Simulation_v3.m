function [Im_UZ, x_grid, y_grid, t_grid, GT_Mask] = Wave_Prop_Simulation_v3( ...
    c_medium, Ts_OCT, T_excitation, alpha, N_pulses, A_Peak, ...
    horizontalMin, horizontalMax, verticalMin, verticalMax, ...
    nPoints, N_excitations, prop_type, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 23-04-2026
% Date of last modification (DD/MM/YY): 04-05-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Description %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wave_Prop_Simulation_v3
%
% Parent function for the generation of canonical wave propagation
% simulations used in the MB-mode OCE simulation framework.
%
% This version keeps only three canonical cases:
%
%   1) 'Wave_Eq_Harmonic_Cylindrical'
%       Analytical harmonic cylindrical propagation in a homogeneous medium.
%       This is kept as a theoretical reference case.
%
%   2) 'K_Wave_Quasi_Harmonic_CustomMedium'
%       k-Wave numerical simulation using the custom medium builder and a
%       quasi-harmonic pulsed-ensemble excitation.
%
%   3) 'K_Wave_Harmonic_CustomMedium'
%       k-Wave numerical simulation using the custom medium builder and a
%       continuous harmonic excitation.
%
% IMPORTANT:
%   The custom-medium k-Wave cases can represent both cylindrical and
%   non-cylindrical media:
%
%       - cylindrical / homogeneous case:
%           use a custom medium with no inclusions.
%
%       - non-cylindrical / heterogeneous case:
%           use a custom medium with inclusions.
%
% IMPORTANT SEMANTICS:
%   - N_pulses:
%       Number of periods of the excitation waveform inside each burst.
%
%   - N_excitations:
%       Number of excitation pulse ensembles to be simulated.
%
% OPTIONAL CONTROL:
%   The function accepts an optional final input:
%
%       sim_opts (struct)
%
%   Supported parent fields:
%
%       sim_opts.visualize
%           true  -> visualize generated field
%           false -> skip visualization
%           default = true
%
%       sim_opts.KWaveCustomMedium
%           Options passed directly to the k-Wave custom-medium child
%           solvers. This should include GUI/headless controls such as:
%
%               sim_opts.KWaveCustomMedium.useGUI = false;
%               sim_opts.KWaveCustomMedium.forceDefaultConfig = true;
%               sim_opts.KWaveCustomMedium.MediumConfig = cfg;
%               sim_opts.KWaveCustomMedium.M_repetitions = M_repetitions;
%               sim_opts.KWaveCustomMedium.nSyncSamples = nSyncSamples;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%   - c_medium (scalar, double):
%       Background/reference sound speed [m/s].
%
%   - Ts_OCT (scalar, double):
%       OCT sampling period [s].
%
%   - T_excitation (scalar, double):
%       Temporal period of the excitation waveform [s].
%
%   - alpha (scalar, double):
%       Attenuation coefficient [dB/m].
%
%   - N_pulses (scalar, integer):
%       Number of periods of the excitation waveform inside each burst.
%
%   - A_Peak (scalar, double):
%       Peak amplitude of the oscillatory field.
%
%   - horizontalMin, horizontalMax (scalar, double):
%       Horizontal simulation bounds [mm].
%
%   - verticalMin, verticalMax (scalar, double):
%       Vertical simulation bounds [mm].
%
%   - nPoints (scalar, integer):
%       Number of spatial points along each dimension.
%
%   - N_excitations (scalar, integer):
%       Number of excitation pulse ensembles.
%
%   - prop_type (char array or string):
%       Supported options:
%
%           * 'Wave_Eq_Harmonic_Cylindrical'
%           * 'K_Wave_Quasi_Harmonic_CustomMedium'
%           * 'K_Wave_Harmonic_CustomMedium'
%
%   - varargin{1} (optional, struct):
%       sim_opts structure for visualization and child-solver control.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%
%   - Im_UZ (3D matrix, double):
%       Simulated field organized as:
%           [length(x_grid) x length(y_grid) x length(t_grid)]
%
%   - x_grid (vector, double):
%       Horizontal spatial axis [m].
%
%   - y_grid (vector, double):
%       Vertical spatial axis [m].
%
%   - t_grid (vector, double):
%       Temporal axis [s].
%
%   - GT_Mask (2D matrix, logical or double):
%       Ground-truth sound-speed map returned by the custom-medium k-Wave
%       solver. Empty for analytical cases.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse optional controls
sim_opts = parse_parent_optional_inputs(varargin{:});

%% Initialization
GT_Mask = [];

%% Error Checks

ratio_Texc_Ts = T_excitation / Ts_OCT;
tol = 1e-12;

if abs(ratio_Texc_Ts - round(ratio_Texc_Ts)) > tol
    error('T_excitation must be an integer multiple of Ts_OCT.');
end

if (horizontalMin >= horizontalMax) || (verticalMin >= verticalMax)
    error(['Invalid simulation dimensions: ensure horizontalMin < horizontalMax ' ...
           'and verticalMin < verticalMax.']);
end

if N_pulses <= 0 || floor(N_pulses) ~= N_pulses
    error('N_pulses must be a strictly positive integer.');
end

if N_excitations <= 0 || floor(N_excitations) ~= N_excitations
    error('N_excitations must be a strictly positive integer.');
end

if alpha < 0 || ~isfinite(alpha)
    error('alpha must be a finite non-negative scalar.');
end

%% Preliminary Calculations

f_excitation = 1 / T_excitation;

% Analytical complex wavenumber.
% Convert alpha from [dB/m] to [Np/m] for the analytical wave-equation case.
alpha_Np_per_m = alpha / (20*log10(exp(1)));

k = 2*pi*f_excitation/c_medium + 1i*alpha_Np_per_m;

%% Build Base Computational Grid for Analytical Case

x_grid = linspace(horizontalMin*1e-3, horizontalMax*1e-3, nPoints);
y_grid = linspace(verticalMin*1e-3, verticalMax*1e-3, nPoints);

t_end = N_pulses / f_excitation;
t_grid = 0:Ts_OCT:(t_end - Ts_OCT);

[X_ndgrid, Y_ndgrid, T_ndgrid] = ndgrid(x_grid, y_grid, t_grid);

%% Compute Wave Propagation Field

switch char(prop_type)

    case 'Wave_Eq_Harmonic_Cylindrical'

        R = sqrt(X_ndgrid.^2 + Y_ndgrid.^2);
        R(R == 0) = eps;

        Uz = (A_Peak ./ R) .* exp(1i * ((R .* k) - ...
            (2*pi*f_excitation*T_ndgrid)));

        Im_UZ = imag(Uz);

    case 'K_Wave_Quasi_Harmonic_CustomMedium'

        [Im_UZ, x_grid, y_grid, t_grid, GT_Mask] = ...
            CustomMedium_QuasiHarmonic_kWave_Acoustic( ...
                c_medium, Ts_OCT, T_excitation, alpha, ...
                N_pulses, A_Peak, ...
                horizontalMin, horizontalMax, verticalMin, verticalMax, ...
                nPoints, N_excitations, sim_opts);

    case 'K_Wave_Harmonic_CustomMedium'

        [Im_UZ, x_grid, y_grid, t_grid, GT_Mask] = ...
            CustomMedium_Harmonic_kWave_Acoustic( ...
                c_medium, Ts_OCT, T_excitation, alpha, ...
                N_pulses, A_Peak, ...
                horizontalMin, horizontalMax, verticalMin, verticalMax, ...
                nPoints, N_excitations, sim_opts);

    otherwise
        error(['Unsupported prop_type. Valid options are: ', ...
               '''Wave_Eq_Harmonic_Cylindrical'', ', ...
               '''K_Wave_Quasi_Harmonic_CustomMedium'', ', ...
               '''K_Wave_Harmonic_CustomMedium''.']);
end

%% Visualize the Simulated Field (optional)

if sim_opts.visualize

    % % Robust color scale based on mean ± 2 std
    % mean_Im_Uz = mean(Im_UZ(:), 'omitnan');
    % std_Im_Uz  = std(Im_UZ(:), 0, 'omitnan');
    % 
    % clim_low  = mean_Im_Uz - 2*std_Im_Uz;
    % clim_high = mean_Im_Uz + 2*std_Im_Uz;
    % clim_robust = max(abs([clim_low, clim_high]));
    clim_robust = 1.3;

    figure;

    for i = 1:length(t_grid)

        imagesc(x_grid, y_grid, Im_UZ(:,:,i));
        colormap("turbo");
        axis image;

        switch char(prop_type)

            case 'Wave_Eq_Harmonic_Cylindrical'
                title({ ...
                    'Particle Velocity Field - Harmonic Cylindrical Propagation', ...
                    'Analytical Wave Equation / Helmholtz Reference Case',''});

            case 'K_Wave_Quasi_Harmonic_CustomMedium'
                title({ ...
                    'Acoustic Pressure Field - Custom Medium', ...
                    'k-Wave Quasi-Harmonic Excitation Pulse Ensembles',''});

            case 'K_Wave_Harmonic_CustomMedium'
                title({ ...
                    'Acoustic Pressure Field - Custom Medium', ...
                    'k-Wave Harmonic Excitation',''});
        end

        xlabel('x [m]');
        ylabel('y [m]');
        axis image;

        if clim_robust > 0
            clim([-clim_robust, clim_robust]);
        end

        colorbar;
        pause(1/120);
    end
end

end

%% ===================================================================== %%
function sim_opts = parse_parent_optional_inputs(varargin)

sim_opts = struct();
sim_opts.visualize = true;

if nargin == 0 || isempty(varargin)
    return;
end

candidate = varargin{1};

if ~isstruct(candidate)
    error('Optional input to Wave_Prop_Simulation_v3 must be a struct named sim_opts.');
end

fn = fieldnames(candidate);
for k = 1:numel(fn)
    sim_opts.(fn{k}) = candidate.(fn{k});
end

end