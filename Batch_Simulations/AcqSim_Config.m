function cfg = AcqSim_Config()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcqSim_Config
%
% Single source of truth for the experiment-matched simulation batches
% (Multimeridian MB mode and SBM). Everything that is not listed here is
% taken unchanged from the thesis framework in ../Codes.
%
% Physical parameters
%   - Homogeneous medium (cylindrical propagation, no inclusions) with the
%     Rayleigh speed of the PC reference phantom: c_ref = 10.54 m/s
%     (quasi-static compression E0 = 0.365 MPa, rho = 1000 kg/m^3).
%   - 18 x 18 mm domain (new acquisitions, 20-10-2026).
%   - Gaussian source FWHM = new needle width measured on the crosshair
%     B-mode of 20-10-2026 (10 x 10 mm scan, 300 A-scans per B-scan):
%     shadow edges at A-scans 89.5 and 195.8 -> 106.3 A-scans -> 3.54 mm.
%
% Run lists
%   - MB : Experiment_Log_06_10_2026.xlsx (sheet PH_1150_503_06_10_2026),
%          one noise-free repetition per (period, meridians) pair.
%   - SBM: Experiment_Log_SBM_19_10_2026.xlsx (sheet PH_1150_503_SBM),
%          N_A = 37, N_P = 4, N_F = P, N_phi = 4, N_turn = smallest
%          k*P - N_A >= 20, one noise-free repetition per row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = struct();

%% Paths
cfg.root_path    = string(fileparts(fileparts(mfilename("fullpath"))));   % repository root
cfg.codes_path   = cfg.root_path + "/Codes";
cfg.kwave_path   = cfg.codes_path + "/k-Wave";   % download k-Wave (http://www.k-wave.org) and place it here
cfg.batch_path   = cfg.root_path + "/Batch_Simulations";
cfg.results_path = cfg.root_path + "/Results/Experiment_Matched_18x18mm";
cfg.display_path = cfg.batch_path + "/display_helpers";   % colormaps / plot helpers

%% Medium / source
cfg.c_ref        = 10.54;    % [m/s] PC phantom reference (Rayleigh)
cfg.rho          = 1000;     % [kg/m^3]
cfg.alpha        = 0;        % [dB/m] no attenuation (thesis cylindrical case)
cfg.A_Peak       = 1;
cfg.SourceFWHM_mm = 3.54;    % [mm] new needle width (see header)

%% Grid / clock
cfg.Ts_OCT  = 5e-6;          % [s] 200 kHz A-line clock
cfg.FOV_mm  = 18;
cfg.horizontalMin = -cfg.FOV_mm/2; cfg.horizontalMax = cfg.FOV_mm/2;
cfg.verticalMin   = -cfg.FOV_mm/2; cfg.verticalMax   = cfg.FOV_mm/2;
cfg.nPoints = 500;           % same grid size as the thesis (dx = 36 um)

%% Estimators (thesis settings)
cfg.kernel_size = 11;
cfg.applyTemporalWindow = true;
cfg.R2Threshold = 0.70;
cfg.MuTV = 250; cfg.TVIterations = 500; cfg.TVTau = 0.25; cfg.TVSigma = 0.25; cfg.TVTheta = 1.0;
cfg.MinK = 500; cfg.MaxK = 8000;
cfg.clim = [4 16];
cfg.MarginPx = 10;           % background disk margin (thesis)
% Near-field blanking, as in the experimental MMB / SBM processing
% (Cfg.NeedleRadius_mm): samples closer than R to the needle tip are set to
% NaN BEFORE the estimators, the same disk is removed from every map, and
% it is drawn gray in S-T maps / videos and darkest-turbo in elastograms.
% R = max(1.2 mm, 2*sqrt(2)*sigma_s) = source support (4.25 mm for the
% 3.54 mm needle, 1.2 mm = experiments' old-needle value for 0.30 mm).
cfg.Source_mm          = [0 0];                     % needle tip (scan centre)
cfg.NeedleRadius_mm    = AcqSim_NeedleRadius(cfg.SourceFWHM_mm);
cfg.NeedleExclusion_mm = cfg.NeedleRadius_mm;       % annulus ROI = blanked region (experimental metric ROI)
% Far-field ROI = R + half the 11 px PG kernel (~1 mm at 18 mm / 100 px):
% the first ~1 mm outside the blanked disk is biased by the interpolation
% bridge and by kernels straddling the disk edge (rim check, and the MB
% session's experimental-chain check: 17.2 -> 12.6 -> 10.85 m/s in 0.5 mm rings).
cfg.EdgeMargin_mm      = 1.0;

%% Excitation periods of the experiments [us]
cfg.Periods_us = [285 500 220 240 260 280 300 320 340];

%% ---------------- Multimeridian MB (QH, 3 cycles) ----------------------
cfg.MB.N_pulses      = 3;
cfg.MB.N_excitations = 2;
cfg.MB.M_repetitions = 450;  % A-lines per position (log)
cfg.MB.nSyncSamples  = 50;   % sync samples (log)
cfg.MB.nAscans       = 100;  % positions per meridian (log)
cfg.MB.ExcPulseWidth = 10;
cfg.MB.AcqPulseWidth = 1;
cfg.MB.startDelay_samples = 0;
% [period_us, meridians] pairs from the 06-10-2026 log
cfg.MB.Runs = [ ...
    500 2; 220 2; 285 2; ...
    285 4; 285 8; 285 16; 285 32; ...
    500 4; 500 8; 500 16; 500 32; ...
    220 8; 240 8; 260 8; 280 8; 300 8; 320 8; 340 8];

%% ---------------- SBM (harmonic, causal) --------------------------------
cfg.SBM.nAscans  = 37;
cfg.SBM.N_P      = 4;
cfg.SBM.MinTurn  = 20;
cfg.SBM.ExcPulseWidth = 10;
cfg.SBM.AcqPulseWidth = 10;
cfg.SBM.N_excitations = 2;   % sim duration only (harmonic source is continuous)
cfg.SBM.M_repetitions = 200; % sim-duration bookkeeping (thesis defaults)
cfg.SBM.nSyncSamples  = 50;
cfg.SBM.SteadyState_s = 1.5e-3; % clip start: after the wave filled the 18 mm domain
% [period_us, B-scans] pairs from the 19-10-2026 SBM log
cfg.SBM.Runs = [ ...
    285 21; ...
    285 11; 285 41; 285 81; ...
    500 21; 220 21; 240 21; 260 21; 280 21; 300 21; 320 21; 340 21];

end
