function [SampleValues, SampleInfo] = Main_SamplingfromAcqvecs_LoopFunction( ...
    WaveProp_Vid, x_grid, y_grid, t_grid, ...
    GalvoX, GalvoY, Acq_Trigger, Exc_Trigger, ...
    Ts_OCT, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 18-05-2026
% Last modified: 21-05-2026
%
% Generalized low-level sampling loop.
%
% This function samples WaveProp_Vid(x,y,t) according to arbitrary
% acquisition vectors. It supports:
%
%   - arbitrary GalvoX/GalvoY trajectories,
%   - 3D interpolation in x, y and time,
%   - non-integer temporal delays,
%   - query-level saccadic offsets,
%   - temporal ensemble cycling.
%
% Saccades affect only the interpolation query:
%
%   x_query = x_nominal + delta_x
%   y_query = y_nominal + delta_y
%
% The nominal coordinates are stored separately and must be used for spatial
% reconstruction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options =============================

P = inputParser;

addParameter(P, 'SamplesPerSimEnsemble', []);
addParameter(P, 'InterpMethod', 'linear');
addParameter(P, 'ExtrapMethod', 'none');
addParameter(P, 'TemporalDelay_s', 0);
addParameter(P, 'QueryDeltaX', []);
addParameter(P, 'QueryDeltaY', []);

parse(P, varargin{:});
Opt = P.Results;

if isempty(Opt.SamplesPerSimEnsemble)
    error(['SamplesPerSimEnsemble must be provided. ', ...
           'For QH MB-like simulations, use M_repetitions + nSyncSamples.']);
end

%% ========================== Normalize inputs =============================

GalvoX      = double(GalvoX(:));
GalvoY      = double(GalvoY(:));
Acq_Trigger = double(Acq_Trigger(:));
Exc_Trigger = double(Exc_Trigger(:));

x_grid = double(x_grid(:));
y_grid = double(y_grid(:));
t_grid = double(t_grid(:));

N_vec = numel(GalvoX);

if numel(GalvoY) ~= N_vec || ...
   numel(Acq_Trigger) ~= N_vec || ...
   numel(Exc_Trigger) ~= N_vec
    error('GalvoX, GalvoY, Acq_Trigger and Exc_Trigger must have equal length.');
end

[Nx, Ny, Nt_sim] = size(WaveProp_Vid);

if numel(x_grid) ~= Nx
    error('numel(x_grid) must match size(WaveProp_Vid,1).');
end

if numel(y_grid) ~= Ny
    error('numel(y_grid) must match size(WaveProp_Vid,2).');
end

if numel(t_grid) ~= Nt_sim
    error('numel(t_grid) must match size(WaveProp_Vid,3).');
end

%% ========================== Unit harmonization ===========================
% Existing acquisition vectors are usually in mm, while k-Wave grids are
% usually in meters. If the simulation grid appears to be in meters, convert
% GalvoX/GalvoY from mm to m.

gridLikelyMeters = max(abs(x_grid)) < 1 && max(abs(y_grid)) < 1;

if gridLikelyMeters
    GalvoX_nominal = GalvoX * 1e-3;
    GalvoY_nominal = GalvoY * 1e-3;
    spatialScaleApplied = 1e-3;
else
    GalvoX_nominal = GalvoX;
    GalvoY_nominal = GalvoY;
    spatialScaleApplied = 1;
end

%% ========================== Query offsets ================================

if isempty(Opt.QueryDeltaX)
    QueryDeltaX = zeros(N_vec,1);
else
    QueryDeltaX = double(Opt.QueryDeltaX(:));
end

if isempty(Opt.QueryDeltaY)
    QueryDeltaY = zeros(N_vec,1);
else
    QueryDeltaY = double(Opt.QueryDeltaY(:));
end

if numel(QueryDeltaX) ~= N_vec || numel(QueryDeltaY) ~= N_vec
    error('QueryDeltaX and QueryDeltaY must be empty or match GalvoX length.');
end

%% ===================== Temporal ensemble cycling =========================
% Use:
%   - ensemble 1 if only one ensemble exists,
%   - ensemble 2 if two ensembles exist,
%   - ensembles 2:end cyclically if more than two ensembles exist.

samplesPerEnsemble_est = round(Opt.SamplesPerSimEnsemble);

if samplesPerEnsemble_est <= 0
    error('SamplesPerSimEnsemble must be positive.');
end

nSimEnsembles = floor(Nt_sim / samplesPerEnsemble_est);

if nSimEnsembles < 1
    error('The simulation movie does not contain a complete temporal ensemble.');
end

if nSimEnsembles == 1
    ensembleCycle = 1;
elseif nSimEnsembles == 2
    ensembleCycle = 2;
else
    ensembleCycle = 2:nSimEnsembles;
end

ensembleDuration_s = samplesPerEnsemble_est * Ts_OCT;

%% =========================== Build interpolant ===========================

F_wave = griddedInterpolant( ...
    {x_grid, y_grid, t_grid}, ...
    double(WaveProp_Vid), ...
    Opt.InterpMethod, ...
    Opt.ExtrapMethod);

%% ====================== Active acquisition samples =======================

activeIdx = find(Acq_Trigger > 0);
nSamples  = numel(activeIdx);

SampleValues = nan(nSamples,1);

%% ===================== Contiguous acquisition groups =====================

AcqGroupID_full = zeros(N_vec,1);

groupCounter = 0;
insideGroup  = false;

for k = 1:N_vec

    if Acq_Trigger(k) > 0 && ~insideGroup
        groupCounter = groupCounter + 1;
        insideGroup  = true;
    elseif Acq_Trigger(k) <= 0
        insideGroup = false;
    end

    if Acq_Trigger(k) > 0
        AcqGroupID_full(k) = groupCounter;
    end
end

%% =============================== Metadata ================================

SampleInfo = struct();

SampleInfo.ActiveIdx             = activeIdx;
SampleInfo.SampleSerial          = (1:nSamples).';

SampleInfo.Xq_nominal            = nan(nSamples,1);
SampleInfo.Yq_nominal            = nan(nSamples,1);
SampleInfo.Xq_moved              = nan(nSamples,1);
SampleInfo.Yq_moved              = nan(nSamples,1);
SampleInfo.DeltaX                = nan(nSamples,1);
SampleInfo.DeltaY                = nan(nSamples,1);

SampleInfo.Tq                    = nan(nSamples,1);
SampleInfo.AcqTime_s             = nan(nSamples,1);
SampleInfo.LocalEnsembleTime_s   = nan(nSamples,1);
SampleInfo.SimEnsembleID         = nan(nSamples,1);
SampleInfo.AcqGroupID            = nan(nSamples,1);

SampleInfo.SamplesPerSimEnsemble = samplesPerEnsemble_est;
SampleInfo.NumSimEnsembles       = nSimEnsembles;
SampleInfo.EnsembleCycle         = ensembleCycle;
SampleInfo.EnsembleDuration_s    = ensembleDuration_s;
SampleInfo.SpatialScaleApplied   = spatialScaleApplied;

%% ================================ Main loop ==============================

for s = 1:nSamples

    idx = activeIdx(s);

    acqTime_s = (idx - 1) * Ts_OCT;

    acqEnsembleNumber = floor(acqTime_s / ensembleDuration_s) + 1;
    localTime_s       = mod(acqTime_s, ensembleDuration_s);

    cycleIdx      = mod(acqEnsembleNumber - 1, numel(ensembleCycle)) + 1;
    simEnsembleID = ensembleCycle(cycleIdx);

    tq = (simEnsembleID - 1) * ensembleDuration_s + ...
          localTime_s + Opt.TemporalDelay_s;

    x_nom = GalvoX_nominal(idx);
    y_nom = GalvoY_nominal(idx);

    x_mov = x_nom + QueryDeltaX(idx);
    y_mov = y_nom + QueryDeltaY(idx);

    x_mov = min(max(x_mov, min(x_grid)), max(x_grid));
    y_mov = min(max(y_mov, min(y_grid)), max(y_grid));

    SampleInfo.Xq_nominal(s)          = x_nom;
    SampleInfo.Yq_nominal(s)          = y_nom;
    SampleInfo.Xq_moved(s)            = x_mov;
    SampleInfo.Yq_moved(s)            = y_mov;
    SampleInfo.DeltaX(s)              = QueryDeltaX(idx);
    SampleInfo.DeltaY(s)              = QueryDeltaY(idx);

    SampleInfo.Tq(s)                  = tq;
    SampleInfo.AcqTime_s(s)           = acqTime_s;
    SampleInfo.LocalEnsembleTime_s(s) = localTime_s;
    SampleInfo.SimEnsembleID(s)       = simEnsembleID;
    SampleInfo.AcqGroupID(s)          = AcqGroupID_full(idx);

    SampleValues(s) = F_wave(x_mov, y_mov, tq);
end

end