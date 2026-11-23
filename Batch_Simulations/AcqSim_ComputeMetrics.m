function Metrics = AcqSim_ComputeMetrics(Est, GT_SpeedMap, Inc_Mask, Bkg_Mask, Ann_Mask, Far_Mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcqSim_ComputeMetrics
%
% Thesis metrics for every estimator in Est (fields with .E masked
% elastogram): Compute_AccuracyError_Masks_v2 and
% Compute_PrecisionError_Masks on the thesis background disk (Bkg_Mask),
% plus the same Ea (Eq. 28, mean based) and Ep (Eq. 29) on the annulus
% Ann_Mask that excludes the needle near field (experiment-style
% blanking) and on Far_Mask (outside the k-Wave source support,
% r >= R + 1 mm: blank radius plus PG-kernel half-width, edge-free). No CNR (no inclusion).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_ref = mode(round(GT_SpeedMap(Bkg_Mask),6));
Metrics = struct();
fn = fieldnames(Est);
for k = 1:numel(fn)
    E = Est.(fn{k}).E;
    [Ea_Bkg, Ea_Inc, Ea_Global, Ea_Map, Ea_Stats] = evalc_quiet(@() ...
        Compute_AccuracyError_Masks_v2(GT_SpeedMap, Inc_Mask, Bkg_Mask, E), 5);
    [Ep_Bkg, Ep_Inc, Ep_Global, Ep_ROI, Ep_Stats] = evalc_quiet(@() ...
        Compute_PrecisionError_Masks(E, Inc_Mask, Bkg_Mask), 5);
    M = struct();
    M.AccuracyError  = struct('Ea_Bkg', Ea_Bkg, 'Ea_Inc', Ea_Inc, 'Ea_Global', Ea_Global, ...
        'Ea_Map_LocalPixelwise', Ea_Map, 'Ea_Stats', Ea_Stats);
    M.PrecisionError = struct('Ep_Bkg', Ep_Bkg, 'Ep_Inc', Ep_Inc, 'Ep_Global', Ep_Global, ...
        'Ep_ROI_Masks', Ep_ROI, 'Ep_Stats', Ep_Stats);
    M.Bkg     = roi_stats(E, Bkg_Mask, c_ref);
    M.Annulus = roi_stats(E, Ann_Mask, c_ref);
    M.FarField = roi_stats(E, Far_Mask, c_ref);
    Metrics.(fn{k}) = M;
end
end

function S = roi_stats(E, Mask, c_ref)
if ~isequal(size(E), size(Mask))
    Mask = imresize(logical(Mask), size(E), 'nearest');
end
v = E(Mask & isfinite(E));
S = struct();
S.N        = numel(v);
S.Coverage = numel(v) / max(1, nnz(Mask)) * 100;   % [%] of ROI with R2 >= thr
S.mean     = mean(v);
S.std      = std(v);
S.median   = median(v);
S.Ea       = abs(S.mean - c_ref) / c_ref * 100;
S.Ep       = S.std / S.mean * 100;
S.c_ref    = c_ref;
end

function varargout = evalc_quiet(fh, nout)
% Run a framework function silencing its command-window report.
varargout = cell(1, nout);
[~, varargout{:}] = evalc_wrap(fh, nout);
end

function [txt, varargout] = evalc_wrap(fh, nout)
out = cell(1, nout);
txt = evalc('[out{:}] = fh();');
varargout = out;
end
