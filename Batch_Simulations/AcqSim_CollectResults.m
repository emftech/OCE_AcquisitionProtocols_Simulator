function [T, Maps] = AcqSim_CollectResults(set_path, protocol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collects every *_SimResults_*.mat of one source set and protocol into a
% long table (one row per run x estimator) plus the masked elastograms.
%
%   T columns: Period_us, f_Hz, Pattern (meridians or B-scans), nTurn,
%   Estimator, <ROI>_mean/_std/_Ea/_Ep/_Cov for ROI = Bkg, Ann, Far,
%   Ea_thesis, Ep_thesis (Compute_* outputs on the thesis disk),
%   OG_Far (2D PG on the full wave field), f_bin (FFT bin used), File.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
protocol = upper(string(protocol));
switch protocol
    case "MB",  root = fullfile(set_path, "Multimeridian_MB"); pat = "MB_SimResults_*.mat";
                est = ["MB_1DPG","MB_1DPG_TV","MB_2DPG","MB_2DPG_TV"];
    case "SBM", root = fullfile(set_path, "SBM"); pat = "SBM_SimResults_*.mat";
                est = ["SBM_2DPG","SBM_2DPG_TV"];
end
files = dir(fullfile(root, "**", pat));
rows = {}; Maps = struct('key', {}, 'Estimator', {}, 'E', {}, 'x', {}, 'y', {});
for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name), 'Results');
    R = S.Results;
    Tus = round(R.Parameters.T_excitation*1e6);
    if protocol == "MB"
        pattern = R.Metadata.Number_Meridians; nTurn = NaN;
        fs = 1/R.Parameters.Ts_OCT;
    else
        pattern = R.Metadata.nBscans; nTurn = R.Metadata.nTurnSamples;
        fs = 1/R.Parameters.Ts_SBM;
    end
    OG = R.OG_WaveProp;
    ogFar = NaN; if isfield(OG, 'FarField_median'), ogFar = OG.FarField_median; end
    for e = est
        if ~isfield(R, e), continue; end
        M = R.(e);
        bk = local_roi(M.Elastogram_Masked, R.GroundTruth.Bkg_Mask, R.Parameters.c_medium);
        if contains(e, "1DPG"), Nfft = 2^nextpow2(1024); else, Nfft = 2^10 + 1; end
        if contains(e, "1DPG"), fsb = 1/R.Parameters.Ts_OCT; else, fsb = fs; end
        fb = (0:Nfft-1)*fsb/Nfft; [~, kb] = min(abs(fb - 1/R.Parameters.T_excitation)); fbin = fb(kb);
        rows(end+1, :) = {Tus, 1e6/Tus, pattern, nTurn, erase(e, ["MB_","SBM_"]), ...
            bk.mean, bk.std, bk.Ea, bk.Ep, bk.Coverage, ...
            M.Annulus.mean, M.Annulus.std, M.Annulus.Ea, M.Annulus.Ep, M.Annulus.Coverage, ...
            M.FarField.mean, M.FarField.std, M.FarField.Ea, M.FarField.Ep, M.FarField.Coverage, ...
            M.AccuracyError.Ea_Bkg, M.PrecisionError.Ep_Bkg, ogFar, fbin, string(files(i).name)}; %#ok<AGROW>
        Maps(end+1) = struct('key', sprintf('%d_%d', Tus, pattern), 'Estimator', erase(e, ["MB_","SBM_"]), ...
            'E', M.Elastogram_Masked, 'x', M.Xaxis_mm, 'y', M.Yaxis_mm); %#ok<AGROW>
    end
end
vn = {'Period_us','f_Hz','Pattern','nTurn','Estimator', ...
      'Bkg_mean','Bkg_std','Bkg_Ea','Bkg_Ep','Bkg_Cov', ...
      'Ann_mean','Ann_std','Ann_Ea','Ann_Ep','Ann_Cov', ...
      'Far_mean','Far_std','Far_Ea','Far_Ep','Far_Cov', ...
      'Ea_thesis','Ep_thesis','OG_Far','f_bin','File'};
if isempty(rows), T = table(); return; end
T = cell2table(rows, 'VariableNames', vn);
T = sortrows(T, {'Period_us','Pattern','Estimator'});
end

function S = local_roi(E, Mask, c_ref)
Mask = imresize(logical(Mask), size(E), 'nearest');
v = E(Mask & isfinite(E));
S = struct('mean', mean(v), 'std', std(v), 'Ea', abs(mean(v)-c_ref)/c_ref*100, ...
    'Ep', std(v)/mean(v)*100, 'Coverage', numel(v)/max(1,nnz(Mask))*100);
end
