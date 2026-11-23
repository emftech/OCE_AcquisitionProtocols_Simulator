%% Run_AcqSim_Batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk, experiment-matched simulations of the Multimeridian MB mode (QH,
% 3 cycles) and SBM (causal harmonic) protocols on an 18 x 18 mm
% homogeneous medium (c = 10.54 m/s). One k-Wave video per period and
% protocol, then every acquisition pattern of the logs for that period.
%
% Source sets:
%   - FWHM 3.54 mm : new needle (main results)
%   - FWHM 0.30 mm : thesis source (control: isolates the needle-size effect)
%
% After this script: Analyze_AcqSim_Batch.m, then Build_Deck_MB_Sim.m and
% Build_Deck_SBM_Sim.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all force;
batch_dir = string(fileparts(mfilename("fullpath")));
addpath(batch_dir);
cfg0 = AcqSim_Config();
AcqSim_SetupPaths(cfg0);
set(groot, 'DefaultFigureVisible', 'off');

% Optional split into parallel workers (environment variables):
%   ACQSIM_SETS="3.54" | "0.30" | "3.54,0.30"   ACQSIM_PROTOCOLS="MB" | "SBM" | "MB,SBM"
% Runs whose last output already exists are skipped (resume after a crash).
SourceSets_mm = [3.54 0.30];
RunProtocols  = ["MB" "SBM"];
if ~isempty(getenv('ACQSIM_SETS')),      SourceSets_mm = str2double(split(string(getenv('ACQSIM_SETS')), ',')).'; end
if ~isempty(getenv('ACQSIM_PROTOCOLS')), RunProtocols  = split(string(getenv('ACQSIM_PROTOCOLS')), ',').'; end
workerTag = strjoin(string(SourceSets_mm), '_') + "_" + strjoin(RunProtocols, '_');

t_all = tic;
Log = table();
for s = 1:numel(SourceSets_mm)
    cfg = cfg0;
    cfg.SourceFWHM_mm = SourceSets_mm(s);
    cfg.NeedleRadius_mm = AcqSim_NeedleRadius(cfg.SourceFWHM_mm);
    cfg.NeedleExclusion_mm = cfg.NeedleRadius_mm;
    cfg.results_path  = cfg0.results_path + "/" + sprintf("FWHM_%s_mm", strrep(sprintf('%.2f', SourceSets_mm(s)), '.', 'p'));
    if ~isfolder(cfg.results_path), mkdir(cfg.results_path); end
    save(cfg.results_path + "/AcqSim_Config_used.mat", "cfg");

    for proto = RunProtocols
        Runs = cfg.(proto).Runs;
        periods = unique(Runs(:,1), 'stable');
        for ip = 1:numel(periods)
            T_us = periods(ip);
            sel = Runs(Runs(:,1) == T_us, 2);
            done = arrayfun(@(m) runDone(cfg, proto, T_us, m), sel);
            if all(done)
                fprintf('[%s] %d us: all %d runs already done, skipped\n', proto, T_us, numel(sel));
                continue;
            end
            sel = sel(~done);
            t_w = tic;
            W = AcqSim_WaveProp(cfg, T_us*1e-6, proto);
            fprintf('[%s] wave video %d us (FWHM %.2f mm) in %.1f s, OG far field %.2f m/s\n', ...
                proto, T_us, cfg.SourceFWHM_mm, toc(t_w), W.OG.FarField_median);
            for r = 1:numel(sel)
                try
                    if proto == "MB"
                        R = AcqSim_MB_SingleRun(cfg, W, sel(r), "batch");
                    else
                        R = AcqSim_SBM_SingleRun(cfg, W, sel(r), "batch");
                    end
                    ok = true; msg = string(R.Metadata.results_fullpath);
                catch ME
                    ok = false; msg = string(ME.message);
                    fprintf(2, '[%s] %d us / %d FAILED: %s\n', proto, T_us, sel(r), ME.message);
                end
                Log = [Log; table(SourceSets_mm(s), proto, T_us, sel(r), ok, msg, ...
                    'VariableNames', {'FWHM_mm','Protocol','Period_us','Pattern','OK','File'})]; %#ok<AGROW>
                close all force;
            end
            clear W R;
        end
    end
end
writetable(Log, cfg0.results_path + "/Batch_RunLog_" + workerTag + ".csv");
fprintf('\nBatch finished in %.1f min (%d runs, %d failed)\n', toc(t_all)/60, height(Log), sum(~Log.OK));

%% ------------------------------------------------------------------------
function tf = runDone(cfg, proto, T_us, pattern)
% A run is complete when its LAST output exists (MB: traces file; SBM: the
% acquisition-phase / S-T figure), both written after the results .mat and videos.
if proto == "MB"
    f = fullfile(cfg.results_path, "Multimeridian_MB", "Cylindrical", "Quasi_Harmonic", ...
        sprintf("%d Meridians", pattern), sprintf("Quasi_Harmonic_Ideal_%dus", T_us), ...
        sprintf("MB%dMer_Traces_%dus.mat", pattern, T_us));
    tf = isfile(f);
else
    P = round(T_us*1e-6/cfg.Ts_OCT);
    nT = ceil((cfg.SBM.nAscans + cfg.SBM.MinTurn)/P)*P - cfg.SBM.nAscans;
    d = fullfile(cfg.results_path, "SBM", sprintf("%d_Ascans_%d_BScans_%d_TurnSamples", cfg.SBM.nAscans, pattern, nT), ...
        sprintf("Harmonic_Ideal_Cylindrical_%dus", T_us));
    tf = ~isempty(dir(fullfile(d, "AcqPhase_and_ST_*.png")));
end
end
