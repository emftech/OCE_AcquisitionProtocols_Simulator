%% Compare_NeedleSets_SimVsExp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 25-11-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Needle-size comparison and simulation-vs-experiment check for the
% thin-needle batch. Run after Run_AcqSim_Batch and Analyze_AcqSim_Batch.
%
% Source sets (every one is read from its own AcqSim_Config_used.mat):
%   - 12 x 12 mm, FWHM 0.65 mm (thin needle, R = 1.5 mm)   <- main
%   - 12 x 12 mm, FWHM 0.30 mm (point source,  R = 1.5 mm)
%   - 18 x 18 mm, FWHM 3.54 mm (legacy thick needle, R = 4.25 mm)
%   - 18 x 18 mm, FWHM 0.30 mm (legacy point source, R = 1.2 mm; same
%     geometry as the 14-11 experiments)
%
% 1) Speed bias vs frequency (sweep pattern: MB 8 meridians, SBM 21
%    B-scans), every estimator, c_ref = 10.54 m/s, on three radial ROIs
%    computed here on the stored (already blanked) maps, inside the thesis
%    background disk: r >= 1.5 mm, r >= 2.5 mm and the set's own far field
%    (r >= R + 1 mm). For the 3.54 mm set r < 4.25 mm is blanked, so its
%    1.5 / 2.5 mm ROIs equal r >= 4.25 mm.
% 2) Simulation vs experiment (14-11-2026, 0.65 mm needle, 18 x 18 mm,
%    2 meridians, 1 Vpp r1 + r2): 1D PG c(r) along the meridians (R2 >= 0.7,
%    0.5 mm bins, median, as Analyze_Tests_14_11_2026) and its mean over
%    the common radius 1.5-6 mm. Simulated runs: 2 meridians where the run
%    list has them (220 / 285 us at 12 mm; all three at 18 mm), else the
%    lowest meridian count (500 us at 12 mm: 4 meridians).
%
% Paths (environment variables, optional):
%   ACQSIM_LEGACY_RESULTS : folder holding FWHM_3p54_mm / FWHM_0p30_mm of
%                           the 18 x 18 mm batch (default <repo>/Results/Experiment_Matched_18x18mm)
%   ACQSIM_EXP_RESULTS    : semicolon-separated experimental Results folders
%                           with Summary/Summary_Tests_*.mat (c(r) table)
% Output: <results_path>/Comparison/ (.xlsx, .png, .fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all force;
batch_dir = string(fileparts(mfilename("fullpath")));
addpath(batch_dir);
cfg0 = AcqSim_Config(); AcqSim_SetupPaths(cfg0);
set(groot, 'DefaultFigureVisible', 'off');
c_ref = cfg0.c_ref;
out = cfg0.results_path + "/Comparison"; if ~isfolder(out), mkdir(out); end

legacy = string(getenv('ACQSIM_LEGACY_RESULTS'));
if legacy == "", legacy = cfg0.root_path + "/Results/Experiment_Matched_18x18mm"; end
expRoots = string(getenv('ACQSIM_EXP_RESULTS'));
if expRoots == "", expRoots = strings(0); else, expRoots = split(expRoots, ';').'; end

Sets = struct('Tag', {"0.65 mm, 12 mm", "0.30 mm, 12 mm", "3.54 mm, 18 mm", "0.30 mm, 18 mm"}, ...
    'Path', {cfg0.results_path + "/FWHM_0p65_mm", cfg0.results_path + "/FWHM_0p30_mm", ...
             legacy + "/FWHM_3p54_mm", legacy + "/FWHM_0p30_mm"});
Sets = Sets(arrayfun(@(s) isfile(s.Path + "/AcqSim_Config_used.mat"), Sets));
setCol = [0.80 0.15 0.15; 0.20 0.45 0.80; 0.45 0.45 0.45; 0.55 0.75 0.95];
rEdges = 0:0.5:10; rMid = rEdges(1:end-1) + 0.25;

%% 1) Bias vs frequency, all estimators, three ROIs
rows = {}; Prof = struct('Set', {}, 'Proto', {}, 'Period_us', {}, 'Pattern', {}, 'Estimator', {}, 'c_r', {});
for s = 1:numel(Sets)
    C = load(Sets(s).Path + "/AcqSim_Config_used.mat", "cfg"); cs = C.cfg;
    if ~isfield(cs, 'FOV_mm'), cs.FOV_mm = cs.horizontalMax - cs.horizontalMin; end   % configs older than the FOV field
    if ~isfield(cs, 'EdgeMargin_mm'), cs.EdgeMargin_mm = 1.0; end
    for proto = ["MB" "SBM"]
        if proto == "MB"
            files = dir(fullfile(Sets(s).Path, "Multimeridian_MB", "**", "MB_SimResults_*.mat"));
            est = ["MB_1DPG" "MB_1DPG_TV" "MB_2DPG" "MB_2DPG_TV"];
        else
            files = dir(fullfile(Sets(s).Path, "SBM", "**", "SBM_SimResults_*.mat"));
            est = ["SBM_2DPG" "SBM_2DPG_TV"];
        end
        for i = 1:numel(files)
            S = load(fullfile(files(i).folder, files(i).name), 'Results'); R = S.Results;
            Tus = round(R.Parameters.T_excitation*1e6);
            if proto == "MB", pat = R.Metadata.Number_Meridians; else, pat = R.Metadata.nBscans; end
            for e = est
                if ~isfield(R, e), continue; end
                M = R.(e); E = M.Elastogram_Masked; x = M.Xaxis_mm(:).'; y = M.Yaxis_mm(:);
                if isequal(size(E), [numel(y) numel(x)]), [X, Y] = meshgrid(x, y); else, [X, Y] = ndgrid(x, y); end
                r = hypot(X - cs.Source_mm(1), Y - cs.Source_mm(2));
                inB = r <= cs.FOV_mm/2 - cs.MarginPx*cs.FOV_mm/cs.nPoints;      % thesis background disk
                ok = isfinite(E) & inB;
                v15 = E(ok & r >= 1.5); v25 = E(ok & r >= 2.5); vF = E(ok & r >= cs.NeedleRadius_mm + cs.EdgeMargin_mm);
                v156 = E(ok & r >= 1.5 & r <= 6);
                rows(end+1, :) = {Sets(s).Tag, proto, Tus, 1e6/Tus, pat, erase(e, ["MB_" "SBM_"]), cs.NeedleRadius_mm, ...
                    mean(v15), std(v15), mean(v25), std(v25), mean(vF), std(vF), median(vF), mean(v156), ...
                    100*(mean(v15)/c_ref - 1), 100*(mean(v25)/c_ref - 1), 100*(mean(vF)/c_ref - 1)}; %#ok<AGROW>
                % c(r): 1D PG on the meridian samples (as the experiments), else radial median of the map
                if e == "MB_1DPG" && isfield(M, 'KVectorField') && isfield(M.KVectorField, 'Kmag_Meridian')
                    K = M.KVectorField; cm = 2*pi*(1e6/Tus)./K.Kmag_Meridian;
                    if median(cm(isfinite(cm)), 'omitnan') < 0.1, cm = cm*1e3; end   % rad/mm -> m/s guard
                    rr = hypot(K.X_mer_mm - cs.Source_mm(1), K.Y_mer_mm - cs.Source_mm(2));
                    okm = isfinite(cm) & K.R2_Meridian >= cs.R2Threshold & cm > 0.5 & cm < 40 & rr >= cs.NeedleRadius_mm;
                    cr = binMedian(rr(okm), cm(okm), rEdges);
                else
                    cr = binMedian(r(ok), E(ok), rEdges);
                end
                Prof(end+1) = struct('Set', Sets(s).Tag, 'Proto', proto, 'Period_us', Tus, 'Pattern', pat, ...
                    'Estimator', erase(e, ["MB_" "SBM_"]), 'c_r', cr); %#ok<AGROW>
            end
        end
    end
end
B = cell2table(rows, 'VariableNames', {'Set','Protocol','Period_us','f_Hz','Pattern','Estimator','R_mm', ...
    'Mean_r1p5','Std_r1p5','Mean_r2p5','Std_r2p5','Mean_Far','Std_Far','Median_Far','Mean_r1p5to6', ...
    'Bias_r1p5_pct','Bias_r2p5_pct','Bias_Far_pct'});
B.Set = string(B.Set); B.Protocol = string(B.Protocol); B.Estimator = string(B.Estimator);
B = sortrows(B, {'Protocol','Set','Estimator','Pattern','Period_us'});
writetable(B, out + "/Bias_AllSets_AllRuns.xlsx");
save(out + "/Comparison_Data.mat", "B", "Prof", "rEdges", "rMid", "Sets");

% Bias vs frequency figures (sweep pattern)
for proto = ["MB" "SBM"]
    if proto == "MB", sp = 8; ests = ["1DPG" "1DPG_TV" "2DPG" "2DPG_TV"]; else, sp = 21; ests = ["2DPG" "2DPG_TV"]; end
    rois = ["Bias_r1p5_pct" "Bias_r2p5_pct" "Bias_Far_pct"]; roiL = ["r \geq 1.5 mm" "r \geq 2.5 mm" "far field (r \geq R + 1 mm)"];
    f = figure('Color', 'w', 'Position', [50 50 420*numel(ests) 330*3]);
    tl = tiledlayout(f, 3, numel(ests), 'TileSpacing', 'compact', 'Padding', 'compact');
    for q = 1:3
        for e = 1:numel(ests)
            ax = nexttile(tl); hold(ax, 'on');
            for s = 1:numel(Sets)
                d = sortrows(B(B.Protocol == proto & B.Set == Sets(s).Tag & B.Pattern == sp & B.Estimator == ests(e), :), 'f_Hz');
                if isempty(d), continue; end
                plot(ax, d.f_Hz/1e3, d.(rois(q)), '-o', 'Color', setCol(s,:), 'MarkerFaceColor', setCol(s,:), ...
                    'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', Sets(s).Tag);
            end
            yline(ax, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            set(ax, 'FontWeight', 'bold', 'FontSize', 10, 'LineWidth', 1.2); grid(ax, 'on'); box(ax, 'on');
            if q == 1, title(ax, strrep(ests(e), '_', '-')); end
            if e == 1, ylabel(ax, {char(roiL(q)), 'bias vs c_{ref} [%]'}); end
            if q == 3, xlabel(ax, 'Excitation frequency [kHz]'); end
            if q == 1 && e == numel(ests), legend(ax, 'Location', 'best', 'FontSize', 8); end
        end
    end
    if proto == "MB", pn = "meridians"; else, pn = "B-scans"; end
    title(tl, sprintf('%s simulation - speed bias vs frequency (%d %s), c_{ref} = %.2f m/s', proto, sp, pn, c_ref), 'FontWeight', 'bold');
    savefig(f, out + "/Bias_vs_Frequency_" + proto + ".fig"); exportgraphics(f, out + "/Bias_vs_Frequency_" + proto + ".png", 'Resolution', 150); close(f);
end

% Compact bias table: sweep pattern, mean over the 9 frequencies and at 220 / 285 / 500 us
trow = {};
for proto = ["MB" "SBM"]
    if proto == "MB", sp = 8; else, sp = 21; end
    for s = 1:numel(Sets)
        for e = unique(B.Estimator(B.Protocol == proto), 'stable').'
            d = B(B.Protocol == proto & B.Set == Sets(s).Tag & B.Pattern == sp & B.Estimator == e, :);
            if isempty(d), continue; end
            at = @(p, fld) mean(d.(fld)(d.Period_us == p));
            trow(end+1, :) = {proto, Sets(s).Tag, e, height(d), mean(d.Bias_r1p5_pct), mean(d.Bias_r2p5_pct), mean(d.Bias_Far_pct), ...
                min(d.Bias_Far_pct), max(d.Bias_Far_pct), at(220, 'Bias_Far_pct'), at(285, 'Bias_Far_pct'), at(500, 'Bias_Far_pct')}; %#ok<AGROW>
        end
    end
end
BT = cell2table(trow, 'VariableNames', {'Protocol','Set','Estimator','nFreq','MeanBias_r1p5_pct','MeanBias_r2p5_pct', ...
    'MeanBias_Far_pct','MinBias_Far_pct','MaxBias_Far_pct','Bias_Far_220us_pct','Bias_Far_285us_pct','Bias_Far_500us_pct'});
writetable(BT, out + "/Bias_Table_SweepPattern.xlsx");
disp(BT);

%% 2) Simulation vs experiment: 1D PG c(r)
E = loadExperiments(expRoots);
per = [220 285 500];
f = figure('Color', 'w', 'Position', [50 50 1500 470]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
srow = {};
for p = 1:3
    ax = nexttile(tl); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    for s = 1:numel(Sets)
        P = Prof(string({Prof.Proto}) == "MB" & [Prof.Period_us] == per(p) & string({Prof.Estimator}) == "1DPG" & string({Prof.Set}) == Sets(s).Tag);
        if isempty(P), continue; end
        [~, k] = min([P.Pattern]); P = P(k);
        plot(ax, rMid, P.c_r, '-', 'Color', setCol(s,:), 'LineWidth', 2.2, ...
            'DisplayName', sprintf('sim %s (%d mer.)', Sets(s).Tag, P.Pattern));
        m = rMid >= 1.5 & rMid <= 6;
        srow(end+1, :) = {"Simulation", Sets(s).Tag, per(p), P.Pattern, mean(P.c_r(m), 'omitnan'), nnz(isfinite(P.c_r(m)))}; %#ok<AGROW>
    end
    for j = 1:numel(E)
        k = find(E(j).Period_us == per(p) & E(j).Vpp == 1);
        for i = k(:).'
            plot(ax, E(j).rMid, E(j).cr(i,:), '--o', 'Color', [0 0 0] + 0.15*(E(j).Rep(i) - 1), 'MarkerFaceColor', 'k', ...
                'MarkerSize', 4, 'LineWidth', 1.4, 'DisplayName', sprintf('exp %s 1 Vpp r%d', E(j).Tag, E(j).Rep(i)));
            m = E(j).rMid >= 1.5 & E(j).rMid <= 6;
            srow(end+1, :) = {"Experiment " + E(j).Tag, sprintf('r%d', E(j).Rep(i)), per(p), 2, mean(E(j).cr(i, m), 'omitnan'), nnz(isfinite(E(j).cr(i, m)))}; %#ok<AGROW>
        end
    end
    yline(ax, c_ref, 'k:', 'c_R = 10.54', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xline(ax, 1.5, ':', 'R_{nf} = 1.5 mm', 'HandleVisibility', 'off'); xline(ax, 6, ':', 'r = 6 mm (12 mm FOV edge)', 'HandleVisibility', 'off');
    xlim(ax, [0 10]); ylim(ax, [4 16]);
    set(ax, 'FontWeight', 'bold', 'FontSize', 11, 'LineWidth', 1.2);
    title(ax, sprintf('%d \\mus (%.2f kHz)', per(p), 1e3/per(p)));
    xlabel(ax, 'Distance from needle tip (mm)'); ylabel(ax, '1D PG speed (m/s)');
    if p == 3, legend(ax, 'Location', 'southwest', 'FontSize', 8); end
end
title(tl, 'Simulation vs experiment - 1D PG speed vs distance from the needle (R^2 \geq 0.7, 0.5 mm bins, median)', 'FontWeight', 'bold');
savefig(f, out + "/SimVsExp_c_of_r.fig"); exportgraphics(f, out + "/SimVsExp_c_of_r.png", 'Resolution', 150); close(f);
SvE = cell2table(srow, 'VariableNames', {'Source','Set','Period_us','Meridians','Mean_c_1p5to6mm','nBins'});
writetable(SvE, out + "/SimVsExp_Mean_1p5to6mm.xlsx");
disp(SvE);
fprintf('Comparison written to %s\n', out);

%% ------------------------------------------------------------------------
function cr = binMedian(r, c, edges)
cr = nan(1, numel(edges) - 1);
for b = 1:numel(edges) - 1
    m = r >= edges(b) & r < edges(b+1);
    if nnz(m) >= 3, cr(b) = median(c(m)); end
end
end

function E = loadExperiments(roots)
E = struct('Tag', {}, 'Period_us', {}, 'Vpp', {}, 'Rep', {}, 'cr', {}, 'rMid', {});
for k = 1:numel(roots)
    f = dir(fullfile(roots(k), "Summary", "Summary_Tests_*.mat"));
    if isempty(f), continue; end
    L = load(fullfile(f(1).folder, f(1).name));
    tag = string(regexp(f(1).name, '\d+_\d+', 'match', 'once'));
    E(end+1) = struct('Tag', strrep(tag, '_', '-'), 'Period_us', L.Period_us(:), 'Vpp', L.Vpp(:), ...
        'Rep', L.Rep(:), 'cr', L.cr, 'rMid', L.rMid); %#ok<AGROW>
end
end
