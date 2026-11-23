%% Analyze_AcqSim_Batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metrics (accuracy Ea, precision Ep; no CNR: homogeneous medium) and
% dispersion analysis of the experiment-matched simulation batch.
%
% For every source set (FWHM_3p54_mm = new needle, FWHM_0p30_mm = thesis
% source) and protocol it writes to <set>/Summary_<protocol>/:
%   - <protocol>_Summary_All.xlsx  (one row per run x estimator, 3 ROIs)
%   - Metrics_vs_Meridians / Metrics_vs_BScans  (.xlsx + thesis-style plots)
%   - Elastogram grids (rows = pattern, columns = estimator), shared turbo range
%   - Dispersion.xlsx + plot: c(f) at the 9 experiment frequencies with the
%     two-range linear fits used for the experiments (2.0-3.6 and 3.5-4.6 kHz)
%     and the FFT-bin prediction c_ref*f_exc/f_bin (estimator-only baseline)
%   - NearField_Radius.xlsx + plot: radius beyond which the full-field 2D PG
%     stays within 5 % of c_ref (per frequency), i.e. a blanking radius
% and, at the batch root, FWHM_Comparison plots (needle vs thesis source).
%
% ROIs: Bkg = thesis disk (R = N/2 - 10 px); Ann = Far = r >= NeedleRadius_mm
% = max(1.2, 2.83 sigma_s) (experimental near-field blanking radius = source
% support). Maps are already blanked inside that disk by the single runs;
% the near-field radius analysis uses the UNBLANKED Results.OG_WaveProp.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all force;
batch_dir = string(fileparts(mfilename("fullpath")));
addpath(batch_dir);
cfg = AcqSim_Config(); AcqSim_SetupPaths(cfg);
set(groot, 'DefaultFigureVisible', 'off');

sets = dir(fullfile(cfg.results_path, "FWHM_*_mm"));
AllT = struct();
for s = 1:numel(sets)
    set_path = fullfile(sets(s).folder, sets(s).name);
    Sc = load(fullfile(set_path, "AcqSim_Config_used.mat"), "cfg"); cfg = Sc.cfg;   % set-specific NeedleRadius_mm / FWHM
    setTag = string(sets(s).name);
    setLbl = strrep(strrep(setTag, '_', ' '), 'p', '.');      % e.g. "FWHM 3.54 mm"
    for proto = ["MB" "SBM"]
        [T, Maps] = AcqSim_CollectResults(set_path, proto);
        if isempty(T), continue; end
        out = fullfile(set_path, "Summary_" + proto);
        if ~isfolder(out), mkdir(out); end
        writetable(T, fullfile(out, proto + "_Summary_All.xlsx"));
        AllT.(matlab.lang.makeValidName(setTag + "_" + proto)) = T;

        if proto == "MB"
            patList = [4 8 16 32]; patName = "Meridians"; sweepPattern = 8;
            gridPeriods = [285 500];
        else
            patList = [11 21 41 81]; patName = "BScans"; sweepPattern = 21;
            gridPeriods = 285;
        end
        ests = unique(T.Estimator, 'stable');

        %% Metrics vs pattern (thesis Fig. 42 style)
        for Tus = gridPeriods
            sub = T(T.Period_us == Tus & ismember(T.Pattern, patList), :);
            if isempty(sub), continue; end
            writetable(sub, fullfile(out, sprintf("Metrics_vs_%s_%dus.xlsx", patName, Tus)));
            for roi = ["Far" "Ann" "Bkg"]
                f = figure('Color','w','Position',[100 100 1300 420]);
                tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
                qty = [roi+"_mean", roi+"_Ea", roi+"_Ep"];
                ylab = ["Mean phase speed [m/s]", "E_a [%]", "E_p [%]"];
                C = lines(numel(ests));
                for q = 1:3
                    ax = nexttile(tl); hold(ax, 'on');
                    for e = 1:numel(ests)
                        se = sub(sub.Estimator == ests(e), :); se = sortrows(se, 'Pattern');
                        plot(ax, se.Pattern, se.(qty(q)), '-o', 'Color', C(e,:), 'LineWidth', 2.5, ...
                            'MarkerSize', 8, 'DisplayName', strrep(ests(e), '_', '-'));
                    end
                    if q == 1, yline(ax, cfg.c_ref, 'k--', 'LineWidth', 1.5, 'DisplayName', 'c_{ref}'); end
                    set(ax, 'XScale', 'log', 'XTick', patList, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
                    xlabel(ax, "Number of " + patName, 'FontWeight', 'bold'); ylabel(ax, ylab(q), 'FontWeight', 'bold');
                    grid(ax, 'on'); box(ax, 'on');
                    if q == 1, legend(ax, 'Location', 'best'); end
                end
                title(tl, sprintf('%s simulation - %d \\mus - ROI: %s - %s', proto, Tus, roiName(roi, cfg), setLbl), 'FontWeight', 'bold');
                fb = fullfile(out, sprintf("Metrics_vs_%s_%dus_%s", patName, Tus, roi));
                savefig(f, fb + ".fig"); exportgraphics(f, fb + ".png", 'Resolution', 150); close(f);
            end

            %% Elastogram grid (experiment helper: shared range, darkest turbo = NaN, white '+' at the tip)
            Mk = Maps(startsWith(string({Maps.key}), sprintf('%d_', Tus)));
            xc = linspace(cfg.horizontalMin, cfg.horizontalMax, 201); yc = linspace(cfg.verticalMin, cfg.verticalMax, 201);
            MapCell = cell(numel(patList), numel(ests));
            for ip = 1:numel(patList)
                for e = 1:numel(ests)
                    k = find(string({Mk.key}) == sprintf('%d_%d', Tus, patList(ip)) & string({Mk.Estimator}) == string(ests(e)), 1);
                    if isempty(k), continue; end
                    MapCell{ip, e} = toCommonGrid(Mk(k).E, Mk(k).x, Mk(k).y, xc, yc);
                end
            end
            f = MMB_PlotElastogramGrid(MapCell, xc, yc, cellstr(string(patList(:))), cellstr(strrep(ests(:).', '_', '-')), ...
                MMB_SpeedCLim(MapCell), sprintf('%s simulation elastograms (R^2 \\geq 0.7, r < %.2f mm blanked) - %d \\mus - %s', ...
                proto, cfg.NeedleRadius_mm, Tus, setLbl), 'Source_mm', cfg.Source_mm, 'RowSuffix', char(patName));
            set(f, 'Visible', 'off');
            fb = fullfile(out, sprintf("ElastogramGrid_%dus", Tus));
            savefig(f, fb + ".fig"); exportgraphics(f, fb + ".png", 'Resolution', 150); close(f);
        end

        %% Dispersion (pattern of the experiment frequency sweep)
        D = T(T.Pattern == sweepPattern, :);
        D = sortrows(D, {'Estimator','f_Hz'});
        D.Bin_prediction = cfg.c_ref * D.f_Hz ./ D.f_bin;
        fitRows = {};
        ranges = [2000 3600; 3500 4600; 0 1e5];
        rangeNames = ["2.0-3.6 kHz", "3.5-4.6 kHz", "all"];
        for e = 1:numel(ests)
            for roi = ["Far" "Ann"]
                se = D(D.Estimator == ests(e), :);
                for r = 1:3
                    m = se.f_Hz >= ranges(r,1) & se.f_Hz <= ranges(r,2) & isfinite(se.(roi+"_mean"));
                    if nnz(m) < 3, continue; end
                    [slope, R2, p] = linfit(se.f_Hz(m)/1e3, se.(roi+"_mean")(m));
                    fitRows(end+1,:) = {ests(e), roi, rangeNames(r), nnz(m), slope, R2, p}; %#ok<AGROW>
                end
            end
        end
        Fits = cell2table(fitRows, 'VariableNames', {'Estimator','ROI','Range','N','Slope_ms_per_kHz','R2','p'});
        writetable(D, fullfile(out, "Dispersion.xlsx"), 'Sheet', 'Data');
        writetable(Fits, fullfile(out, "Dispersion.xlsx"), 'Sheet', 'Fits');
        for roi = ["Far" "Ann"]
            f = figure('Color','w','Position',[100 100 900 560]); ax = axes(f); hold(ax, 'on');
            C = lines(numel(ests));
            for e = 1:numel(ests)
                se = D(D.Estimator == ests(e), :);
                errorbar(ax, se.f_Hz/1e3, se.(roi+"_mean"), se.(roi+"_std"), '-o', 'Color', C(e,:), ...
                    'LineWidth', 2.5, 'MarkerSize', 8, 'CapSize', 6, 'DisplayName', strrep(ests(e), '_', '-'));
            end
            og = unique(D(:, {'f_Hz','OG_Far'}), 'rows');
            plot(ax, og.f_Hz/1e3, og.OG_Far, 's--', 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'MarkerSize', 9, ...
                'DisplayName', '2D PG on full wave field (far field)');
            if proto == "MB"
                bp = D(D.Estimator == "2DPG", :);
                plot(ax, bp.f_Hz/1e3, bp.Bin_prediction, ':', 'Color', [0.6 0 0.6], 'LineWidth', 2, ...
                    'DisplayName', 'c_{ref} f_{exc}/f_{bin} (2D FFT-bin baseline)');
            end
            yline(ax, cfg.c_ref, 'k--', 'LineWidth', 1.5, 'DisplayName', 'c_{ref} = 10.54 m/s');
            xlabel(ax, 'Excitation frequency [kHz]', 'FontWeight', 'bold'); ylabel(ax, 'Phase speed [m/s]', 'FontWeight', 'bold');
            title(ax, {sprintf('%s simulation dispersion (%d %s) - ROI: %s', proto, sweepPattern, patName, roiName(roi, cfg)), ...
                sprintf('Homogeneous non-dispersive medium: any slope is estimator / sampling / source bias - %s', setLbl)}, ...
                'FontWeight', 'bold');
            set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5); grid(ax, 'on'); box(ax, 'on');
            legend(ax, 'Location', 'best');
            fb = fullfile(out, "Dispersion_" + roi);
            savefig(f, fb + ".fig"); exportgraphics(f, fb + ".png", 'Resolution', 150); close(f);
        end

        %% Near-field radius from the full-field 2D PG (one per frequency)
        NF = nearFieldRadius(T, cfg);
        writetable(NF, fullfile(out, "NearField_Radius.xlsx"));
        f = figure('Color','w','Position',[100 100 760 480]); ax = axes(f); hold(ax,'on');
        plot(ax, NF.f_Hz/1e3, NF.R5_mm, 'o', 'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', 'r where |c-c_{ref}| < 5 % (full-field 2D PG)');
        plot(ax, NF.f_Hz/1e3, NF.Lambda_mm, '--', 'LineWidth', 2, 'DisplayName', '\lambda = c_{ref}/f');
        yline(ax, 2*sqrt(2)*cfg.SourceFWHM_mm/(2*sqrt(2*log(2))), ':', 'LineWidth', 1.5, ...
            'DisplayName', 'source support radius 2.83\sigma_s = blanking radius');
        xlabel(ax, 'Excitation frequency [kHz]', 'FontWeight', 'bold'); ylabel(ax, 'Radius [mm]', 'FontWeight', 'bold');
        title(ax, {sprintf('%s - near-field (blanking) radius - %s', proto, setLbl), ...
            'missing points: never within 5 % up to r = 7.6 mm (global bias, see dispersion)'}, 'FontWeight', 'bold');
        legend(ax, 'Location', 'best'); grid(ax, 'on'); set(ax, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 1.5);
        fb = fullfile(out, "NearField_Radius");
        savefig(f, fb + ".fig"); exportgraphics(f, fb + ".png", 'Resolution', 150); close(f);
    end
end

%% Needle vs thesis source (2D PG-TV, far field) for both protocols
cfg = AcqSim_Config();
fn = fieldnames(AllT);
for proto = ["MB" "SBM"]
    k = fn(endsWith(fn, "_" + proto));
    if numel(k) < 2, continue; end
    sweepPattern = 8; if proto == "SBM", sweepPattern = 21; end
    f = figure('Color','w','Position',[100 100 900 520]); ax = axes(f); hold(ax,'on');
    C = lines(numel(k));
    for i = 1:numel(k)
        T = AllT.(k{i}); se = sortrows(T(T.Pattern == sweepPattern & T.Estimator == "2DPG_TV", :), 'f_Hz');
        errorbar(ax, se.f_Hz/1e3, se.Far_mean, se.Far_std, '-o', 'Color', C(i,:), 'LineWidth', 2.5, ...
            'MarkerSize', 8, 'DisplayName', strrep(strrep(erase(k{i}, "_" + proto), '_', ' '), 'p', '.'));
    end
    yline(ax, cfg.c_ref, 'k--', 'LineWidth', 1.5, 'DisplayName', 'c_{ref}');
    xlabel(ax, 'Excitation frequency [kHz]', 'FontWeight', 'bold'); ylabel(ax, 'Phase speed [m/s]', 'FontWeight', 'bold');
    title(ax, sprintf('%s - 2D PG-TV far-field speed: new needle vs thesis source', proto), 'FontWeight', 'bold');
    legend(ax, 'Location', 'best'); grid(ax, 'on'); set(ax, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 1.5);
    fb = fullfile(cfg.results_path, "FWHM_Comparison_" + proto);
    savefig(f, fb + ".fig"); exportgraphics(f, fb + ".png", 'Resolution', 150); close(f);
end
%% Edge (rim) check: far field r >= R vs r >= R + 1 mm (pattern of the sweep)
% With the near field blanked BEFORE the estimators, the Delaunay
% interpolation of MB bridges the hole and the 2D PG kernels at the edge
% read high speeds (phase-flat bridge). Excluding one extra kernel width
% (~1 mm) shows how much of the 2D far-field mean is this edge effect.
rows = {};
for set = string({sets.name})
    C = load(fullfile(cfg.results_path, set, "AcqSim_Config_used.mat"), "cfg"); R0 = C.cfg.NeedleRadius_mm;
    for proto = ["MB" "SBM"]
        if proto == "MB"
            f = dir(fullfile(cfg.results_path, set, "Multimeridian_MB", "**", "MB_SimResults_8_*.mat")); estN = ["MB_1DPG" "MB_2DPG" "MB_2DPG_TV"];
        else
            f = dir(fullfile(cfg.results_path, set, "SBM", "37_Ascans_21_*", "**", "SBM_SimResults_*.mat")); estN = ["SBM_2DPG" "SBM_2DPG_TV"];
        end
        for i = 1:numel(f)
            S = load(fullfile(f(i).folder, f(i).name), 'Results'); Tus = round(S.Results.Parameters.T_excitation*1e6);
            for e = estN
                M = S.Results.(e); E = M.Elastogram_Masked; x = M.Xaxis_mm(:).'; y = M.Yaxis_mm(:);
                if isequal(size(E), [numel(y) numel(x)]), [X, Y] = meshgrid(x, y); else, [X, Y] = ndgrid(x, y); end
                r = hypot(X - cfg.Source_mm(1), Y - cfg.Source_mm(2)); inB = r <= cfg.FOV_mm/2 - cfg.MarginPx*cfg.FOV_mm/cfg.nPoints;   % thesis background disk (same outer edge as the ROIs)
                v1 = E(r >= R0 & inB & isfinite(E)); v2 = E(r >= R0 + 1 & inB & isfinite(E));
                rows(end+1, :) = {set, proto, Tus, erase(e, ["MB_" "SBM_"]), R0, mean(v1), std(v1), median(v1), mean(v2), std(v2), median(v2)}; %#ok<AGROW>
            end
        end
    end
end
Rim = sortrows(cell2table(rows, 'VariableNames', {'Set','Protocol','Period_us','Estimator','R_mm', ...
    'Mean_R','Std_R','Median_R','Mean_Rplus1','Std_Rplus1','Median_Rplus1'}), {'Set','Protocol','Estimator','Period_us'});
writetable(Rim, fullfile(cfg.results_path, "RimCheck_FarField_R_vs_Rplus1mm.xlsx"));

save(fullfile(cfg.results_path, "AcqSim_AllSummaries.mat"), "AllT");
fprintf('Analysis finished.\n');

%% ------------------------------------------------------------------------
function s = roiName(roi, cfg)
switch roi
    case "Far", s = sprintf("far field r >= R + 1 mm = %.2f mm (edge-free)", cfg.NeedleRadius_mm + cfg.EdgeMargin_mm);
    case "Ann", s = sprintf("r >= %.2f mm (= near-field blanking radius)", cfg.NeedleExclusion_mm);
    otherwise,  s = "thesis background disk";
end
end

function [slope, R2, p] = linfit(x, y)
x = x(:); y = y(:);
X = [ones(size(x)) x]; b = X\y; yh = X*b;
SSr = sum((y-yh).^2); SSt = sum((y-mean(y)).^2);
R2 = 1 - SSr/SSt; slope = b(2);
n = numel(x); se = sqrt(SSr/(n-2)/sum((x-mean(x)).^2));
tstat = slope/se; p = 2*(1 - tcdf(abs(tstat), n-2));
end

function NF = nearFieldRadius(T, cfg)
% Radial profile of the 2D PG on the full wave field (stored in every run
% as Results.OG_WaveProp); one profile per period.
periods = unique(T.Period_us);
rows = zeros(numel(periods), 4);
for i = 1:numel(periods)
    file = T.File(find(T.Period_us == periods(i), 1));
    hit = dir(fullfile(cfg.results_path, "**", file));
    S = load(fullfile(hit(1).folder, hit(1).name), 'Results');
    OG = S.Results.OG_WaveProp;
    [X, Y] = ndgrid(OG.Xaxis_mm, OG.Yaxis_mm); r = hypot(X, Y);
    edges = 0:0.25:max(OG.Xaxis_mm);
    prof = nan(numel(edges)-1, 1);
    for b = 1:numel(edges)-1
        m = r >= edges(b) & r < edges(b+1);
        prof(b) = median(OG.Elastogram_Masked(m), 'omitnan');
    end
    rc = edges(1:end-1) + 0.125;
    bad = abs(prof - cfg.c_ref)/cfg.c_ref > 0.05 | isnan(prof);
    valid = rc <= 0.85*max(OG.Xaxis_mm);            % ignore the PML/edge ring
    lastBad = find(bad(:) & valid(:), 1, 'last');
    lastValid = find(valid, 1, 'last');
    if isempty(lastBad), R5 = 0;
    elseif lastBad >= lastValid - 1, R5 = NaN;      % never within 5 % up to 0.85*FOV/2 (global bias)
    else, R5 = rc(lastBad) + 0.125; end
    rows(i,:) = [periods(i), 1e6/periods(i), R5, cfg.c_ref/(1e6/periods(i))*1e3];
end
NF = array2table(rows, 'VariableNames', {'Period_us','f_Hz','R5_mm','Lambda_mm'});
end

function Eq = toCommonGrid(E, x, y, xc, yc)
% Resample a map (imagesc(x, y, E) convention, [Ny x Nx]) onto the common
% grid of the elastogram figure (nearest neighbour, NaN outside), so every
% panel of MMB_PlotElastogramGrid shares the same mm axes.
x = x(:).'; y = y(:);
if ~isequal(size(E), [numel(y) numel(x)]), E = E.'; end
[Xq, Yq] = meshgrid(xc, yc);
Eq = interp2(x, y, double(E), Xq, Yq, 'nearest', NaN);
end
