%% Make_SourceProfile_Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Needle-width measurement on the crosshair B-mode of 20-10-2026
%    (10 x 10 mm scan, 300 A-scans per B-scan = one meridian; the axis of
%    the photo is perspective-skewed, so A-scan positions are mapped with
%    the local yellow cursor lines at A-scans 0, 150 and 300):
%    depth band rows 380-470 of the photo, shadow edges at photo columns
%    475 and 788 -> A-scans 89.5 and 195.8 -> 106.3 A-scans -> 3.54 mm.
% 2) Spatial (Gaussian, FWHM 3.54 vs 0.30 mm) and temporal (QH 3 cycles,
%    harmonic causal) excitation profiles used by the batch.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all force;
addpath(fileparts(mfilename("fullpath")));
cfg = AcqSim_Config(); AcqSim_SetupPaths(cfg);
out = cfg.results_path + "/Source_Profile";
I = imread(out + "/Needle_Crosshair_Bmode_10x10mm_20-10-2026.png"); I = I(:,:,1:3);
G = double(rgb2gray(I));

% Cursor columns (A-scan 0, 150, 300) at the analysed depth band
c0 = 197; c150 = 663; c300 = 1072; rows = 380:470;
prof = movmean(mean(G(rows,:),1), 9);
thr = 0.5*(median(prof(400:460)) + median(prof(500:760)));      % half-way between phantom and shadow
left  = find(prof(430:560) < thr, 1, 'first') + 429;
right = find(prof(700:830) > thr, 1, 'first') + 699;
toA = @(c) (c <= c150).*((c - c0)/(c150 - c0)*150) + (c > c150).*(150 + (c - c150)/(c300 - c150)*150);
aL = toA(left); aR = toA(right); widthA = aR - aL; width_mm = widthA*10/300;
fprintf('Shadow edges: cols %d / %d -> A-scans %.1f / %.1f -> %.1f A-scans = %.2f mm\n', left, right, aL, aR, widthA, width_mm);

f = figure('Color','w','Position',[50 50 1500 620]);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); imshow(I, 'Parent', ax); hold(ax, 'on');
plot(ax, [left left], [120 1150], 'c-', 'LineWidth', 2.5); plot(ax, [right right], [120 1150], 'c-', 'LineWidth', 2.5);
rectangle(ax, 'Position', [1 rows(1) size(I,2)-1 numel(rows)], 'EdgeColor', 'm', 'LineWidth', 1.5, 'LineStyle', '--');
text(ax, (left+right)/2, 100, sprintf('%.1f A-scans = %.2f mm', widthA, width_mm), 'Color', 'c', ...
    'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'center', 'BackgroundColor', 'k');
title(ax, {'Crosshair B-mode, 10 x 10 mm, 300 A-scans per B-scan', 'needle shadow edges (cyan), analysed depth band (magenta)'}, 'FontWeight', 'bold');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, toA(1:size(I,2)), prof, 'k-', 'LineWidth', 2);
xline(ax, [aL aR], 'c-', 'LineWidth', 2.5); yline(ax, thr, 'm--', 'LineWidth', 1.5);
xlim(ax, [0 300]); xlabel(ax, 'A-scan # (300 = 10 mm)', 'FontWeight', 'bold'); ylabel(ax, 'Mean intensity (photo, a.u.)', 'FontWeight', 'bold');
title(ax, 'Lateral intensity profile in the depth band', 'FontWeight', 'bold'); grid(ax, 'on');
set(ax, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 1.5);
exportgraphics(f, out + "/Needle_Width_Measurement.png", 'Resolution', 150); savefig(f, out + "/Needle_Width_Measurement.fig"); close(f);

% Spatial and temporal profiles
x = linspace(-9, 9, 1001);
f = figure('Color','w','Position',[50 50 1500 430]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
for F = [cfg.SourceFWHM_mm 0.30]
    s = F/(2*sqrt(2*log(2))); g = exp(-x.^2/(2*s^2)); g(g < exp(-4)) = 0;
    plot(ax, x, g, 'LineWidth', 2.5, 'DisplayName', sprintf('FWHM %.2f mm (support r < %.2f mm)', F, 2*sqrt(2)*s));
end
xlabel(ax, 'r [mm]', 'FontWeight', 'bold'); ylabel(ax, 'G_{xy}', 'FontWeight', 'bold'); legend(ax, 'Location', 'southoutside');
title(ax, 'Spatial source profile', 'FontWeight', 'bold'); grid(ax, 'on'); set(ax, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 1.5);
ax = nexttile(tl); hold(ax, 'on');
k = linspace(0, 6000, 500); fr = [2000 2941 3509 4545];
for F = [cfg.SourceFWHM_mm 0.30]
    s = F/(2*sqrt(2*log(2)))*1e-3;
    plot(ax, k, 20*log10(exp(-k.^2*s^2/2)), 'LineWidth', 2.5, 'DisplayName', sprintf('FWHM %.2f mm', F));
end
for q = fr, xline(ax, 2*pi*q/cfg.c_ref, ':', sprintf('%.1f kHz', q/1e3), 'LineWidth', 1.5, 'HandleVisibility', 'off', 'FontWeight', 'bold'); end
ylim(ax, [-60 3]); xlabel(ax, 'Wavenumber k [rad/m]', 'FontWeight', 'bold'); ylabel(ax, 'Source k-spectrum [dB]', 'FontWeight', 'bold');
legend(ax, 'Location', 'southwest'); title(ax, 'Radiation weight at k_0 = \omega/c_{ref}', 'FontWeight', 'bold'); grid(ax, 'on');
set(ax, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 1.5);
ax = nexttile(tl); hold(ax, 'on');
T = 285e-6; t = (0:499)*cfg.Ts_OCT;
qh = sin(2*pi*t/T).*(t < 3*T); hc = sin(2*pi*t/T);
plot(ax, t*1e3, qh + 2.5, 'LineWidth', 2, 'DisplayName', 'MB: QH, 3 cycles per 450+50 ensemble');
plot(ax, t*1e3, hc, 'LineWidth', 2, 'DisplayName', 'SBM: harmonic, causal (steady state)');
xlabel(ax, 't [ms]', 'FontWeight', 'bold'); yticks(ax, []); legend(ax, 'Location', 'southoutside');
title(ax, 'Temporal profiles (285 \mus shown)', 'FontWeight', 'bold'); grid(ax, 'on'); set(ax, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 1.5);
exportgraphics(f, out + "/Excitation_Profiles.png", 'Resolution', 150); savefig(f, out + "/Excitation_Profiles.fig"); close(f);
