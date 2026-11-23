function AcqSim_SaveElastogram(E, x_mm, y_mm, clim_ms, ttl, filebase, Source_mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masked elastogram (turbo, masked pixels in the darkest colour), saved as
% .fig + .png, bold thesis style. Optional Source_mm = [x y]: white '+' at
% the needle tip (experiment convention for near-field-blanked maps).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 560 500]);
ax = axes(f);
im = imagesc(ax, x_mm, y_mm, E);
set(im, 'AlphaData', isfinite(E));
set(ax, 'Color', [0.19 0.07 0.23]);           % turbo(1): darkest colour for masked pixels
axis(ax, 'image'); set(ax, 'YDir', 'normal');
colormap(ax, turbo); clim(ax, clim_ms);
cb = colorbar(ax); cb.Label.String = 'Phase speed [m/s]'; cb.Label.FontWeight = 'bold';
xlabel(ax, 'x [mm]', 'FontWeight', 'bold'); ylabel(ax, 'y [mm]', 'FontWeight', 'bold');
if nargin >= 7 && ~isempty(Source_mm)
    hold(ax, 'on'); plot(ax, Source_mm(1), Source_mm(2), '+', 'Color', 'w', 'MarkerSize', 12, 'LineWidth', 2); hold(ax, 'off');
end
title(ax, ttl, 'FontWeight', 'bold');
set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
savefig(f, filebase + ".fig");
exportgraphics(f, filebase + ".png", 'Resolution', 150);
close(f);
end
