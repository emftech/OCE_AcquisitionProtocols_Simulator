function AcqSim_WriteVideo(V, x_mm, y_mm, ttl, filebase, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcqSim_WriteVideo
%
% Writes the frame videos that the live scripts displayed (sampled,
% interpolated, raw, recovered and wave-propagation frames) to .mp4, using
% the wave-display convention of 20-10-2026: blue-white-red (MMB_RedBlue)
% with white = 0, symmetric limits at +-P97 of |v|, NaN / unsampled pixels
% in gray [0.6 0.6 0.6].
%
%   V        : [N1 x N2 x Nt]; displayed as imagesc(x_mm, y_mm, frame)
%              after the optional transpose (opts.Transpose = true when
%              V is indexed (x, y, t)). Pass the UNBLANKED data.
%   ttl      : title prefix; the frame number / time is appended.
%   opts     : .Transpose (false), .YDir ('normal'), .FrameStep (1),
%              .FrameRate (15), .dt_s ([] -> frame number in title),
%              .NaNMask ([] or [N1 x N2] logical, V orientation): near-field
%              disk / band shown gray in every frame,
%              .CLim ([]): imposed symmetric limit (overrides P97), e.g. the
%              limit of a cropped video when writing its full-field twin,
%              .WriteFull (true): when NaNMask is given, also write
%              <filebase>_Full.mp4 with the near field kept, using the SAME
%              colour limits as the cropped video (P97 of the far field),
%              so the near field saturates and the far field stays readable.
%
% Outputs: <filebase>.mp4 (cropped when NaNMask is given) and, optionally,
% <filebase>_Full.mp4.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 6, opts = struct(); end
if ~isfield(opts, 'Transpose'), opts.Transpose = false; end
if ~isfield(opts, 'YDir'),      opts.YDir = 'normal'; end
if ~isfield(opts, 'FrameStep'), opts.FrameStep = 1; end
if ~isfield(opts, 'FrameRate'), opts.FrameRate = 15; end
if ~isfield(opts, 'dt_s'),      opts.dt_s = []; end
if ~isfield(opts, 'NaNMask'),   opts.NaNMask = []; end
if ~isfield(opts, 'WriteFull'), opts.WriteFull = true; end
if ~isfield(opts, 'CLim'),      opts.CLim = []; end

V = double(V);   % local copy: the caller's array is never modified
if isempty(opts.NaNMask)
    lim = climOf(V);
    if ~isempty(opts.CLim), lim = opts.CLim; end
    if isempty(lim), return; end
    writeOne(V, lim, ttl, filebase + ".mp4", x_mm, y_mm, opts);
    return;
end

Vc = V;
Vc(repmat(logical(opts.NaNMask), 1, 1, size(V,3))) = NaN;
lim = climOf(Vc);                         % far-field (cropped) limits for both videos
if isempty(lim), return; end
writeOne(Vc, lim, ttl, filebase + ".mp4", x_mm, y_mm, opts);
if opts.WriteFull
    ttlFull = regexprep(char(ttl), '\s*\(r < [^)]*\)', '');
    writeOne(V, lim, [ttlFull ' - full field (near field saturated)'], filebase + "_Full.mp4", x_mm, y_mm, opts);
end
end

%% ------------------------------------------------------------------------
function lim = climOf(V)
a = abs(V(isfinite(V)));
if isempty(a), lim = []; return; end
lim = prctile(a, 97); if lim <= 0, lim = max(a); end
end

function writeOne(V, lim, ttl, file, x_mm, y_mm, opts)
vw = VideoWriter(char(file), 'MPEG-4');
vw.FrameRate = opts.FrameRate; vw.Quality = 95;
open(vw);
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 620 540]);
ax = axes(f);
for k = 1:opts.FrameStep:size(V,3)
    F = V(:,:,k);
    if opts.Transpose, F = F.'; end
    cla(ax);
    im = imagesc(ax, x_mm, y_mm, F);
    set(im, 'AlphaData', isfinite(F));
    set(ax, 'Color', [0.6 0.6 0.6], 'YDir', opts.YDir);
    axis(ax, 'image'); colormap(ax, MMB_RedBlue(255)); clim(ax, [-lim lim]);
    cb = colorbar(ax); cb.Label.String = 'Particle velocity (a.u.)';
    xlabel(ax, 'x [mm]', 'FontWeight', 'bold'); ylabel(ax, 'y [mm]', 'FontWeight', 'bold');
    if isempty(opts.dt_s)
        title(ax, {ttl, sprintf('Frame %d / %d', k, size(V,3))}, 'FontWeight', 'bold');
    else
        title(ax, {ttl, sprintf('t = %.3f ms', (k-1)*opts.dt_s*1e3)}, 'FontWeight', 'bold');
    end
    set(ax, 'FontSize', 11, 'FontWeight', 'bold');
    drawnow;
    writeVideo(vw, getframe(f));
end
close(vw); close(f);
end
