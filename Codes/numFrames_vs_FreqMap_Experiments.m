function numFrames_vs_FreqMap_Experiments(Periods_us, Ts_OCT_us, SaveBase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 05-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version of numFrames_vs_FreqMap for the excitation periods used in the
% validation experiments (200, 250, 285 and 500 us).
%
% For the Synchronized BM-mode (SBM) protocol, the delay between
% consecutive BM frames must be an integer multiple of the OCT sampling
% period. For an excitation period T = M*T_OCT, the feasible numbers of
% synchronized frames are the divisors of M.
%
% Inputs (optional):
%   Periods_us : excitation periods [us]   (default [200 250 285 500])
%   Ts_OCT_us  : OCT sampling period [us]  (default 5, i.e. 200 kHz)
%   SaveBase   : path without extension to save .png and .fig (default: none)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(Periods_us), Periods_us = [200 250 285 500]; end
if nargin < 2 || isempty(Ts_OCT_us),  Ts_OCT_us = 5; end
if nargin < 3, SaveBase = ''; end

colors = lines(numel(Periods_us));

fig = figure('Color', 'w', 'Position', [100 100 1000 650]);
hold on;
h = gobjects(numel(Periods_us), 1);
legText = cell(numel(Periods_us), 1);

for k = 1:numel(Periods_us)

    M = Periods_us(k)/Ts_OCT_us;
    if abs(M - round(M)) > 1e-9
        warning('%g us is not an integer multiple of T_OCT = %g us; skipped.', Periods_us(k), Ts_OCT_us);
        continue;
    end

    N_frames = divisors(round(M));
    f = 1/(Periods_us(k)*1e-6)*ones(size(N_frames));

    h(k) = plot(f, N_frames, 'o', 'MarkerSize', 7, 'LineWidth', 1.5, ...
        'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', colors(k,:));

    % Inline labels for the larger counts (small ones are listed in the legend)
    for n = N_frames(N_frames >= 8)
        text(f(1) + 60, n, sprintf('%d', n), 'FontSize', 11, 'Color', colors(k,:));
    end
    legText{k} = sprintf('%g \\mus (%.1f Hz): N = %s', Periods_us(k), 1e6/Periods_us(k), ...
        strjoin(string(N_frames), ', '));
end

hold off; grid on; box on;
xlabel('Frequency [Hz]');
ylabel('Number of Frames');
title('Frequency vs Number of Frames for the Experimental Excitation Periods');
valid = isgraphics(h);
legend(h(valid), legText(valid), 'Location', 'northeast', 'FontSize', 12);
fmax = 1e6/min(Periods_us); fmin = 1e6/max(Periods_us);
xlim([fmin - 0.15*(fmax-fmin), fmax + 0.15*(fmax-fmin)]);
ylim([-5, max(Periods_us)/Ts_OCT_us + 10]);
set(gca, 'FontSize', 14);

if ~isempty(SaveBase)
    savefig(fig, [SaveBase '.fig']);
    exportgraphics(fig, [SaveBase '.png'], 'Resolution', 200);
end

end
