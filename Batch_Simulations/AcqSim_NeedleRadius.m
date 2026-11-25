function R = AcqSim_NeedleRadius(FWHM_mm, floor_mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Near-field blanking radius agreed with the experimental processing
% sessions: R = max(floor_mm, 2*sqrt(2)*sigma_s), i.e. the compact support
% of the k-Wave Gaussian source (G > exp(-4)), never below the experimental
% blanking radius floor_mm.
%   floor_mm = 1.5 mm : thin 0.65 mm needle (14-11 / 25-11-2026 processing)
%   floor_mm = 1.2 mm : default, old needle of the 06-10-2026 experiments
%                       (legacy 18 x 18 mm batch of 20-10-2026)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2 || isempty(floor_mm), floor_mm = 1.2; end
sigma_s = FWHM_mm/(2*sqrt(2*log(2)));
R = max(floor_mm, 2*sqrt(2)*sigma_s);
end
