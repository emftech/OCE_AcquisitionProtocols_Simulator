function R = AcqSim_NeedleRadius(FWHM_mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (batch layer written with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Near-field blanking radius agreed with the experimental processing
% sessions: R = max(1.2 mm, 2*sqrt(2)*sigma_s), i.e. the compact support of
% the k-Wave Gaussian source (G > exp(-4)), never below the 1.2 mm used for
% the old needle in the 06-10-2026 experiments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_s = FWHM_mm/(2*sqrt(2*log(2)));
R = max(1.2, 2*sqrt(2)*sigma_s);
end
