function numFrames_vs_FreqMap()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 15-10-2025
% Date of last modification (DD/MM/YY): 15-10-2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Function's Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function 'wave_prop_cylindrical_WE' generates the simulation of a
% Propagating Wave using the Monochromatic Spherical Solutions to the
% Helmholtz Equation, as described in the main live script. This function 
% uses directly a numerical implementation of the aforementioned solutions 
% to the Wave Equation (hence the name 'WE' at the end) and it does so for
% the harmonic, one frequency case (i.e. monochromatic solutions), effectively
% implementing the particle velocity field.
% Since it builds the most simple kind of propagation (cylindrical wavefront,
% without inclusions), the medium's sound speed is a scalar constant (c).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: None
% Outputs: Implicit - Figure with Number of Frames vs Frequency Graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation fo Number of Frames vs Frequency Graph

N_frames_250us = divisors(250/5);
N_frames_260us = divisors(260/5);
N_frames_270us = divisors(270/5);
N_frames_285us = divisors(285/5);
N_frames_300us = divisors(300/5);

f_250us = 1/250e-6*ones(size(N_frames_250us));
f_260us = 1/260e-6*ones(size(N_frames_260us));
f_270us = 1/270e-6*ones(size(N_frames_270us));
f_285us = 1/285e-6*ones(size(N_frames_285us));
f_300us = 1/300e-6*ones(size(N_frames_300us));

figure,
plot(f_285us,N_frames_285us, 'bo','MarkerSize',5)
hold on
plot(f_250us, N_frames_250us, 'ro','MarkerSize',5)
plot(f_260us, N_frames_260us, 'mo','MarkerSize',5)
plot(f_270us, N_frames_270us, 'ko','MarkerSize',5)
plot(f_300us, N_frames_300us, 'go','MarkerSize',5)
xlabel('Frequency [Hz]')
ylabel('Number of Frames')
title('Frequency vs Number of Frames for Different Excitation Periods')
legend('285 µs', '250 µs', '260 µs', '270 µs', '300 µs')
hold off
grid on;
xlim([3200 4200])
ylim([-10 80])

end