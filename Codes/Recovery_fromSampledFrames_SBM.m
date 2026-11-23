function [Acq_Phase_Matrix, Recovered_Frames] = Recovery_fromSampledFrames_SBM( ...
                            nTurnSamples, T_excitation,...
                            Ts_OCT, nFrames, SampledFrames_raw,phaseStepSamples)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 08-11-2025
% Date of last modification (DD/MM/YY): 08-11-2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Function's Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform the phase compensation in order to recover the original wave
% propagation video (unaffected by doppler shifts) using the timeframes for
% acquisition in order to compute the phase that's characteristic of the
% galvanometer's movement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   - nAscans (scalar, integer): A-lines per horizontal line
%   - nBscans (scalar, integer): number of horizontal lines -> must be an 
%                                odd number to go through the middle
%   - nFrames (scalar, integer): Number of BM Frames
%   - startDelay_samples (scalar, integer): small delay so excitation rising 
%                                           edge leads acquisition/scan
%   - nTurnSamples (scalar, integer): number of samples for each transition
%                                     turn between B-scans
%   - AcqPulseWidth (scalar, integer): acquisition trigger's 'TTL' 
%                                      pulse-width [samples]
%   - ExcPulseWidth (scalar, integer): excitation trigger's 'TTL' 
%                                      pulse-width [samples]
%   - T_excitation (scalar, double): Excitation waveform's temporal period [s]
%   - Ts_OCT (scalar, double): OCT's sampling period, used to check if the
%                              excitation's period is a multiple of the
%                              OCT's universal clock [s]
%   - SampledFrames_raw (matrix, double): 
% Outputs:
%   - Acq_Phase_Matrix
%   - Recovered_Frames (matrix, double): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute matrix of phases to compensate

% Phase stepping: how many OCT samples to shift acquisition per frame
% phaseStepSamples = (T_excitation/nFrames)/Ts_OCT;

% Compute spatial samples' distribution:
[nBscans, nAscans, nFrames] = size(SampledFrames_raw);

% Total samples for which the galvo travels in a single frame
N_singleframe = nAscans*nBscans + (nBscans-1)*nTurnSamples;

f_excitation = 1/T_excitation;
omega_excitation = 2*pi*f_excitation;

% Single frame timeline:
timeline_singleframe = Ts_OCT:Ts_OCT:(N_singleframe)*Ts_OCT;

% Compute mask that discards the turn samples:
mask_isAline = false(1, N_singleframe);   % initialize all false

if nTurnSamples == 0
    mask_isAline(:) = true;
else
    stride = nAscans + nTurnSamples;      % distance between the *start* of successive rows
    % compute start indices of each row (1-based)
    rowStarts = 1 + (0:(nBscans-1)) * stride;   % length nBscans

    % produce all a-line indices per row
    % use kron to broadcast the offsets [0:(nAscans-1)]
    offsets = 0:(nAscans-1);
    allIdx = rowStarts(:) + offsets(:).';   % produces nBscans x nAscans matrix (implicit broadcasting in recent MATLAB)
    allIdx = allIdx.';                      % nAscans x nBscans -> transpose to vectorize row-major
    allIdx = allIdx(:).';                   % flatten in row-major order

    % keep only those indices falling inside N_singleframe (last row may produce indices beyond end)
    allIdx = allIdx(allIdx <= N_singleframe);

    mask_isAline(allIdx) = true;

end

% allocate
Acq_Phase_Matrix = NaN(nBscans, nAscans);

% prepare L->R offsets and R->L flips
offsets_lr = [0:nAscans/2-1 floor(nAscans/2-1) floor(nAscans/2):nAscans-2];          % 1 x nAscans
offsets_rl = fliplr(offsets_lr);    % 1 x nAscans (for serpentine rows)

% make matrix of offsets per row (nBscans x nAscans)
rowIdx = (1:nBscans).';
% <-- CHANGE: flip the opposite set of rows compared to previous version.
% Previously you flipped even rows; here we flip odd rows so serpentine matches
% the acquisition ordering from your diagnostic. If you prefer flipping even
% rows instead, revert this line to mod(rowIdx,2)==0.
isFlipRow = mod(rowIdx,2) == 1;     % true for odd-numbered rows (will apply fliplr)

offsets_matrix = repmat(offsets_lr, nBscans, 1);
if any(isFlipRow)
    offsets_matrix(isFlipRow, :) = repmat(offsets_rl, sum(isFlipRow), 1);
end

% compute candidate indices into the single-frame timeline
% idx_mat(r,c) = rowStarts(r) + offsets_matrix(r,c)
idx_mat = rowStarts(:) + offsets_matrix;   % nBscans x nAscans

% Build a validity mask that requires:
%  - index within [1, N_singleframe]
%  - mask_isAline at that index is true (not a turn sample)
valid_idx_mask = (idx_mat >= 1) & (idx_mat <= N_singleframe);   % bounds
if any(valid_idx_mask(:))
    % for indices within bounds, check mask_isAline
    % create logical of same shape: true if mask_isAline(index)==true
    linearIdx = idx_mat(valid_idx_mask);         % vector of candidate indices
    valid_mask_after_turns = mask_isAline(linearIdx); % see which are real a-lines

    % insert into the full valid matrix
    tmp = false(size(idx_mat));
    tmp(valid_idx_mask) = valid_mask_after_turns;   % now tmp marks final valid A-lines
    valid_idx_mask = tmp;                           % nBscans x nAscans logical
end

% fill phase matrix using the final valid_idx_mask
if any(valid_idx_mask(:))
    idx_vals = idx_mat(valid_idx_mask);            % the actual timeline indices used
    t_vals = timeline_singleframe(idx_vals);       % timestamps (seconds)
    Acq_Phase_Matrix(valid_idx_mask) = omega_excitation * t_vals;  % continuous phase
end

% Wrapped version (principal value (-pi,pi])
Acq_Phase_Matrix = angle( exp(1j * Acq_Phase_Matrix) );

%% Compute FFT along time dimension and compensate the phase

Nfft = nFrames;
% For the acquired OCT frames, frames are separated by phaseStepSamples * Ts_OCT
dt_SampledFrames = double(phaseStepSamples) * Ts_OCT;
% Compute FFT along frames for the SampledFrames_raw matrix
F_SampledFrames  = fft(double(SampledFrames_raw), Nfft, 3);
% Compute frequency bins
freqs_SampledFrames = (0:(Nfft-1))' * (1/dt_SampledFrames / Nfft);
% Get exact frequency bin at the excitiation frequency
[~, kpos_SampledFrames] = min(abs(freqs_SampledFrames - f_excitation));
% Extract phasor at excitation frequency bin
Ph_SampledFrames = squeeze(F_SampledFrames(:,:,kpos_SampledFrames)); % nB x nA complex

% Compensate for the acquisition phase
Ph_Compensated = Ph_SampledFrames .* exp(1j * Acq_Phase_Matrix);

% Replace in the FFT the compensated phasor:
F_compensated = F_SampledFrames;
F_compensated(:,:,kpos_SampledFrames) = Ph_Compensated;

% compute correct conjugate index (robust)
Nfft = size(F_compensated, 3);
kpos = kpos_SampledFrames;
kneg = mod(Nfft - kpos + 2 - 1, Nfft) + 1;   % robust 1-based mapping

% If kpos is not a DC/Nyquist special case then set conjugate to conj()
isDC = (kpos == 1);
isNyquist = (mod(Nfft,2)==0) && (kpos == (Nfft/2 + 1));
if ~isDC && ~isNyquist
    F_compensated(:,:,kneg) = conj(Ph_Compensated);
else
    % DC or Nyquist: nothing to mirror, keep single real-valued bin as needed
end

% Take the inverse FFT in the third dimension in order to compute the Recovered_Frames
% volume:
Recovered_Frames = ifft(F_compensated, [], 3, 'symmetric');


end