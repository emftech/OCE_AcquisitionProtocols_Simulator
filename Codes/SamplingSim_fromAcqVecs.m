function [SampledFrames_raw, SampledFrames_vis] = SamplingSim_fromAcqVecs(...
    WavePropVid, x_grid, y_grid, t_grid, ...
    GalvoX, GalvoY, AcqTrigger, Exc_Aux_Trigger, ...
    nFrames, nAscans, nBscans, nTurnSamples, T_exc, Ts_OCT)
% SamplingSim_fromAcqVecs - sample WavePropVid obeying acquisition vectors
% Outputs:
%   SampledFrames_raw (nBscans x nAscans x nFrames)
%   SampledFrames_vis (same but rows scanned R->L are flipped for visualization)
%
% NOTE: This version maps each acquired sample's time to the *phase relative
% to that frame's excitation rising edge* (t_rel = t_sample - t_exc_start).
% That makes phaseStepSamples == 0 produce static frames.

%% Input normalization / checks
GalvoX = GalvoX(:); GalvoY = GalvoY(:);
AcqTrigger = double(AcqTrigger(:));
Exc_Aux_Trigger = double(Exc_Aux_Trigger(:));

N_total = numel(GalvoX);
if numel(GalvoY) ~= N_total || numel(AcqTrigger) ~= N_total || numel(Exc_Aux_Trigger) ~= N_total
    error('GalvoX/GalvoY/AcqTrigger/Exc_Aux_Trigger must have same length.');
end

if nargin < 14 || isempty(Ts_OCT)
    error('Ts_OCT (OCT sampling period [s]) is required as the last argument.');
end
dt_sample = Ts_OCT;

% gridded interpolant
F_wave = griddedInterpolant({x_grid, y_grid, t_grid}, WavePropVid, 'cubic', 'linear');

% preallocate
SampledFrames_raw = nan(nBscans, nAscans, nFrames);
SampledFrames_vis = nan(nBscans, nAscans, nFrames);

% detect excitation rising edges (frame starts)
exc_rising = find(diff([0; double(Exc_Aux_Trigger)]) == 1);

% global sample times (s)
t_sample_global = (0:(N_total-1))' * dt_sample;

% dynamic tolerance for plateau detection in Y
tol_abs = 1e-9;
tol_rel = 1e-6;
tol = tol_abs + tol_rel * max(1, max(abs(GalvoY)));

% loop frames
for m = 1:nFrames
    % determine global start/end indices for this BM frame using exc_rising (if present)
    if m <= numel(exc_rising)
        start_idx = exc_rising(m);
        if m < numel(exc_rising)
            end_idx = exc_rising(m+1) - 1;
        else
            end_idx = N_total;
        end
    else
        if ~isempty(exc_rising)
            start_idx = exc_rising(end);
            end_idx = N_total;
        else
            start_idx = 1;
            end_idx = N_total;
        end
    end
    start_idx = max(1, start_idx);
    end_idx = min(N_total, end_idx);
    if end_idx < start_idx
        continue;
    end

    idx_range = (start_idx:end_idx)';
    Xseg = GalvoX(idx_range);
    Yseg = GalvoY(idx_range);

    % plateau detection (where Y is approximately constant)
    dY = [0; abs(diff(Yseg))];
    is_plateau = dY < tol;
    run_starts = find(diff([0; is_plateau]) == 1);
    run_ends   = find(diff([is_plateau; 0]) == -1);
    if numel(run_starts) ~= numel(run_ends)
        nmin = min(numel(run_starts), numel(run_ends));
        run_starts = run_starts(1:nmin);
        run_ends = run_ends(1:nmin);
    end
    run_lengths = run_ends - run_starts + 1;

    % choose good runs
    minLen = max(1, round(nAscans * 0.5));
    good_mask = run_lengths >= minLen;
    good_starts = run_starts(good_mask);
    good_ends   = run_ends(good_mask);

    % relax if not enough runs
    thr = minLen - 1;
    while numel(good_starts) < nBscans && thr >= 1
        good_mask = run_lengths >= thr;
        good_starts = run_starts(good_mask);
        good_ends   = run_ends(good_mask);
        thr = thr - 1;
    end

    if numel(good_starts) < nBscans
        [~, idx_sorted] = sort(run_lengths, 'descend');
        pick = idx_sorted(1:min(nBscans, numel(idx_sorted)));
        pick = sort(pick);
        good_starts = run_starts(pick);
        good_ends = run_ends(pick);
    end

    nTaken = min(nBscans, numel(good_starts));
    good_starts = good_starts(1:nTaken);
    good_ends   = good_ends(1:nTaken);

    % prepare queries
    N_needed = nAscans * nBscans;
    Xq_all = nan(N_needed,1);
    Yq_all = nan(N_needed,1);
    Tq_all = nan(N_needed,1);
    writePos = 1;

    for row = 1:nBscans
        if row <= nTaken
            rs = good_starts(row);  % these are indices relative to Yseg / idx_range
            re = good_ends(row);
            % ------------------ NEW ROBUST SELECTION -----------------------
            % Select nAscans unique integer indices inside idx_range(rel) using:
            %  - if run length L >= nAscans : pick floor-based evenly spaced indices
            %  - if L < nAscans  : take all inside run and expand left/right inside frame
            L = re - rs + 1;
            if L >= nAscans
                % evenly spaced integer positions inside [rs..re]
                if nAscans == 1
                    rel_pos = rs;
                else
                    % positions in [1..L], evenly spaced, integer indices chosen by:
                    % floor( (0:(nAscans-1)) * (L-1)/(nAscans-1) ) + 1
                    rel_pos = floor((0:(nAscans-1)) * (L-1) / max(1,(nAscans-1))) + 1;
                    % convert to indices in idx_range coordinates
                    rel_pos = rs - 1 + rel_pos;  % now in [rs..re]
                end
                global_idxs = idx_range(rel_pos(:));
            else
                % run is shorter than requested -> take whole run and expand to neighbors
                avail = (rs:re)';                  % positions inside idx_range coords
                global_idxs = idx_range(avail);   % initial selection
                % expand alternately left/right inside idx_range window
                left = rs - 1;
                right = re + 1;
                % bounds for available positions in relative coordinates: 1 .. length(idx_range)
                left_bound = 1;
                right_bound = length(idx_range);
                while numel(global_idxs) < nAscans && (left >= left_bound || right <= right_bound)
                    if left >= left_bound
                        global_idxs = [idx_range(left); global_idxs]; %#ok<AGROW>
                        left = left - 1;
                        if numel(global_idxs) == nAscans
                            break;
                        end
                    end
                    if right <= right_bound
                        global_idxs = [global_idxs; idx_range(right)]; %#ok<AGROW>
                        right = right + 1;
                    end
                end
                % if still short, pad with NaN (should be rare)
                if numel(global_idxs) < nAscans
                    global_idxs = [global_idxs; nan(nAscans - numel(global_idxs),1)];
                end
                % finally trim/ensure length
                global_idxs = global_idxs(1:nAscans);
            end
            % ------------------ END NEW SELECTION --------------------------

            % unit conversion if x_grid small -> assume meters
            if max(abs(x_grid)) < 1
                Xq_all(writePos:writePos+nAscans-1) = GalvoX(global_idxs) * 1e-3;
                Yq_all(writePos:writePos+nAscans-1) = GalvoY(global_idxs) * 1e-3;
            else
                Xq_all(writePos:writePos+nAscans-1) = GalvoX(global_idxs);
                Yq_all(writePos:writePos+nAscans-1) = GalvoY(global_idxs);
            end
            % **absolute times** of those global indices
            Tq_all(writePos:writePos+nAscans-1) = (global_idxs - 1) * dt_sample;
        else
            % leave NaNs
        end
        writePos = writePos + nAscans;
    end

    % ----------------- IMPORTANT CHANGE -----------------
    % Map times into excitation-phase *relative to this frame's excitation start*.
    % This makes phaseStepSamples==0 produce static frames (same relative phases).
    t_exc_start = (start_idx - 1) * dt_sample;     % seconds, time of excitation rising edge for this frame
    t_rel = Tq_all - t_exc_start;                  % relative times (may be negative for pre-exc samples)
    % now map relative time into [0,T_exc)
    tq_phase = mod(t_rel, T_exc);
    % avoid exact zero-phase values (numerical)
    zero_idx = (tq_phase == 0) & (~isnan(Tq_all) & Tq_all > 0);
    tq_phase(zero_idx) = eps;
    % ----------------------------------------------------

    % sample interpolant
    valid_q = ~isnan(Xq_all) & ~isnan(Yq_all) & ~isnan(tq_phase);
    sampled_vals = nan(size(Xq_all));
    if any(valid_q)
        sampled_vals(valid_q) = F_wave(Yq_all(valid_q), Xq_all(valid_q), tq_phase(valid_q));
    end

    frame_img = reshape(sampled_vals, nAscans, nBscans).';
    
    SampledFrames_raw(:,:,m) = frame_img;

    % VIS version: flip rows scanned right->left
    frame_vis = frame_img;
    for row = 1:nBscans
        if row <= nTaken
            rs = good_starts(row);
            re = good_ends(row);
            % reuse robust selection to determine actual X_row ordering
            % compute same global_idxs as above logic for this row to inspect X ordering
            % (here keep it simple and use a safe approximation)
            rel_idxs = round(linspace(rs, re, nAscans));
            global_idxs = idx_range(rel_idxs);
            X_row = GalvoX(global_idxs);
            if mean(diff(X_row)) < 0
                frame_vis(row,:) = fliplr(frame_vis(row,:));
            end
        end
    end
    SampledFrames_vis(:,:,m) = frame_vis;
end

fprintf('Sampling done: produced SampledFrames_raw (%dx%dx%d) and SampledFrames_vis (%dx%dx%d)\n', ...
    size(SampledFrames_raw,1), size(SampledFrames_raw,2), size(SampledFrames_raw,3), ...
    size(SampledFrames_vis,1), size(SampledFrames_vis,2), size(SampledFrames_vis,3));

end