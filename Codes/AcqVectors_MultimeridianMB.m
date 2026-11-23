function [Exc_Trigger,Acq_Trigger,GalvoX,GalvoY] = AcqVectors_MultimeridianMB(horizontalMin,horizontalMax,verticalMin,...
                                                        verticalMax,nAscans,nBscans,nFrames,...
                                                        startDelay_samples,nSyncSamples,AcqPulseWidth,...
                                                        ExcPulseWidth,T_excitation,Ts_OCT,M_repetitions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 20-03-2026
% Date of last modification (DD/MM/YY): 20-03-2026
% Description:
% Generate acquisition vectors for the Multimeridian MB Mode scanning
% protocol, adapted to the simulation framework.
%
% Timing logic per spatial point:
%   - Excitation starts at the first sample of the point dwell
%   - Acquisition also starts at that same first sample
%   - Acquisition is repeated M_repetitions times (one A-line per OCT sample)
%   - Then there is an idle hold of startDelay_samples (kept for compatibility)
%   - Then nSyncSamples of idle smooth transition toward the next point
%
% Notes:
%   - If AcqPulseWidth > 1, the acquisition trigger may appear continuous
%     during the acquisition dwell, because acquisitions occur every sample.
%   - Set AcqPulseWidth = 1 if you want one-sample TTL pulses per A-line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    horizontalMin (1,1) double
    horizontalMax (1,1) double
    verticalMin   (1,1) double
    verticalMax   (1,1) double
    nAscans       (1,1) double {mustBeInteger,   mustBePositive}
    nBscans       (1,1) double {mustBeInteger,   mustBePositive}
    nFrames       (1,1) double {mustBeInteger,   mustBePositive}
    startDelay_samples (1,1) double {mustBeInteger, mustBeNonnegative}
    nSyncSamples  (1,1) double {mustBeInteger,   mustBeNonnegative}
    AcqPulseWidth (1,1) double {mustBeInteger,   mustBePositive}
    ExcPulseWidth (1,1) double {mustBeInteger,   mustBePositive}
    T_excitation  (1,1) double {mustBePositive}
    Ts_OCT        (1,1) double {mustBePositive}
    M_repetitions (1,1) double {mustBeInteger,   mustBePositive}
end

%% -------------------------- Derived quantities ---------------------------
samplesPerExcPeriod_raw = T_excitation / Ts_OCT;
samplesPerExcPeriod = round(samplesPerExcPeriod_raw);

if abs(samplesPerExcPeriod_raw - samplesPerExcPeriod) > 1e-10
    warning(['T_excitation / Ts_OCT is not an integer. It was rounded from %.12f to %d samples.'], ...
             samplesPerExcPeriod_raw, samplesPerExcPeriod);
end

fprintf('samplesPerExcPeriod = %d (T_exc = %.3f ms)\n', samplesPerExcPeriod, 1e3*T_excitation);

% Per-point dwell:
% acquisition dwell + optional compatibility hold
pointDwellSamples = M_repetitions + startDelay_samples;
% fprintf('pointDwellSamples = %d samples (%d acquisition + %d compatibility hold)\n', ...
%         pointDwellSamples, M_repetitions, startDelay_samples);

%% --------------------- Build multimeridian point pattern -----------------
% Same normalized geometry as LabVIEW code:
% r in [-1,1], theta in [0,pi), alternating radial order between meridians

r = linspace(-1,1,nAscans);
theta = linspace(0,pi,nBscans+1);
theta = theta(1:end-1);

[rr,tt] = meshgrid(r,theta);
rr = rr.';
tt = tt.';

% Alternate the radial direction for every second meridian
rr(:,2:2:end) = rr(end:-1:1,2:2:end);

[x_norm,y_norm] = pol2cart(tt,rr);

% Column-wise ordering reproduces the meridian-by-meridian acquisition order
x_norm_pts = x_norm(:);
y_norm_pts = y_norm(:);

% Map normalized coordinates [-1,1] -> physical mm
x_pts = ((x_norm_pts + 1)/2) * (horizontalMax - horizontalMin) + horizontalMin;
y_pts = ((y_norm_pts + 1)/2) * (verticalMax   - verticalMin)   + verticalMin;

nSpatialPoints = numel(x_pts);
fprintf('Total spatial points per frame = %d (%d meridians x %d points/meridian)\n', ...
        nSpatialPoints, nBscans, nAscans);

%% ----------------------- Preallocate full vectors ------------------------
% Samples per frame:
%   every point contributes pointDwellSamples
%   every jump between points contributes nSyncSamples
samplesPerFrame = nSpatialPoints*pointDwellSamples + (nSpatialPoints-1)*nSyncSamples;

% Between frames, add one extra smooth return to the first point
samplesBetweenFrames = nSyncSamples;

totalSamples = nFrames*samplesPerFrame + (nFrames-1)*samplesBetweenFrames;

GalvoX  = zeros(totalSamples,1);
GalvoY  = zeros(totalSamples,1);
Exc_aux = zeros(totalSamples,1);
Acq_aux = zeros(totalSamples,1);

%% ----------------------- Build the full waveform -------------------------
cursor = 0;

for iFrame = 1:nFrames

    % If not first frame, add smooth return from last point of previous frame
    % to the first point of the new frame
    if iFrame > 1 && nSyncSamples > 0
        p0 = [GalvoX(cursor), GalvoY(cursor)];
        p1 = [x_pts(1), y_pts(1)];
        [Xtrans,Ytrans] = localSmoothTransition(p0,p1,nSyncSamples);

        idx = cursor + (1:nSyncSamples);
        GalvoX(idx)  = Xtrans;
        GalvoY(idx)  = Ytrans;
        % triggers remain zero during transitions
        cursor = cursor + nSyncSamples;
    end

    for iPt = 1:nSpatialPoints

        x0 = x_pts(iPt);
        y0 = y_pts(iPt);

        % ------------------- Dwell at current spatial point ----------------
        idx_dwell = cursor + (1:pointDwellSamples);

        GalvoX(idx_dwell) = x0;
        GalvoY(idx_dwell) = y0;

        dwellStart = idx_dwell(1);
        dwellEnd   = idx_dwell(end);

        % Excitation TTL starts at the same time as acquisition
        exc_idx = dwellStart : min(dwellStart + ExcPulseWidth - 1, dwellEnd);
        Exc_aux(exc_idx) = 1;

        % Acquisition TTL repeated every A-line for M_repetitions samples
        acqStart = dwellStart;
        acqEnd   = dwellStart + M_repetitions - 1;

        % Pulse every OCT sample (one A-line per sample)
        for k = 0:(M_repetitions-1)
            pulseStart = acqStart + k;
            pulseEnd   = min(pulseStart + AcqPulseWidth - 1, acqEnd);
            Acq_aux(pulseStart:pulseEnd) = 1;
        end

        cursor = cursor + pointDwellSamples;

        % ---------------- Transition to next spatial point -----------------
        isLastPointOfFrame = (iPt == nSpatialPoints);

        if ~isLastPointOfFrame && nSyncSamples > 0
            p0 = [x_pts(iPt),   y_pts(iPt)];
            p1 = [x_pts(iPt+1), y_pts(iPt+1)];

            [Xtrans,Ytrans] = localSmoothTransition(p0,p1,nSyncSamples);

            idx_trans = cursor + (1:nSyncSamples);
            GalvoX(idx_trans) = Xtrans;
            GalvoY(idx_trans) = Ytrans;
            % triggers remain zero during transitions

            cursor = cursor + nSyncSamples;
        end
    end
end

% Safety trim in case of any mismatch
GalvoX  = GalvoX(1:cursor);
GalvoY  = GalvoY(1:cursor);
Exc_aux = Exc_aux(1:cursor);
Acq_aux = Acq_aux(1:cursor);

Exc_Trigger = double(Exc_aux(:));
Acq_Trigger = double(Acq_aux(:));
GalvoX      = double(GalvoX(:));
GalvoY      = double(GalvoY(:));

N_total = numel(GalvoX);
fprintf('Built Multimeridian MB waveform timeline: total samples = %d\n', N_total);

%% ----------------------------- Diagnostics -------------------------------
% Spatial point indices (frame 1 only) for visualization
frame1_point_start_idx = zeros(nSpatialPoints,1);
frame1_point_end_idx   = zeros(nSpatialPoints,1);
tmpCursor = 0;
for iPt = 1:nSpatialPoints
    frame1_point_start_idx(iPt) = tmpCursor + 1;
    frame1_point_end_idx(iPt)   = tmpCursor + pointDwellSamples;
    tmpCursor = tmpCursor + pointDwellSamples;
    if iPt < nSpatialPoints
        tmpCursor = tmpCursor + nSyncSamples;
    end
end

% 1) Spatial acquisition pattern
figure('Name','Multimeridian MB Pattern','Units','normalized','Position',[0.05 0.1 0.6 0.7]);
hold on;
plot(GalvoX, GalvoY, '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 1.0);

% plot nominal spatial points of first frame
scatter(x_pts, y_pts, 28, 1:nSpatialPoints, 'filled');
plot(x_pts, y_pts, 'k:', 'LineWidth', 0.8);

% highlight start and end
plot(x_pts(1), y_pts(1), 'go', 'MarkerFaceColor','g', 'MarkerSize',8);
plot(x_pts(end), y_pts(end), 'ro', 'MarkerFaceColor','r', 'MarkerSize',8);

axis equal;
grid on;
xlabel('X (mm)');
ylabel('Y (mm)');
title(sprintf('Multimeridian MB Scan Pattern (%d meridians, %d points/meridian)', nBscans, nAscans));
cb = colorbar;
colormap("jet")
ylabel(cb,'Spatial point order');
legend('Full galvo path','Nominal points','Point connection','Start','End','Location','best');
hold off;

% 2) Global timeline
t_ms = ((0:N_total-1)' * Ts_OCT) * 1e3;

figure('Name','Multimeridian MB Timeline','Units','normalized','Position',[0.05 0.05 0.85 0.35]);
yyaxis left;
plot(t_ms, GalvoX, '-b', 'LineWidth', 1.1); hold on;
plot(t_ms, GalvoY, '-r', 'LineWidth', 1.1);
ylabel('Galvo position (mm)');

yyaxis right;
smax = max(abs([GalvoX; GalvoY]));
if smax == 0
    smax = 1;
end
plot(t_ms, Exc_Trigger * 0.9*smax, '-k', 'LineWidth', 1.5);
plot(t_ms, Acq_Trigger * 0.6*smax, '-g', 'LineWidth', 1.2);
ylabel('Trigger scale (a.u.)');

xlabel('Time (ms)');
title('Multimeridian MB Acquisition Vectors - Global Timeline');
legend('Galvo X','Galvo Y','Excitation Trigger','Acquisition Trigger','Location','best');
grid on;
box on;
hold off;

% 3) Zoom into first few spatial points
nPtsToShow = min(4,nSpatialPoints);
zoomEnd = frame1_point_end_idx(nPtsToShow);
if nPtsToShow < nSpatialPoints
    zoomEnd = min(N_total, zoomEnd + nSyncSamples);
end

figure('Name','Multimeridian MB Timeline - First Spatial Points','Units','normalized','Position',[0.08 0.08 0.85 0.45]);
yyaxis left;
plot(t_ms(1:zoomEnd), GalvoX(1:zoomEnd), '-b', 'LineWidth', 1.2); hold on;
plot(t_ms(1:zoomEnd), GalvoY(1:zoomEnd), '-r', 'LineWidth', 1.2);
ylabel('Galvo position (mm)');

yyaxis right;
plot(t_ms(1:zoomEnd), Exc_Trigger(1:zoomEnd) * 0.9*smax, '-k', 'LineWidth', 1.6);
stem(t_ms(1:zoomEnd), Acq_Trigger(1:zoomEnd) * 0.6*smax, '-g', 'LineWidth', 1.2);
ylabel('Trigger scale (a.u.)');

for iPt = 1:nPtsToShow
    xline(t_ms(frame1_point_start_idx(iPt)), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.9);
end

xlabel('Time (ms)');
title(sprintf('First %d Spatial Points - MB Timing Detail', nPtsToShow));
legend('Galvo X','Galvo Y','Excitation Trigger','Acquisition Trigger','Location','best');
grid on;
box on;
hold off;

end

%% ============================ Local function =============================
function [Xtrans,Ytrans] = localSmoothTransition(p0,p1,L)
% Smooth cubic transition with zero slope at both ends.
% Returns L samples excluding the exact endpoints.

if L <= 0
    Xtrans = zeros(0,1);
    Ytrans = zeros(0,1);
    return;
end

s = (1:L)' ./ (L + 1);   % excludes endpoints
w = 3*s.^2 - 2*s.^3;     % smoothstep, C1-continuous with zero end slopes

Xtrans = p0(1) + (p1(1) - p0(1)) .* w;
Ytrans = p0(2) + (p1(2) - p0(2)) .* w;
end