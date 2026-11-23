function U = TVDenoise_ChambollePock2D(F, mu, Mask, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TVDenoise_ChambollePock2D
%
% Solves approximately:
%
%       min_U 0.5 ||U - F||_2^2 + mu * TV(U)
%
% using a Chambolle-Pock primal-dual scheme.
%
% This file is intentionally kept separate so it can later be reused for:
%   - PG-TV
%   - PG-TNV extensions
%   - speed-map post-regularization tests
%   - wave-number-map denoising experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.FunctionName = 'TVDenoise_ChambollePock2D';

addParameter(p, 'NIterations', 250);
addParameter(p, 'Tau', 0.25);
addParameter(p, 'Sigma', 0.25);
addParameter(p, 'Theta', 1.0);

parse(p, varargin{:});
opts = p.Results;

if mu == 0
    U = F;
    return;
end

Mask = logical(Mask);

Ffilled = FillNaNs_Nearest2D(F);
Ffilled(~Mask) = 0;

U = Ffilled;
U(~Mask) = 0;

Ubar = U;

Px = zeros(size(F));
Py = zeros(size(F));

for it = 1:opts.NIterations
    [Gx, Gy] = ForwardGradient2D(Ubar);

    Px = Px + opts.Sigma * Gx;
    Py = Py + opts.Sigma * Gy;

    NormP = max(1, sqrt(Px.^2 + Py.^2) / mu);

    Px = Px ./ NormP;
    Py = Py ./ NormP;

    DivP = Divergence2D(Px, Py);

    Uold = U;

    U = (U + opts.Tau * DivP + opts.Tau * Ffilled) / (1 + opts.Tau);

    U(~Mask) = Ffilled(~Mask);

    Ubar = U + opts.Theta * (U - Uold);
end

U(~Mask) = NaN;

end

%% ======================================================================= %
%                          LOCAL FUNCTIONS
% ======================================================================== %

function [Gx, Gy] = ForwardGradient2D(U)

Gx = zeros(size(U));
Gy = zeros(size(U));

Gx(1:end-1,:) = U(2:end,:) - U(1:end-1,:);
Gy(:,1:end-1) = U(:,2:end) - U(:,1:end-1);

end

function Div = Divergence2D(Px, Py)

Div = zeros(size(Px));

Div(1,:) = Px(1,:);
Div(2:end-1,:) = Px(2:end-1,:) - Px(1:end-2,:);
Div(end,:) = -Px(end-1,:);

Div(:,1) = Div(:,1) + Py(:,1);
Div(:,2:end-1) = Div(:,2:end-1) + Py(:,2:end-1) - Py(:,1:end-2);
Div(:,end) = Div(:,end) - Py(:,end-1);

end

function Ufilled = FillNaNs_Nearest2D(U)

Ufilled = U;
bad = ~isfinite(U);

if ~any(bad(:))
    return;
end

good = isfinite(U);

if ~any(good(:))
    Ufilled(:) = 0;
    return;
end

[~, idx] = bwdist(good);
Ufilled(bad) = U(idx(bad));

end