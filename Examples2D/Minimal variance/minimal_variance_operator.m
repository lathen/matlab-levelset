%% Level set operator defining a PDE for minimal variance
% Implements the level set PDE for minimal variance used in the example:
% <minimal_variance.html minimal_variance.m>

function [dphi_dt,dt] = minimal_variance_operator(ls, varargin)

% Assume that the target image is passed as the first argument
A = double(varargin{1});

% Check if alpha parameter for regularization is given as a parameter
alpha = 0;
if length(varargin) > 1
    alpha = varargin{2};
end

% Compute the heaviside function
H = heaviside(ls);

% Compute average inside of the contour (c1)
nominator = A.*H;
denominator = H;
c1 = sum(nominator(:)) / sum(denominator(:));

% Compute average outside of the contour (c2)
nominator = A.*(1-H);
denominator = (1-H);
c2 = sum(nominator(:)) / sum(denominator(:));

% Restrict remaining computations to narrowband
A = A(narrowband(ls));

% Compute variance inside and outside
A1 = A-c1;
A2 = A-c2;
F = (A1.^2 - A2.^2);

% Determine stable timestep for explicit time integration
if alpha == 0
    dt = 0.9/max(abs(F(:)));
else
    dt = min(0.9/max(abs(F(:))), 0.9/(4*abs(alpha)));
end

% Evaluate upwind gradient
[Dx2,Dy2] = godunov(ls,-F);
normgrad = sqrt(Dx2 + Dy2);

% Compute and return level set PDE (dphi_dt)
if alpha == 0 % without regularization
    dphi_dt = F.*normgrad;
else % with regularization (add curvature flow term)
    [Dx,Dy] = diff_central(ls);
    k = curvature(ls);
    dphi_dt = F.*normgrad + alpha * k .* sqrt(Dx.^2 + Dy.^2);
    dphi_dt(isnan(k)) = 0;
    dphi_dt(isinf(k)) = 0;
end
