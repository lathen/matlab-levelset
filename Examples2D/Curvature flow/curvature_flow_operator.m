%% Level set operator defining a PDE for curvature flow
% Implements the level set PDE for curvature flow used in the example:
% <curvature_flow.html curvature_flow.m>

function [dphi_dt,dt] = curvature_flow_operator(ls, varargin)

% Assume that the alpha parameter is passed as an argument
alpha = varargin{1};

% Determine stable timestep for explicit time integration
dt = 0.9/(4*abs(alpha));

% Compute central differences and curvature
[Dx,Dy] = diff_central(ls);
k = curvature(ls);

% Compute and return level set PDE (dphi_dt)
% Check for shocks/singularities (when curvature is Inf)
dphi_dt = alpha * k .* sqrt(Dx.^2 + Dy.^2);
dphi_dt(isnan(k)) = 0;
