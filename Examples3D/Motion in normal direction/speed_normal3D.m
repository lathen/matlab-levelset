%% Level set operator defining a PDE for motion in normal direction 
% Implements the level set PDE for dilation used in the example:
% <dilation3D.html dilation3D.m>
function [dphi_dt,dt] = speed_normal3D(ls, varargin)

% Assume that the speed is passed as argument
F = varargin{1};
if isscalar(F) % expand if given as a scalar
    F = F*ones(size(narrowband(ls)));
else % restrict to narrowband if given as a field
    F = F(narrowband(ls));
end

% Determine stable timestep for explicit time integration
dt = 0.9/max(abs(F(:)));

% Compute upwind differences using Godunov's method
[Dx2,Dy2,Dz2] = godunov(ls,F);

% Compute and return level set PDE (dphi_dt)
dphi_dt = -F .* sqrt(Dx2 + Dy2 + Dz2);
