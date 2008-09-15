
function [ddt,dt] = advect(ls, varargin)

if size(varargin) < 2
    error('Need to specify vector field for advection');
elseif ~isnumeric(varargin{1}) || ~isnumeric(varargin{2})
    error('Vector field for advection is not numeric');
end

% Expand to array if scalar
Fx = varargin{1};
if isscalar(Fx)
    Fx = Fx*ones(size(narrowband(ls)));
else
    Fx = Fx(narrowband(ls)); % extract components only within narrowband
end

% Expand to array if scalar
Fy = varargin{2};
if isscalar(Fy)
    Fy = Fy*ones(size(narrowband(ls)));
else
    Fy = Fy(narrowband(ls)); % extract components only within narrowband
end

% Determine stable timestep
normgrad = sqrt(Fx.^2 + Fy.^2);
dt = 0.9 / max(normgrad(:));

% Compute differences for upwinding
[Dx_m,Dx_p,Dy_m,Dy_p] = diff_upwind(ls);

% Determine correct difference according to upwind direction
mask = Fx<0;
Dx = zeros(size(narrowband(ls)));
Dx(mask) = Dx_p(mask);
Dx(~mask) = Dx_m(~mask);

% Determine correct difference according to upwind direction
mask = Fy<0;
Dy = zeros(size(narrowband(ls)));
Dy(mask) = Dy_p(mask);
Dy(~mask) = Dy_m(~mask);

ddt = -(Dx.*Fx + Dy.*Fy);
