
function [ddt,dt] = advect_operator(ls, varargin)

if size(varargin) < 3
    error('Need to specify vector field for advection');
elseif ~isnumeric(varargin{1}) || ~isnumeric(varargin{2}) || ~isnumeric(varargin{3})
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

% Expand to array if scalar
Fz = varargin{3};
if isscalar(Fz)
    Fz = Fz*ones(size(narrowband(ls)));
else
    Fz = Fz(narrowband(ls)); % extract components only within narrowband
end

% Determine stable timestep
normgrad = sqrt(Fx.^2 + Fy.^2 + Fz.^2);
dt = 0.9 / max(normgrad(:));

% Compute one-sided differences
[Dx_m,Dx_p,Dy_m,Dy_p,Dz_m,Dz_p] = diff_upwind(ls);

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

% Determine correct difference according to upwind direction
mask = Fz<0;
Dz = zeros(size(narrowband(ls)));
Dz(mask) = Dz_p(mask);
Dz(~mask) = Dz_m(~mask);

ddt = -(Dx.*Fx + Dy.*Fy + Dz.*Fz);
