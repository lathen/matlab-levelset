
function [ddt,dt] = min_curvature_flow_operator(ls, varargin)

if isempty(varargin{1})
    error('Need to specify magnitude of curvature flow');
elseif ~isnumeric(varargin{1})
    error('Magnitude of curvature flow is not numeric');
end

a = varargin{1};
if isscalar(a)
    a = a*ones(size(narrowband(ls)));
end

dt = 0.9/(6*max(abs(a(:))));

[Dx,Dy,Dz] = diff_central(ls);
k = min_curvature(ls);

ddt = a .* k .* sqrt(Dx.^2 + Dy.^2 + Dz.^2);
ddt(isnan(k)) = 0;
ddt(isinf(k)) = 0;
