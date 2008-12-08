
function [ddt,dt] = speed_normal(ls, varargin)

if isempty(varargin)
    error('Need to specify magnitude of speed in normal direction');
elseif ~isnumeric(varargin{1})
    error('Magnitude of speed in normal direction is not numeric');
end

F = varargin{1};
if isscalar(F)
    F = F*ones(size(narrowband(ls)));
end

dt = 0.9/max(abs(F(:)));

[Dx2,Dy2] = godunov(ls,F);

ddt = -F .* sqrt(Dx2 + Dy2);
