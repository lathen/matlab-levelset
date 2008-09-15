
function [ddt,dt] = minimal_variance(ls, varargin)

if isempty(varargin)
    error('Need to specify image to segment');
elseif ~isnumeric(varargin{1}) || ndims(varargin{1}) ~= 2
    error('Image is not a 2-D array (require a grayscale input)');
end

A = double(varargin{1});
A = A(narrowband(ls)); % extract image only within the narrowband
% TODO: we really want the inside of ls union with narrowband...

% Compute the heaviside function
H = 0.5*(1 + 2/pi*atan(-ls));

% Compute average inside (c1) and outside (c2)
nominator = A.*H;
denominator = H;
c1 = sum(nominator(:)) / sum(denominator(:));
nominator = A.*(1-H);
denominator = (1-H);
c2 = sum(nominator(:)) / sum(denominator(:));

% Compute variance inside and outside
A1 = A-c1;
A2 = A-c2;
F = (A1.^2 - A2.^2);

% Determine a stable timestep
dt = 0.9/max(F(:));

% Evaluate upwind gradient
[Dx2,Dy2] = godunov(ls,-F);
normgrad = sqrt(Dx2 + Dy2);

ddt = F.*normgrad;
