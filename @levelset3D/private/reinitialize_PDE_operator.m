function [ddt,dt] = reinitialize_PDE_operator(ls, varargin)

% Compute the sign function
[Dx,Dy,Dz] = diff_central(ls);
normgrad2 = Dx.^2 + Dy.^2 + Dz.^2;
sign = ls ./ sqrt(ls.^2 + normgrad2);

[Dx2,Dy2,Dz2] = godunov(ls,sign);

dt = 0.5;
ddt = sign .* (1 - sqrt(Dx2 + Dy2 + Dz2));
