function [CPx,CPy,CPz] = CPT(ls, points)

[Dx,Dy,Dz] = ls.diff_central(ls.phi, points);

[y,x,z] = ind2sub(size(ls), points);

CPx = double(x) - Dx.*ls.phi(points);
CPy = double(y) - Dy.*ls.phi(points);
CPz = double(z) - Dz.*ls.phi(points);
