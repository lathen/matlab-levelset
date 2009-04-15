function [CPx,CPy] = CPT(ls, points)

[Dx,Dy] = ls.diff_central(ls.phi, points);

[y,x] = ind2sub(size(ls), points);

CPx = double(x) - Dx.*ls.phi(points);
CPy = double(y) - Dy.*ls.phi(points);
