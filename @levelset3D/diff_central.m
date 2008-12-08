function [Dx,Dy,Dz] = diff_central(ls)

[Dx,Dy,Dz] = ls.diff_central(ls.phi, ls.band);
