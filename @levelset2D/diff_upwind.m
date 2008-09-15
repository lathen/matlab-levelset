function [Dxm,Dxp,Dym,Dyp] = diff_upwind(ls)

[Dxm,Dxp,Dym,Dyp] = ls.diff_upwind(ls.phi, ls.band);
