function [Dxm,Dxp,Dym,Dyp,Dzm,Dzp] = diff_upwind(ls)

[Dxm,Dxp,Dym,Dyp,Dzm,Dzp] = ls.diff_upwind(ls.phi, ls.band);
