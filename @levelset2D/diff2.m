function [D2x,D2y,DxDy] = diff2(ls)

[D2x,D2y,DxDy] = ls.diff2(ls.phi, ls.band);
