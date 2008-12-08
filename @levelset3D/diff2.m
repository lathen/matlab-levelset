function [D2x,D2y,D2z,DxDy,DxDz,DyDz] = diff2(ls)

[D2x,D2y,D2z,DxDy,DxDz,DyDz] = ls.diff2(ls.phi, ls.band);
