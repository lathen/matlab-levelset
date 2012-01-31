function k = curvature_simple(ls)

[Dx,Dy] = ls.diff_central(ls.phi, ls.band);
[D2x,D2y,DxDy] = ls.diff2(ls.phi, ls.band);

Dx2 = Dx.^2;
Dy2 = Dy.^2;

k = (Dx2.*D2y + Dy2.*D2x - 2*Dx.*Dy.*DxDy) ./ (sqrt(Dx2 + Dy2) .* (Dx2 + Dy2));
