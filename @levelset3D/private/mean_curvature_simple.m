function k = mean_curvature_simple(ls)

[Dx,Dy,Dz] = ls.diff_central(ls.phi, ls.band);
[D2x,D2y,D2z,DxDy,DxDz,DyDz] = ls.diff2(ls.phi, ls.band);

Dx2 = Dx.^2;
Dy2 = Dy.^2;
Dz2 = Dz.^2;

denominator = 2*(Dx2 + Dy2 + Dz2).^1.5;
nominator = Dx2.*(D2y + D2z) - 2*Dy.*Dz.*DyDz + ...
            Dy2.*(D2x + D2z) - 2*Dx.*Dz.*DxDz + ...
            Dz2.*(D2x + D2y) - 2*Dx.*Dy.*DxDy;
k = nominator ./ denominator;
