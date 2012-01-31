function k = mean_curvature_gaussian_derivatives(ls)

G0 = gaussiankernel(7,0.7,0);
G1 = gaussiankernel(7,0.7,1);
G2 = gaussiankernel(7,0.7,2);

Dx = imfilter(ls.phi,G0','replicate','same','conv');
Dx = imfilter(Dx,reshape(G0,[1 1 7]),'replicate','same','conv');
Dx = imfilter(Dx,G1,'replicate','same','conv');

Dy = imfilter(ls.phi,G0,'replicate','same','conv');
Dy = imfilter(Dy,reshape(G0,[1 1 7]),'replicate','same','conv');
Dy = imfilter(Dy,G1','replicate','same','conv');

Dz = imfilter(ls.phi,G0,'replicate','same','conv');
Dz = imfilter(Dz,G0','replicate','same','conv');
Dz = imfilter(Dz,reshape(G1,[1 1 7]),'replicate','same','conv');

D2x = imfilter(ls.phi,G0','replicate','same','conv');
D2x = imfilter(D2x,reshape(G0,[1 1 7]),'replicate','same','conv');
D2x = imfilter(D2x,G2,'replicate','same','conv');

D2y = imfilter(ls.phi,G0,'replicate','same','conv');
D2y = imfilter(D2y,reshape(G0,[1 1 7]),'replicate','same','conv');
D2y = imfilter(D2y,G2','replicate','same','conv');

D2z = imfilter(ls.phi,G0,'replicate','same','conv');
D2z = imfilter(D2z,G0','replicate','same','conv');
D2z = imfilter(D2z,reshape(G2,[1 1 7]),'replicate','same','conv');

DxDy = imfilter(ls.phi,reshape(G0,[1 1 7]),'replicate','same','conv');
DxDy = imfilter(DxDy,G1,'replicate','same','conv');
DxDy = imfilter(DxDy,G1','replicate','same','conv');

DxDz = imfilter(ls.phi,G0','replicate','same','conv');
DxDz = imfilter(DxDz,G1,'replicate','same','conv');
DxDz = imfilter(DxDz,reshape(G1,[1 1 7]),'replicate','same','conv');

DyDz = imfilter(ls.phi,G0,'replicate','same','conv');
DyDz = imfilter(DyDz,G1','replicate','same','conv');
DyDz = imfilter(DyDz,reshape(G1,[1 1 7]),'replicate','same','conv');

Dx = Dx(ls.band);
Dy = Dy(ls.band);
Dz = Dz(ls.band);
D2x = D2x(ls.band);
D2y = D2y(ls.band);
D2z = D2z(ls.band);
DxDy = DxDy(ls.band);
DxDz = DxDz(ls.band);
DyDz = DyDz(ls.band);

Dx2 = Dx.^2;
Dy2 = Dy.^2;
Dz2 = Dz.^2;

denominator = 2*(Dx2 + Dy2 + Dz2).^1.5;
nominator = Dx2.*(D2y + D2z) - 2*Dy.*Dz.*DyDz + ...
            Dy2.*(D2x + D2z) - 2*Dx.*Dz.*DxDz + ...
            Dz2.*(D2x + D2y) - 2*Dx.*Dy.*DxDy;
k = nominator ./ denominator;
