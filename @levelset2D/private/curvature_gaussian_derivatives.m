function k = curvature_gaussian_derivatives(ls)

G0 = gaussiankernel(7,0.7,0);
G1 = gaussiankernel(7,0.7,1);
G2 = gaussiankernel(7,0.7,2);

Dx = imfilter(ls.phi,G0','replicate','same','conv');
Dx = imfilter(Dx,G1,'replicate','same','conv');

Dy = imfilter(ls.phi,G0,'replicate','same','conv');
Dy = imfilter(Dy,G1','replicate','same','conv');

D2x = imfilter(ls.phi,G0','replicate','same','conv');
D2x = imfilter(D2x,G2,'replicate','same','conv');

D2y = imfilter(ls.phi,G0,'replicate','same','conv');
D2y = imfilter(D2y,G2','replicate','same','conv');

DxDy = imfilter(ls.phi,G1,'replicate','same','conv');
DxDy = imfilter(DxDy,G1','replicate','same','conv');

Dx = Dx(ls.band);
Dy = Dy(ls.band);
D2x = D2x(ls.band);
D2y = D2y(ls.band);
DxDy = DxDy(ls.band);

Dx2 = Dx.^2;
Dy2 = Dy.^2;

k = (Dx2.*D2y + Dy2.*D2x - 2*Dx.*Dy.*DxDy) ./ (sqrt(Dx2 + Dy2) .* (Dx2 + Dy2));
