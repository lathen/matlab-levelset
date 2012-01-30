function k = gaussiankernel(n,sigma,order)

radius = (n-1)/2;
x = -radius:radius;

G = exp(-x.^2/(2*sigma^2));
G = G/sum(G);

switch order
    case 0
        k = G;
    case 1
        k = -x/sigma^2 .* G;
    case 2
        k = (x-sigma).*(x+sigma)/sigma^4 .* G;
    otherwise
        disp(['Gaussian derivative of order ',int2str(order),' not defined']);
        k = [];
end
