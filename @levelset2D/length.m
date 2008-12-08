function L = length(ls)

L = 0.5*(1.0 + cos(pi*ls.phi));
L(abs(ls.phi) > 1) = 0;

L = sum(L(:));
