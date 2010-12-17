function L = length(ls)

ind = (abs(ls.phi) <= 1);
L = 0.5*(1.0 + cos(pi*ls.phi(ind)));
L = sum(L(:));
