function H = heaviside2(ls)

eps = 1.5;

H = 0.5*(1 + 2*atan(ls.phi(ls.band)/eps)/pi);
