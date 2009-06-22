function H = heaviside(ls)

eps = 1.5;

H = 0.5*(1 + ls.phi(ls.band)/eps + sin(pi*ls.phi(ls.band)/eps)/pi);
H(ls.phi(ls.band) < -eps) = 0;
H(ls.phi(ls.band) > eps) = 1;
