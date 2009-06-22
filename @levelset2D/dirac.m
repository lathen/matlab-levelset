function d = dirac(ls)

eps = 1.5;

d = 0.5/eps*(1 + cos(pi*ls.phi(ls.band)/eps));
d(ls.phi(ls.band) < -eps) = 0;
d(ls.phi(ls.band) > eps) = 0;
