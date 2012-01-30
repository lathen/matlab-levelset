function H = heaviside_full(ls)

eps = 1.5;

H = 0.5*(1 + ls.phi/eps + sin(pi*ls.phi/eps)/pi);
H(ls.phi < -eps) = 0;
H(ls.phi > eps) = 1;
