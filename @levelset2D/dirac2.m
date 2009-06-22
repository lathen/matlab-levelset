function d = dirac2(ls)

eps = 1.5;

d = eps./(pi*(eps^2 + ls.phi(ls.band).^2));
