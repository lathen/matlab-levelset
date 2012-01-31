function k = mean_curvature_principal(ls)

[k1 k2] = principal_curvatures(ls.phi, ls.band);
k = (k1 + k2)/2;
