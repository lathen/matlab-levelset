function E = extend(ls,F)

% Closest-point-transform
E = F;
[X,Y,Z] = CPT(ls,ls.band);

% Clamp to boundary
[rows,cols,slices] = size(ls.phi);
X(X<1) = 1;
X(X>cols) = cols;
Y(Y<1) = 1;
Y(Y>rows) = rows;
Z(Z<1) = 1;
Z(Z>slices) = slices;

E(ls.band) = interp3(F,X,Y,Z);
