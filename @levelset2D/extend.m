function E = extend(ls,F)

% Closest-point-transform
E = F;
[X,Y] = CPT(ls,ls.band);

% Clamp to boundary
[rows,cols] = size(ls.phi);
X(X<1) = 1;
X(X>cols) = cols;
Y(Y<1) = 1;
Y(Y>rows) = rows;

E(ls.band) = interp2(F,X,Y);
