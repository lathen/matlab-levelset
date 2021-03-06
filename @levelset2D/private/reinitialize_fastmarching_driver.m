function ls = reinitialize_fastmarching_driver(ls)

%band = ls.band;
[phi, ls.band] = reinitialize_fastmarching(ls.phi, ls.band, ls.bandwidth);

%A = zeros(size(band));
%A(ls.phi(band) > 0) = ls.bandwidth;
%A(ls.phi(band) < 0) = -ls.bandwidth;
%ls.phi(band) = A;

ls.phi(ls.phi > 0) = ls.bandwidth;
ls.phi(ls.phi < 0) = -ls.bandwidth;
ls.phi(ls.band) = phi;
