function ls = reinitialize_fastmarching_driver(ls)

[phi, ls.band] = reinitialize_fastmarching(ls.phi, ls.band, ls.bandwidth);
ls.phi(ls.band) = phi;
