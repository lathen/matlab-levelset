function [F,V] = geometry(ls)

[F,V] = isosurface(ls.phi,0);
