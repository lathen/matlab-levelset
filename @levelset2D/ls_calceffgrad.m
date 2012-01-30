function [effgrad] = ls_calceffgrad(ls, old_phi, old_band, expband)

persistent zerofield
if(isempty(zerofield))
   zerofield = zeros(size(ls.phi));
end

domain          = intersect(ls.band, old_band);
effgrad         = zerofield;
effgrad(domain) = (ls.phi(domain) - old_phi(domain));

effgrad=ls_expandfield2d(effgrad, domain, expband);
