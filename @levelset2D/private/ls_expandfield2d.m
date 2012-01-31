function [expfield] = ls_expandfield2d(field, domain, expband)

domain_diff = setdiff(expband,domain);

[rid,cid] = ind2sub(size(field),domain_diff);
[rip,cip] = ind2sub(size(field),domain);
[D,I] = pdist2(single([rip' cip']),single([rid' cid']),'euclidean','Smallest',1);

expfield = field;
expfield(domain_diff) = expfield(domain(I));
