function [expfield] = ls_expandfield3d(field, domain, expband)

domain_diff = setdiff(expband,domain);

[rid,cid,sid] = ind2sub(size(field),domain_diff);
[rip,cip,sip] = ind2sub(size(field),domain);
[D,I] = pdist2(single([rip' cip' sip']),single([rid' cid' sid']),'euclidean','Smallest',1);

expfield = field;
expfield(domain_diff) = expfield(domain(I));
