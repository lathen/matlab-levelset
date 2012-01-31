function [ls,iterations,elapsed] = propagate_sgd(ls, time, td, nl, eta, top, first_time,operator, varargin)
%PROPAGATE_SGD  Propagate a levelset3D object a given amount of time 
%   using the stochastic gradient descent approach.
%   [LS,iterations,elapsed] = PROPAGATE_SGD(LS, time, td,nl,eta,firs_time,operator, ...)

persistent timemap
if(isempty(timemap) || first_time)
   timemap = zeros(size(ls.phi));
end

[grad, iterations, elapsed] = ls_calcgrad(ls, time, operator, varargin{:});

tm = timemap(ls.band);
tm = max(tm - 50,0);

step = eta*grad(ls.band) + nl*randn(size(ls.band)).*exp(-td*tm);
step = min(step,top);

timemap(ls.band) = timemap(ls.band) + 1;
%figure(100);imagesc(tm);colorbar

ls.phi(ls.band) = ls.phi(ls.band) + step;
ls = reinitialize(ls);

