function [ls,iterations,elapsed] = propagate_sgd(ls, time, td, nl, eta, operator, varargin)
%PROPAGATE_SGD  Propagate a levelset3D object a given amount of time 
%   using the stochastic gradient descent approach.
%   [LS,iterations,elapsed] = PROPAGATE_SGD(LS, time, td,nl,eta,operator, ...)

% td: time decay
% nl: noise level
% eta: learning rate (step 

'propagate_sgd_3d'

[grad, iterations, elapsed] = ls_calcgrad(ls, time, operator, varargin);

ls.phi(ls.band) = ls.phi(ls.band) + elapsed*eta*grad;
ls = reinitialize(ls);
