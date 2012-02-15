function [ls,iterations,elapsed] = propagate_momentum(ls, time, omega, eta, top, first_time, operator, varargin)
%PROPAGATE_MOMENTUM  Propagate a levelset3D object a given amount of time
%                    using a momentum approach.
%   [LS,iterations,elapsed] = PROPAGATE_MOMENTUM(LS, time, omage, eta, top,
%   operator, ...) propagates LS a given amount of time using a level set
%   PDE specified by operator. If time is zero, the time integrator
%   propagates the solution by one stable time step (assuming explicit
%   integration). The level set PDE is preferably specified in a separate
%   function, and operator is set to the name of this function. An
%   arbitrary number of arguments can be given following the operator
%   parameter. These arguments are passed to the operator function.
%
%   The operator function must adhere to the prototype:
%     function [dphi_dt,dt] = levelsetPDE(LS, varargin)
%   where dphi_dt is the evaluated level set PDE, dt is a stable time step
%   for explicit time integration, LS is a levelset3D object and varargin
%   specifies the parameters of the level set PDE.
%
%   The parameter omega specifies the momentum parameter, while eta
%   specifies the step length as described in the paper at
%   http://dmforge.itn.liu.se/ssvm09/
%   This strategy is mainly used for segmentation applications.


% Set some persistent variables to store momentum for next call
persistent step_previous;
persistent domain_previous;
if(isempty(step_previous) || first_time)
    step_previous = zeros(size(ls));
    domain_previous = 1:numel(ls);
end

% Save previous level set function and narrowband
phi_previous  = ls.phi;
domain        = ls.band;

[grad, iterations, elapsed] = ls_calcgrad(ls, time, operator, varargin{:});

% The domain of previous and current steps are different (since the
% narrowband has moved). To fix this, find the nearest neighbour
step_previous = ls_expandfield3d(step_previous, domain_previous, domain);

% Compute the step, incorporating momentum and the previous step
step = eta*(1-omega)*grad(domain) + omega*step_previous(domain);

% Cut the rate of change so we don't move too fast
step = min(step,top);

% Update level set function and reinitialize
ls.phi(domain) = phi_previous(domain) + step;
ls = reinitialize(ls);

% Save the current step and domain for next iteration
%step_previous(domain) = step;
step_previous = ls_calceffgrad(ls, phi_previous, domain, domain);
domain_previous = domain;
