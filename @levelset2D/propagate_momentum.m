function [ls,iterations,elapsed] = propagate_momentum(ls, time, omega, eta, top, operator, varargin)
%PROPAGATE_MOMENTUM  Propagate a levelset2D object a given amount of time
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
%   for explicit time integration, LS is a levelset2D object and varargin
%   specifies the parameters of the level set PDE.
%
%   The parameter omega specifies the momentum parameter, while eta
%   specifies the step length as described in the paper at
%   http://dmforge.itn.liu.se/ssvm09/
%   This strategy is mainly used for segmentation applications.

%   Author: Gunnar Läthén (gunnar.lathen@itn.liu.se)
%   $Date: 2008/10/01


% Set some persistent variables to store momentum for next call
persistent step_previous;
persistent domain_previous;
if isempty(step_previous)
    step_previous = zeros(size(ls));
    domain_previous = 1:numel(ls);
end

% Start counters
elapsed = 0;
iterations = 0;

% Save previous level set function and narrowband
phi_previous  = ls.phi;
band_previous = ls.band;

% Do one iteration if the requested time is 0
if (time == 0)
    
    % Propagate the level set function in time
    [phi,dt] = ls.integrate(ls, Inf, operator, varargin{:});
    
    % Update level set function and exit
    ls.phi(ls.band) = phi;
    elapsed = dt;
    iterations = 1;

% Else, iterate until requested time is reached
else
    while (elapsed < time)

        % Propagate the level set function in time
        [phi,dt] = ls.integrate(ls, time-elapsed, operator, varargin{:});

        % Update level set function and continue
        ls.phi(ls.band) = phi;
        elapsed = elapsed + dt;
        iterations = iterations + 1;
    end    
end

% Rebuild the distance function and the narrowband
ls = reinitialize(ls);

% Compute the current (approximate) gradient in the common domain,
% given current and previous time instances
domain = intersect(ls.band, band_previous);
grad = (ls.phi(domain) - phi_previous(domain)) /  elapsed;

% The domain of previous and current steps are different (since the
% narrowband has moved). To fix this, first compute the distance transform
% of the previous domain (this is only required in domain_diff, but we use
% the convenience of built-in bwdist) 
domain_diff = setdiff(domain, domain_previous);
BW = false(size(ls));
BW(domain_previous) = true;
[D,L] = bwdist(BW);

% Then, extend the previous step by picking the closest values
step_previous(domain_diff) = step_previous(L(domain_diff));

% Compute the step, incorporating momentum and the previous step
step = eta*(1-omega)*grad + omega*step_previous(domain);

% Cut the rate of change so we don't move too fast
step = 2*top ./ (1 + exp(-2*step/top)) - top;

% Save the current step and domain for next iteration
step_previous(domain) = step;
domain_previous = domain;

% Update level set function and reinitialize
ls.phi(domain) = phi_previous(domain) + elapsed*step;
ls = reinitialize(ls);

% Some plots for debugging
%  figure(44); hold off; clf;
%  subplot(3,2,1);imagesc(phi_previous);colorbar;hold on; plot(ls, 'contour y');
%  subplot(3,2,2);imagesc(ls.phi);colorbar;hold on; plot(ls, 'contour y');
%  subplot(3,2,3);imagesc(step_previous);colorbar;hold on; plot(ls, 'contour y');
%  tmp = zeros(size(ls)); tmp(domain) = grad;
%  subplot(3,2,4);imagesc(eta*tmp);colorbar;hold on; plot(ls, 'contour y');
%  tmp = zeros(size(ls)); tmp(domain_diff) = 1;
%  subplot(3,2,5);imshow(tmp, []);hold on; plot(ls, 'contour y');
%  drawnow;
