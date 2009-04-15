function [ls,iterations,elapsed] = propagate(ls, time, operator, varargin)
%PROPAGATE  Propagate a levelset2D object a given amount of time.
%   [LS,iterations,elapsed] = PROPAGATE(LS, time, operator, ...)
%   propagates LS a given amount of time using a level set PDE specified by
%   operator. If time is zero, the time integrator propagates the solution
%   by one stable time step (assuming explicit integration). The level set
%   PDE is preferably specified in a separate function, and operator is set
%   to the name of this function. An arbitrary number of arguments can be
%   given following the operator parameter. These arguments are passed to
%   the operator function.
%
%   The operator function must adhere to the prototype:
%     function [dphi_dt,dt] = levelsetPDE(LS, varargin)
%   where dphi_dt is the evaluated level set PDE, dt is a stable time step
%   for explicit time integration, LS is a levelset2D object and varargin
%   specifies the parameters of the level set PDE.

%   Author: Gunnar Läthén (gunnar.lathen@itn.liu.se)
%   $Date: 2007/10/17


% Start some counters
elapsed = 0;
iterations = 0;

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
