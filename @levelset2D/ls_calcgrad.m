function [grad,iterations,elapsed] = ls_calcgrad(ls, time, operator, varargin)

% Save previous level set function and narrowband
phi_previous  = ls.phi;
band_previous = ls.band;

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

% Rebuild the distance function and the narrowband
ls = reinitialize(ls);

grad = ls_calceffgrad(ls, phi_previous, band_previous,band_previous);
%grad = grad / elapsed;
