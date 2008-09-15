
function [ls,iterations] = propagate(ls, time, operator, varargin)

elapsed = 0;
iterations = 0;
while (elapsed < time)
    
    % Propagate the level set function in time
    [phi,dt] = ls.integrate(ls, time, elapsed, operator, varargin{:});

    % Check for convergence if requested time is infinite
    if isinf(time)
        difference = abs((phi <= 0) - (ls <= 0));
        converge = sum(difference(:))
        if (converge == 0)
            break;
        end
    end

    % Update level set function and continue
    ls.phi(ls.band) = phi;
    elapsed = elapsed + dt;
    iterations = iterations + 1;

end
