function [phi,dt] = euler(ls, time, elapsed, operator, varargin)

% Evaluate operator and let the operator determine a stable time step
[ddt,dt] = feval(operator, ls, varargin{:});
if dt > time-elapsed
    dt = time-elapsed;
end

% Propagate the solution one time step
phi = ls + ddt.*dt;
