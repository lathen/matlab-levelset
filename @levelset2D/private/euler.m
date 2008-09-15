function [phi,dt] = euler(ls, max_dt, operator, varargin)

% Evaluate operator and let the operator determine a stable time step
[ddt,dt] = feval(operator, ls, varargin{:});

% If stable timestep exceeds max requested time, clamp...
if dt > max_dt
    dt = max_dt;
end

% Propagate the solution one time step
phi = ls + ddt.*dt;
