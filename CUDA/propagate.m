function [phi elapsed iterations] = propagate(phi, operator, integrator, time, dt, varargin)

% Start some counters
elapsed = 0;
iterations = 0;

dphidt = parallel.gpu.GPUArray.zeros(size(phi), 'single');
dim = gpuArray(int32(size(phi)));

% Do one iteration if the requested time is 0
if (time == 0)
    
    tic; dphidt = feval(operator, phi, dphidt, dim, varargin{:}); toc;
    tic; phi = feval(integrator, phi, dphidt, dt, dim); toc;

    elapsed = dt;
    iterations = 1;

% Else, iterate until requested time is reached
else
    while (elapsed < time)

        tic; dphidt = feval(operator, phi, dphidt, dim, varargin{:}); toc;
        tic; phi = feval(integrator, phi, dphidt, dt, dim); toc;

        % Update level set function and continue
        elapsed = elapsed + dt;
        iterations = iterations + 1;
    end
end
