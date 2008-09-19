
function [ls,iterations, elapsed] = propagate(ls, time, operator, varargin)

persistent old_delta_phi;
if(isempty(old_delta_phi))
    old_delta_phi = zeros(size(ls));
    'hepp'
end

alpha = 0;
lr = 10; % Learning rate
elapsed = 0;
iterations = 0;


% Do one iteration if the requested time is 0
if (time == 0)

    old_phi = ls.phi;
    
    % Propagate the level set function in time
    [phi,dt] = ls.integrate(ls, Inf, operator, varargin{:});
    new_phi = ls.phi;
    new_phi(ls.band) = phi;
    
    % Update level set function and exit
    est_delta_phi = new_phi - ls.phi;
    delta_phi     = alpha * old_delta_phi + (1-alpha)*lr*est_delta_phi;

    ls.phi(ls.band) = ls.phi(ls.band) + delta_phi(ls.band);
    
    figure(43); hold off; clf;
    subplot(3,2,1);imagesc(old_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,2);imagesc(ls.phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,3);imagesc(old_delta_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,4);imagesc(delta_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,5);imagesc(abs(delta_phi- old_delta_phi));colorbar;hold on; plot(ls, 'contour y');

    
    %old_delta_phi = ls.phi - old_phi;
    old_delta_phi = delta_phi;
    
    
    ls = rebuild_narrowband(ls);
    %old_delta_phi = ls.phi - old_phi;
    %old_delta_phi = est_delta_phi;

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
