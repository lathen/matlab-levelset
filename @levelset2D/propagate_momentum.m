
function [ls, iterations, elapsed] = propagate_momentum(ls, time, alpha, lr, top, operator, varargin)

% Set some persistent variables to store momentum for next call
persistent old_delta_phi;
persistent XI;
persistent YI;
if(isempty(old_delta_phi))
    old_delta_phi = zeros(size(ls));
    nrows = size(old_delta_phi,1);
    ncols = size(old_delta_phi,2);
    [XI,YI] = meshgrid(1:ncols,1:nrows);
    XI = double(XI);
    YI = double(YI);
end

% Start counters
elapsed = 0;
iterations = 0;

% Save current level set phi and level set band
old_phi  = ls.phi;
old_band = ls.band;

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

iterations

% Rebuild the distance function and the narrowband
ls = rebuild_narrowband(ls);

% Compute the current gradient and extend values to the entire grid (if we
% have a narrowband)
common_band = intersect(ls.band, old_band);
[Y, X] = ind2sub(size(ls), common_band);
curr_delta_phi = griddata(double(X),double(Y),ls.phi(common_band) - old_phi(common_band),XI,YI,'nearest');

% Compute the rate of change given the previous gradient weighted with the
% current gradient (momentum)
delta_phi = alpha * old_delta_phi + (1-alpha)*lr*curr_delta_phi;

% Cut the rate of change so we don't move too fast
delta_phi(delta_phi > top) = top;
delta_phi(delta_phi < (-top)) = -top;

% Save the current gradient for next call
old_delta_phi = delta_phi;

% Update level set function and reinitialize
ls.phi = old_phi + delta_phi;
ls = rebuild_narrowband(ls);

% Some plots for debugging
figure(44); hold off; clf;
subplot(3,2,1);imagesc(old_phi);colorbar;hold on; plot(ls, 'contour y');
subplot(3,2,2);imagesc(ls.phi);colorbar;hold on; plot(ls, 'contour y');
subplot(3,2,3);imagesc(old_delta_phi);colorbar;hold on; plot(ls, 'contour y');
subplot(3,2,4);imagesc(lr*curr_delta_phi);colorbar;hold on; plot(ls, 'contour y');
subplot(3,2,5);imagesc(abs(lr*curr_delta_phi- old_delta_phi));colorbar;hold on; plot(ls, 'contour y');
drawnow;
