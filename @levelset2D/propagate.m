
function [ls,iterations, elapsed] = propagate(ls, time, operator, varargin)

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


%% L
% alpha = 0.9
% lr = 6
% top = 7

%%retina
% alpha = 0.5
% lr = 1
% top = 2

alpha = 0.9;
%virt_m = 0.1;
lr = 6; % Learning rate
elapsed = 0;
iterations = 0;
top = 7;

% tmp = ls.phi;
% tmp(tmp >= 1)  = 1;
% tmp(tmp <= -1) = -1;
% ls.phi = tmp;

old_phi  = ls.phi;
old_band = ls.band;

% Do one iteration if the requested time is 0
if (time == 0)
    
    % Propagate the level set function in time
    [phi,dt] = ls.integrate(ls, Inf, operator, varargin{:});
    ls.phi(ls.band) = phi;
    ls = rebuild_narrowband(ls);
    
    % Update level set function and exit
    curr_delta_phi = ls.phi - old_phi;
    delta_phi      = alpha * old_delta_phi + (1-alpha)*lr*curr_delta_phi;
    delta_phi(delta_phi > top) = top;
    delta_phi(delta_phi < (-top)) = -top;
    
    %delta_phi     = old_delta_phi + curr_delta_phi / virt_m;
    old_delta_phi = delta_phi;

    ls.phi = old_phi + delta_phi;
    ls = rebuild_narrowband(ls);
    
    figure(44); hold off; clf;
    subplot(3,2,1);imagesc(old_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,2);imagesc(ls.phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,3);imagesc(old_delta_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,4);imagesc(lr*curr_delta_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,5);imagesc(abs(lr*curr_delta_phi- old_delta_phi));colorbar;hold on; plot(ls, 'contour y');
    drawnow;
    %pause;
    
    %old_delta_phi = ls.phi - old_phi;
    %old_delta_phi = delta_phi;
    
    
    %ls = rebuild_narrowband(ls);
    %old_delta_phi = ls.phi - old_phi;
    %old_delta_phi = est_delta_phi;

    elapsed = dt;
    iterations = 1;
    %pause;

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

    iterations
    ls = rebuild_narrowband(ls);
    
    % Update level set function and exit
    common_band = intersect(ls.band, old_band);
    [Y, X] = ind2sub(size(ls), common_band);
    %error('hepp')
    curr_delta_phi = griddata(double(X),double(Y),ls.phi(common_band) - old_phi(common_band),XI,YI,'nearest');
    
    
    delta_phi      = alpha * old_delta_phi + (1-alpha)*lr*curr_delta_phi;
    delta_phi(delta_phi > top) = top;
    delta_phi(delta_phi < (-top)) = -top;
    
    %delta_phi     = old_delta_phi + curr_delta_phi / virt_m;
    old_delta_phi = delta_phi;

    ls.phi = old_phi + delta_phi;
    ls = rebuild_narrowband(ls);
    
    figure(44); hold off; clf;
    subplot(3,2,1);imagesc(old_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,2);imagesc(ls.phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,3);imagesc(old_delta_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,4);imagesc(lr*curr_delta_phi);colorbar;hold on; plot(ls, 'contour y');
    subplot(3,2,5);imagesc(abs(lr*curr_delta_phi- old_delta_phi));colorbar;hold on; plot(ls, 'contour y');
    drawnow;
    
end

