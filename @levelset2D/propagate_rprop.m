function [ls, iterations, elapsed] = propagate_rprop(ls, time, LR_MAX, LR_MIN, LR_0, top, first_time, operator, varargin)

acc_factor = 1.2; %Constant
dec_factor = 0.5; %Constant

% Set some persistent variables to store momentum for next call
persistent old_grad_phi;  %Old gradient
%persistent old_grad_band; %Old gradient band
persistent lr;           %Individual learning rates

if(first_time)
   %rand('twister',sum(100*clock));
   old_grad_phi  = zeros(size(ls));
   %old_grad_band = ls.band; 
   lr           = zeros(size(ls)) + LR_MAX;
   lr(ls.band)  = LR_0;
   %lr           = rand(size(ls))*LR_0 + 2 * LR_0;
end

% Start counters
elapsed = 0;
iterations = 0;

% Save current level set phi and level set band
old_phi  = ls.phi;
old_band = ls.band;
%old_D    = bwdist(old_phi <= 0) - bwdist(old_phi >= 0);
%figure;imagesc(old_D);pause
%figure(44);imagesc(lr);colorbar;
%figure(45);imagesc(old_phi);colorbar;
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
ls = reinitialize(ls);

% Compute the current gradient and extend values to the entire grid (if we
% have a narrowband)
%tmpi_band = intersect(ls.band, old_band);
%tmpu_band = union(ls.band, old_band);

%[Y, X] = ind2sub(size(ls), common_band);
%curr_grad_phi = griddata(double(X),double(Y),ls.phi(common_band) - old_phi(common_band),XI,YI,'nearest');
%curr_grad_phi = ls.phi(tmpu_band) - old_phi(tmpu_band);
curr_grad_phi = ls.phi - old_phi;
%curr_grad_phi = ls.phi - dummy;

%Only use the above to calculate the gradient
%temp_phi  = ls.phi;
%temp_band = ls.band;
ls.phi    = old_phi;
ls.band   = old_band;

delta_phi = lr .* sign(curr_grad_phi);

% Cut the rate of change so we don't move too fast
%delta_phi(delta_phi > top) = top;
%delta_phi(delta_phi < (-top)) = -top;
delta_phi = 2*top ./ (1 + exp(-2*delta_phi/top)) - top;

% Update level set function and reinitialize
%ls.phi(tmpu_band) = old_phi(tmpu_band) + delta_phi;
ls.phi = old_phi + delta_phi;

ls = reinitialize(ls);

%The sign we really took
%union_band = union(ls.band, old_band);
real_grad_phi = ls.phi - old_phi;
%figure(100); hold off; clf;
%imagesc(real_grad_phi);colorbar;hold on; plot(ls, 'contour y');

%sign_disagrees = sign(curr_grad_phi) ~= sign(real_grad_phi);
%figure(102); hold off; clf;
%imagesc(sign_disagrees);colorbar;hold on; plot(ls, 'contour y');

eff_grad_phi = real_grad_phi;

%union_time_band = union(union_band, old_grad_band);

%RPROP
grad_sprod = sign(old_grad_phi .* eff_grad_phi);

acc_i  = grad_sprod > 0;
%null_i = grad_sprod == 0;
dec_i  = grad_sprod < 0;

lr(acc_i)  = min(lr(acc_i)  * acc_factor, LR_MAX);
lr(dec_i)  = max(lr(dec_i)  * dec_factor, LR_MIN);
%lr(null_i) = max(lr(null_i) * dec_factor, LR_MIN);
%lr(sign_disagrees) = 2;

%delta_phi(dec_i) = 0; %In original RPROP, do not perform update if sign change
old_grad_phi  = eff_grad_phi;
%old_grad_band = union_band;
%old_grad_phi(dec_i) = 0; %In original RPROP, do not adapt lr in next iteration if sign change.



% Some plots for debugging
%figure(44); hold off; clf;
%subplot(4,2,1);imagesc(old_phi);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,2);imagesc(temp_phi);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,3);imagesc(curr_grad_phi);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,4);imagesc(delta_phi);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,5);imagesc(grad_sprod);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,6);imagesc(ls.phi);colorbar;hold on; plot(ls, 'contour y');
drawnow;
%pause;
