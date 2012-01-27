function [ls, iterations, elapsed] = propagate_rprop(ls, time, LR_MAX, LR_MIN, LR_0, top, first_time, operator, varargin)

acc_factor = 1.2; %Constant
dec_factor = 0.5; %Constant

% Set some persistent variables to store momentum for next call
persistent old_grad_phi;  %Old gradient
persistent zero_field;
persistent old_grad_band; %Old gradient band
persistent lr;           %Individual learning rates

if(isempty(zero_field) || first_time)
   zero_field   = zeros(size(ls));
   old_grad_phi = zero_field;
   old_grad_band = ls.band; 
   lr           = zero_field + LR_0;
   lr(ls.band)  = LR_0;
end


% Save current level set phi and level set band
old_phi  = field(ls);
old_band = narrowband(ls);

[curr_grad_phi, iterations, elapsed] = ls_calcgrad(ls, time, operator, varargin{:});

delta_phi = lr(ls.band) .* sign(curr_grad_phi(ls.band));


% Cut the rate of change so we don't move too fast
delta_phi = min(delta_phi,top);
%delta_phi = 2*top ./ (1 + exp(-2*delta_phi/top)) - top;

% Update level set function and reinitialize
ls.phi(ls.band) = old_phi(ls.band) + delta_phi;

ls = reinitialize(ls);

%The sign we really took
ind = intersect(narrowband(ls),old_band);
eff_grad_phi = zero_field;
eff_grad_phi(ind) = ls.phi(ind) - old_phi(ind);


%Expand old_grad_phi to common domain




%RPROP
grad_sprod = sign(old_grad_phi .* eff_grad_phi);

%null_i = grad_sprod == 0;
%acc_i  = grad_sprod > 0;
%dec_i  = grad_sprod < 0;
sprod = sign(ls.phi .* eff_grad_phi);

acc_i  = (grad_sprod > 0) & (sprod < 0);
   %((ls.phi > 0 & eff_grad_phi < 0) | ...
   %(ls.phi < 0 & eff_grad_phi > 0));
dec_i  = (grad_sprod < 0) & (sprod > 0);
   %((ls.phi > 0 & eff_grad_phi > 0) | ...
   %(ls.phi < 0 & eff_grad_phi < 0));

dec_i_far = dec_i & (abs(ls.phi) > 2.0);

%figure(46);imagesc(grad_sprod);colorbar;
%figure(48);imagesc(eff_grad_phi);colorbar;
%figure(49);imagesc(curr_grad_phi);colorbar;
lr(acc_i)  = min(lr(acc_i)  * acc_factor, LR_MAX);
lr(dec_i)  = max(lr(dec_i)  * dec_factor, LR_MIN);
lr(dec_i_far)  = max(lr(dec_i_far), min(abs(ls.phi(dec_i_far))/4.0, LR_0));
%lr(dec_i_far)  = max(lr(dec_i_far), min(0.1, LR_0));

%lr((abs(ls.phi) <=1) & (lr == LR_MIN)) = 1;
%lr(dec_i)  = max(abs(eff_grad_phi(dec_i))  , LR_MIN);
%lr(null_i) = max(lr(null_i) * dec_factor, LR_MIN);
%lr(sign_disagrees) = 2;

%delta_phi(dec_i) = 0; %In original RPROP, do not perform update if sign change
%ls.phi(dec_i) = old_phi(dec_i);
%ls = reinitialize(ls);

old_grad_phi  = eff_grad_phi;
%old_grad_band = union_band;
old_grad_phi(dec_i) = 0; %In original RPROP, do not adapt lr in next iteration if sign change.



% Some plots for debugging
%figure(44); hold off; clf;
%subplot(4,2,1);imagesc(old_phi);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,2);imagesc(temp_phi);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,3);imagesc(curr_grad_phi);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,4);imagesc(delta_phi);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,5);imagesc(grad_sprod);colorbar;hold on; plot(ls, 'contour y');
%subplot(4,2,6);imagesc(ls.phi);colorbar;hold on; plot(ls, 'contour y');
%drawnow;
%pause;
