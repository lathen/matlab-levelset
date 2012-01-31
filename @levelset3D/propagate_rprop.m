function [ls, iterations, elapsed] = propagate_rprop(ls, time, LR_MAX, LR_MIN, LR_0, top, first_time, operator, varargin)

acc_factor = 1.2; %Constant
dec_factor = 0.5; %Constant

% Set some persistent variables to store momentum for next call
persistent zero_field;
persistent old_grad_phi;  %Old gradient
persistent old_grad_band; %Old gradient band
persistent lr;           %Individual learning rates

if(isempty(zero_field) || first_time)
   zero_field    = zeros(size(ls));
   old_grad_phi  = zero_field;
   old_grad_band = ls.band; 
   lr            = zero_field + LR_0;
   lr(ls.band)   = LR_0;
end

% Save current level set phi and level set band
old_phi  = field(ls);
old_band = narrowband(ls);

[curr_grad_phi, iterations, elapsed] = ls_calcgrad(ls, time, operator, varargin{:});
delta_phi = lr(ls.band) .* sign(curr_grad_phi(ls.band));

% Cut the rate of change so we don't move too fast
delta_phi = min(delta_phi,top);

% Update level set function and reinitialize
ls.phi(ls.band) = old_phi(ls.band) + delta_phi;
ls = reinitialize(ls);

expband      = union(union(old_band, old_grad_band),ls.band);
eff_grad_phi = ls_calceffgrad(ls, old_phi, old_band, expband);
expoldgrad   = ls_expandfield3d(old_grad_phi, old_grad_band, expband);
expnewphi    = ls_expandfield3d(ls.phi, ls.band, expband);

%RPROP
%-----
grad_sprod = zero_field;
grad_sprod(expband) = sign(expoldgrad(expband) .* eff_grad_phi(expband));

sprod = zero_field;
sprod(expband) = sign(expnewphi(expband) .* eff_grad_phi(expband));

acc_i     = (grad_sprod > 0) & (sprod < 0);
dec_i     = (grad_sprod < 0) & (sprod > 0);
dec_i_far = dec_i & (abs(expnewphi) > 2.0);

lr(acc_i)     = min(lr(acc_i)  * acc_factor, LR_MAX);
lr(dec_i)     = max(lr(dec_i)  * dec_factor, LR_MIN);
lr(dec_i_far) = max(lr(dec_i_far), min(abs(expnewphi(dec_i_far))/4.0, LR_0));

old_grad_phi            = zero_field;
old_grad_phi(old_band)  = eff_grad_phi(old_band);
old_grad_band           = old_band;
old_grad_phi(dec_i)     = 0; %In original RPROP, do not adapt lr in next iteration if sign change.
