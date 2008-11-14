
function [ls, iterations, elapsed] = propagate_rprop(ls, time, LR_MAX, LR_MIN, LR_0, top, first_time, operator, varargin)

acc_factor = 1.2; %Constant
dec_factor = 0.5; %Constant

alpha = 1;
beta  = 0.05;  %Gradient at fi==0 for the "anti-contour field"
w     = 5; %Inverse-Width of "anti-contour field" (large w --> field close to contour)
gamma = 4;
delta = 0;

% Set some persistent variables to store momentum for next call
persistent old_grad_phi; %Old gradient
persistent lr;           %Individual learning rates
persistent XI;
persistent YI;
persistent h_av;
persistent f;
persistent fpos;
persistent fneg;
if(first_time)
    rand('twister',sum(100*clock));
    old_grad_phi = zeros(size(ls));
    lr           = zeros(size(ls)) + LR_0;
    %lr           = rand(size(ls))*LR_0 + 2 * LR_0;
    nrows = size(old_grad_phi,1);
    ncols = size(old_grad_phi,2);
    [XI,YI] = meshgrid(1:ncols,1:nrows);
    XI = double(XI);
    YI = double(YI);
    h_av = fspecial('average',[3 3]);
    
    %% NB!! We assume static f here! Must change this if steering
    %% information changes!
    f = -varargin{1}; % Change sign because the filters is positive on the lines and negative outside...
    fpos = f >= 0; %The mask for positive (and zero) f
    fneg = f < 0;  %The mask for negative f
end

% Start counters
elapsed = 0;
iterations = 0;

% Save current level set phi and level set band
old_phi  = ls.phi;
old_band = ls.band;

mean_im = imfilter(ls.phi,h_av,'symmetric');

a_im = zeros(size(ls));
a_im(fpos) = 2*f(fpos).*(1-ls.phi(fpos));%Calculate gradient for positive f: 2*f*(1-fi)
a_im(fneg) = 2*f(fneg).*(1+ls.phi(fneg));%Calculate gradient for negative f: 2*f*(1+fi)
b_im       = 1./(1+(w*ls.phi).^2);
g_im       = -1*(ls.phi-mean_im);
d_im       = -1*ls.phi;

curr_grad_phi = alpha*a_im + beta*b_im + gamma*g_im + delta*d_im;

%Old anti-contour implementation, v2 :-)
%cont_im = colfilt(sign(ls.phi),[3 3],'sliding',@mykern);
%curr_grad_phi = -alpha * varargin{1} + beta*cont_im - gamma * (ls.phi - mean_im) - delta * ls.phi;

%Old anti-contour implementation, v1
%curr_grad_phi = -alpha * varargin{1} - gamma * (ls.phi - mean_im) - delta * ls.phi;
%con_i = (ls.phi >= -0.5) & (ls.phi <= 0.5);
%sum(sum(con_i))
%curr_grad_phi(con_i) = curr_grad_phi(con_i) + beta;

%Original
%curr_grad_phi = griddata(double(X),double(Y),ls.phi(common_band) - old_phi(common_band),XI,YI,'nearest');


%RPROP
grad_sprod = sign(old_grad_phi .* curr_grad_phi);
acc_i  = grad_sprod > 0;
%null_i = grad_sprod == 0;
dec_i  = grad_sprod < 0;

lr(acc_i) = min(lr(acc_i) * acc_factor, LR_MAX);
lr(dec_i) = max(lr(dec_i) * dec_factor, LR_MIN);
delta_phi = lr .* sign(curr_grad_phi);
delta_phi(dec_i) = 0; %In original RPROP, do not perform update if sign change
old_grad_phi = curr_grad_phi;
old_grad_phi(dec_i) = 0; %In original RPROP, do not adapt lr in next iteration if sign change.

% Cut the rate of change so we don't move too fast
delta_phi(delta_phi > top) = top;
delta_phi(delta_phi < (-top)) = -top;


% Update level set function and reinitialize
ls.phi = old_phi + delta_phi;
%ls = rebuild_narrowband(ls);

% Some plots for debugging
figure(44); hold off; clf;
subplot(3,3,1);imagesc(old_phi);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,2);imagesc(a_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,3);imagesc(b_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,4);imagesc(g_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,5);imagesc(d_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,6);imagesc(curr_grad_phi);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,7);imagesc(delta_phi);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,8);imagesc(grad_sprod);colorbar;hold on; %plot(ls, 'contour y');
%subplot(4,2,6);imagesc(cont_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,9);imagesc(ls.phi);colorbar;hold on; %plot(ls, 'contour y');

figure(45); hold off; clf;
subplot(3,3,1);imagesc(old_phi);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,2);imagesc(alpha*a_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,3);imagesc(beta*b_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,4);imagesc(gamma*g_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,5);imagesc(delta*d_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,6);imagesc(curr_grad_phi);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,7);imagesc(delta_phi);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,8);imagesc(grad_sprod);colorbar;hold on; %plot(ls, 'contour y');
%subplot(4,2,6);imagesc(cont_im);colorbar;hold on; %plot(ls, 'contour y');
subplot(3,3,9);imagesc(ls.phi);colorbar;hold on; %plot(ls, 'contour y');

drawnow;
if(first_time)
   pause;
end

function [y] = mykern(x)
nrows = size(x,1);
center = ceil(nrows/2.0);
a = repmat(x(center,:),[nrows 1]);
b = a.*x;
y = sum(b < 0);

