function [grad,iterations,elapsed] = ls_calcgrad(ls, time, operator, varargin)

% Save previous level set function and narrowband
old_ls = ls;
phi_previous  = ls.phi;
band_previous = ls.band;

% Start some counters
elapsed = 0;
iterations = 0;

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

% Rebuild the distance function and the narrowband
ls = reinitialize(ls);

domain = intersect(ls.band, band_previous);
grad = (ls.phi(domain) - phi_previous(domain)) /  elapsed;

domain_diff = setdiff(domain, band_previous);

[rid,cid,sid] = ind2sub(size(ls.phi),domain_diff);
[rip,cip,sip] = ind2sub(size(ls.phi),band_previous);
[D,I] = pdist2(single([rip' cip' sip']),single([rid' cid' sid']),'euclidean','Smallest',1);

grad(domain_diff) = grad(band_previous(I));

figure(98);imagesc(old_ls.phi);colorbar;hold on; plot(old_ls, 'contour y');
figure(99);imagesc(ls.phi);colorbar;hold on; plot(ls, 'contour y');
figure(100);imagesc(grad);colorbar;
pause;
