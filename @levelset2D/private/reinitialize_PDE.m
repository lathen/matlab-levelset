function ls = reinitialize_PDE(ls)

% Create binary mask for dilation of the narrow band
A = zeros(size(ls));
A(narrowband(ls)) = 1;

% Dilate the narrow band if the width of the band is finite
delta = ls.bandwidth;
if (delta == Inf)
    delta = max(size(ls));
else
    SE = strel('disk',delta);
    A = imdilate(A,SE);

    ls.band = uint32(find(A == 1))';
end

% Run the reinitialize operator
[ls,iter,elapsed] = propagate(ls,delta,'reinitialize_PDE_operator');

% Update the band
ls.band = ls.band((abs(ls.phi(ls.band)) <= ls.bandwidth));
