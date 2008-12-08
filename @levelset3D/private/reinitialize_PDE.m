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
    B = permute(A, [3 2 1]);
    B = imdilate(B,SE);
    A = imdilate(A,SE);
    A = A + permute(B, [3 2 1]);

    ls.band = uint32(find(A > 0))';
end

% Run the reinitialize operator
[ls,iter] = propagate(ls,delta,'reinitialize_PDE_operator');

% Update the band
ls.band = ls.band((abs(ls.phi(ls.band)) <= ls.bandwidth));
