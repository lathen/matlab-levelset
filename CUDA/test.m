
%% Create helix to initialize the zero level set
[X,Y,Z] = meshgrid(-100:100, -100:100, -40:200);
phi = Inf*ones(size(X));
for t = 0 : 0.02*pi : 6*pi
    R = 60;
    r = 8;
    x = R*cos(t);
    y = R*sin(t);
    z = 8*t;
    sphere = sqrt((X-x).^2 + (Y-y).^2 + (Z-z).^2);
    sphere(sphere > r+1) = r+1;
    sphere(sphere < r) = r;
    sphere = sphere - r - 0.5;
    
    phi = min(phi,sphere);
end

phi = single(phi);
[rows, cols, slices] = size(phi);


%%
blockDim = [30 6 2];

[vol offset] = extendvolumeGPU(phi, blockDim, [1 1 1]);

k_reinitialize = make_kernel('ls_reinitialize', vol, blockDim, 1); 
k_speednormal = make_kernel('ls_speednormal', vol, blockDim, 1); 
k_integrate = make_kernel('ls_integrateeuler', vol, blockDim, 0); 



%%
vol = propagate(vol, k_reinitialize, k_integrate, 10, 0.5);

V = gather(vol);
figure, isosurface(V,0);
figure, imshow(V(:,:,120), []);

%%
F = parallel.gpu.GPUArray.zeros(size(vol), 'single');
a = 1;
dt = 0.9/(6*a);
vol = propagate(vol, k_speednormal, k_integrate, 20, dt, F, a);

V = gather(vol);
figure, isosurface(V,0);
figure, imshow(V(:,:,120), []);
