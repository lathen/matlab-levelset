

dim = 64;
A = zeros(dim,dim,dim);

radius = dim/3;
loops = 5;
tube_radius = 3;
z_b = tube_radius*2;

phi = 0:pi/180 / loops/2:2*pi;
n = length(phi);
x = round(dim/2 + radius*cos(loops*phi));
y = round(dim/2 + radius*sin(loops*phi));
z = round(z_b:(dim-1-z_b*2)/n:dim-z_b);

plot3(x,y,z(1:n));
axis([0 dim 0 dim 0 dim]);

ind = sub2ind(size(A), x,y,z(1:n));
A(ind) = 1;
A = A>0;

B = bwdist(A);
B = B-tube_radius;
B = smooth3(B, 'gaussian', 9);
figure; isosurface(B,0);

% Create levelset3D and plot
%LS = levelset3D(A,5);
%LS = levelset3D(A,5, 'Euler', 'FirstOrder', 'PDE');
LS = levelset3D(B,5, 'Euler', 'FirstOrder', 'FastMarching');
figure;
plot(LS, 'contour');

% Check largest gradient  magnitude
[Dx,Dy,Dz] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2 + Dz.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))

% Reinitialize
tic; LS = rebuild_narrowband(LS); toc
figure;
plot(LS, 'contour', 'narrowband 5');

% Check largest gradient  magnitude
[Dx,Dy,Dz] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2 + Dz.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))

LS_mean = levelset3D(LS);

% Propagate
time = 5;
for i = 1:30
    tic; [LS,iter] = propagate(LS,time,'min_curvature_flow',1); toc
    tic; LS = rebuild_narrowband(LS); toc

    %tic; [LS_mean,iter] = propagate(LS_mean,time,'mean_curvature_flow',1); toc
    %tic; LS_mean = rebuild_narrowband(LS_mean); toc

    figure(100);
    clf;
    subplot(121), plot(LS); view(28,30);
    %subplot(122), plot(LS_mean); view(28,30);
    drawnow;
end

