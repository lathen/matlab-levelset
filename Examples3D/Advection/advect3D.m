%% Advection
% Example of advection of level set in an external vector field.

%% Create level set object
% First we create a level set object by four spheres which initialize the
% zero level set.

% Create four spheres to initialize the zero level set
[X,Y,Z] = meshgrid(-10:10, -8:8, -9:9);
phi = sqrt(X.^2 + Y.^2 + Z.^2);
phi(phi > 5) = 5;
phi(phi < 4) = 4;
phi = phi - 4.5;
phi = repmat(phi, [2 2 2]);

% Create level set object with narrowband of 4. We reinitialize to
% compute the signed distance function within the narrowband. The plot
% shows the zero level set as well as the gradient in the narrowband.
narrowband = 4;
phi = phi*narrowband*2; % match the range of phi with narrowband
LS = levelset3D(phi, narrowband);
LS = reinitialize(LS);
figure, plot(LS, 'contour', 'gradient');
view(20,20); zoom(0.8);


%% Create external vector field
% This vector field will advect the level set object.
S = size(LS);
rows = S(1); cols = S(2); slices = S(3);
[Fx,Fy,Fz] = meshgrid([1:cols/2 cols/2:-1:1], ...
                      [1:rows/2 rows/2:-1:1], ...
                      [1:slices/2 slices/2:-1:1]);

figure;
d = 5;
[X,Y,Z] = meshgrid(1:d:cols, 1:d:rows, 1:d:slices);
quiver3(X,Y,Z, Fx(1:d:end,1:d:end,1:d:end), ...
               Fy(1:d:end,1:d:end,1:d:end), ...
               Fz(1:d:end,1:d:end,1:d:end), 'g');
hold on;
plot(LS);
view(20,20);


%% Propagate the level set function
% Propagate the level set function a given amount of time using the level
% set PDE specified in the advect_operator function:
% <advect3D_operator.html advect3D_operator.m>
figure;
plot(LS);
time = 0.1;
for i = 1:7
    LS = propagate(LS,time,'advect3D_operator', Fx,Fy,Fz);
    LS = reinitialize(LS);
    figure;
    plot(LS);
    drawnow;
end
