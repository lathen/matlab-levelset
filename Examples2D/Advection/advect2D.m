%% Advection
% Example of advection of level set in an external vector field.

%% Create level set object
% First we create a level set object by four circles which initialize the
% zero level set.

% Create four circles to initialize the zero level set
[X,Y] = meshgrid(-50:50, -40:40);
phi = sqrt(X.^2 + Y.^2);
phi(phi > 25) = 25;
phi(phi < 24) = 24;
phi = phi - 24.5;
phi = repmat(phi, 2,2);

% Create level set object with narrowband of 16. We reinitialize to
% compute the signed distance function within the narrowband. The plot
% shows the zero level set as well as the level set function (green colors
% are positive, while red are negative)
LS = levelset2D(phi, 16);
LS = reinitialize(LS);
figure, plot(LS, 'contour', 'phi');


%% Create external vector field
% This vector field will advect the level set object.
S = size(LS);
rows = S(1); cols = S(2);
[Fx,Fy] = meshgrid([1:cols/2 cols/2:-1:1], [1:rows/2 rows/2:-1:1]);

figure;
d = 10;
quiver([1:d:cols], [1:d:rows], Fx(1:d:end,1:d:end), Fy(1:d:end,1:d:end), 'g');
hold on;
plot(LS);


%% Propagate the level set function
% Propagate the level set function a given amount of time using the level
% set PDE specified in the advect_operator function:
% <advect2D_operator.html advect2D_operator.m>
time = 0.1;
for i = 1:10
    LS = propagate(LS,time,'advect2D_operator', Fx,Fy);
    LS = reinitialize(LS);
    clf;
    quiver([1:d:cols], [1:d:rows], Fx(1:d:end,1:d:end), Fy(1:d:end,1:d:end), 'g');
    hold on;
    plot(LS);
    drawnow;
end
