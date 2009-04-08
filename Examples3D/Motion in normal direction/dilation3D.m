%% Dilation
% Example of dilation of level set.

%% Create level set object
% First we create a level set object by a helix which initializes the
% zero level set.

% Create helix to initialize the zero level set
[X,Y,Z] = meshgrid(-25:25, -26:26, -10:50);
phi = Inf*ones(size(X));
for t = 0 : 0.02*pi : 6*pi
    R = 15;
    r = 2;
    x = R*cos(t);
    y = R*sin(t);
    z = 2*t;
    sphere = sqrt((X-x).^2 + (Y-y).^2 + (Z-z).^2);
    sphere(sphere > r+1) = r+1;
    sphere(sphere < r) = r;
    sphere = sphere - r - 0.5;
    
    phi = min(phi,sphere);
end

% Create level set object with narrowband of 4. We reinitialize to
% compute the signed distance function within the narrowband.
narrowband = 4;
phi = phi*narrowband*2; % match the range of phi with narrowband
LS = levelset3D(phi, narrowband);
LS = reinitialize(LS);

figure;
plot(LS);
view(20,20);


%% Propagate the level set function
% Propagate the level set function a given amount of time using the level
% set PDE specified in the speed_normal3D function:
% <speed_normal3D.html speed_normal3D.m>
time = 1;
for i = 1:5
    LS = propagate(LS,time,'speed_normal3D',1);
    LS = reinitialize(LS);
    clf;
    plot(LS);
    view(20,20);
    drawnow;
end
