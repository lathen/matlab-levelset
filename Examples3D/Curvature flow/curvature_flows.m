%% Mean/min curvature flow
% Example of mean and minimum curvature flow of level set.

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
LS_mean = LS;
LS_min = LS;

figure;
plot(LS);
view(20,20);


%% Propagate the level set function
% Propagate the level set function a given amount of time using the level
% set PDEs specified in the mean_curvature_flow_operator and 
% min_curvature_flow_operator functions:
% <mean_curvature_flow_operator.html mean_curvature_flow_operator.m>
% <min_curvature_flow_operator.html min_curvature_flow_operator.m>
time = 1;
for i = 1:5
    LS_mean = propagate(LS_mean,time,'mean_curvature_flow_operator',1);
    LS_mean = reinitialize(LS_mean);
    LS_min = propagate(LS_min,time,'min_curvature_flow_operator',5);
    LS_min = reinitialize(LS_min);
    figure;
    subplot(1,3,1); plot(LS); view(20,20); title('Initial level set');
    subplot(1,3,2); plot(LS_mean); view(20,20); title('Mean curvature flow');
    subplot(1,3,3); plot(LS_min); view(20,20); title('Min curvature flow');
    drawnow;
end
