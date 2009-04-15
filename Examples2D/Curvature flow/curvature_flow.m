%% Curvature flow
% Example of level set undergoing curvature flow.

%% Create level set object
% First we create a level set object by four circles which initialize the
% zero level set.

% Create four circles to initialize the zero level set
[X,Y] = meshgrid(-80:20, -60:20);
phi = sqrt(X.^2 + Y.^2);
phi(phi > 25) = 25;
phi(phi < 24) = 24;
phi = phi - 24.5;
phi = [phi fliplr(phi)];
phi = [phi; flipud(phi)];

% Create level set object with narrowband of 16. We reinitialize to
% compute the signed distance function within the narrowband. The plot
% shows the zero level set as well as the level set function (green colors
% are positive, while red are negative)
LS = levelset2D(phi, 16);
LS = reinitialize(LS);
figure, plot(LS, 'contour', 'phi');


%% Propagate the level set function
% Propagate the level set function a given amount of time using the level
% set PDE specified in the curvature_flow_operator function:
% <curvature_flow_operator.html curvature_flow_operator.m>
figure;
plot(LS);
time = 5;
for i = 1:5
    LS = propagate(LS,time,'curvature_flow_operator', 1);
    LS = reinitialize(LS);
    figure;
    plot(LS);
    drawnow;
end
