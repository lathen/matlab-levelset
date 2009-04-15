%% Erosion and dilation
% Example of erosion and dilation of level set.

%% Create level set object
% First we create a level set object by four circles which initialize the
% zero level set.

% Create four circles to initialize the zero level set
[X,Y] = meshgrid(-35:35, -40:40);
phi = sqrt(X.^2 + Y.^2);
phi(phi > 25) = 25;
phi(phi < 24) = 24;
phi = phi - 24.5;
phi = repmat(phi, 2,2);

% Create level set object with narrowband of 16. We reinitialize to
% compute the signed distance function within the narrowband. The plot
% shows the zero level set as well as the level set function (green colors
% are positive, while red are negative)
LS = levelset2D(phi, 20);
LS = reinitialize(LS);
figure, plot(LS, 'contour', 'phi');


%% Propagate the level set function
% Propagate the level set function a given amount of time using the level
% set PDE specified in the speed_normal function:
% <speed_normal.html speed_normal.m>
time = 4;
for i = 1:10
    speed = (i <= 5) * 2 - 1;  % do dilation, followed by erosion
    LS = propagate(LS,time,'speed_normal', speed);
    LS = reinitialize(LS);
    figure;
    hold on;
    plot(LS);
    drawnow;
end
