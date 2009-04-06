%% Minimal variance
% Example of image segmentation by "minimal variance" as proposed by Chan
% and Vese in:
% T. Chan and L. Vese. Active contours without edges. IEEE Transactions on
% Image Processing, 10(2):266?277, February 2001.

%% Create level set object
% First we create a level set object by circles which initialize the
% zero level set.

% Create circles to initialize the zero level set
[X,Y] = meshgrid(-10:10, -10:10);
phi = sqrt(X.^2 + Y.^2);
phi(phi > 5) = 5;
phi(phi < 4) = 4;
phi = phi - 4.5;
phi = repmat(phi, 10,10);

% Create level set object and reinitialize to ompute the signed distance
% function. The plot shows the zero level set as well as the level set
% function (green colors re positive, while red are negative)
LS = levelset2D(phi);
LS = reinitialize(LS);
figure, plot(LS, 'contour', 'phi');


%% Read input image for segmentation
A = imread('image.jpg');
A = imresize(A, size(LS));
figure;
imshow(A, []);
hold on;
plot(LS);


%% Propagate the level set function
% Propagate the level set function a given amount of time using the level
% set PDE specified in the minimal_variance_operator function:
% <minimal_variance_operator.html minimal_variance_operator.m>
figure;
imshow(A, []);
hold on;
plot(LS);
time = 0.005;
for i = 1:5
    LS = propagate(LS,time,'minimal_variance_operator', A, 0.5);
    LS = reinitialize(LS);

    figure;
    imshow(A, []);
    hold on;
    plot(LS);
    drawnow;
end
