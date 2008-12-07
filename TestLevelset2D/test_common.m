
%close all;
%clc;
clear;
addpath('..');

% Grid size
width = 200;
height = 200;

% Read shape from image
mask = imread('shape.gif');
mask = imresize(mask, [height width]);
mask = mask > 100; % make binary

% Compute signed distance transform (negative inside, positive outside)
A = 10 - 20*mask; % initialization only around zero-crossing

% Create levelset2D and plot
%LS = levelset2D(A);
%LS = levelset2D(A,10);
%LS = levelset2D(A,10, 'Euler', 'WENO');
%LS = levelset2D(A,Inf, 'Euler', 'FirstOrder', 'FastMarching');
LS = levelset2D(A,Inf, 'Euler', 'FirstOrder', 'FastSweeping');
tic; LS = rebuild_narrowband(LS); toc
figure;
plot(LS, 'contour', 'phi');

% Check largest gradient  magnitude
[Dx,Dy] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))
