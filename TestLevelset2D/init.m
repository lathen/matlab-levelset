
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

% Create zero crossing to initialize the level set function
narrowband = 20;
A = narrowband - 2*narrowband*mask;

% Create levelset2D and plot
%LS = levelset2D(A);
%LS = levelset2D(A,narrowband);
%LS = levelset2D(A,narrowband, 'Euler', 'WENO');
%LS = levelset2D(A,Inf, 'Euler', 'FirstOrder', 'FastMarching');
LS = levelset2D(A,narrowband, 'Euler', 'FirstOrder', 'FastSweeping');
tic; LS = reinitialize(LS); toc
figure;
plot(LS, 'contour', 'phi');

% Check largest gradient  magnitude
[Dx,Dy] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))
