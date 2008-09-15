
close all;
clear;
clc;
addpath('..');

% Grid size
width = 200;
height = 200;

% Read shape from image
mask = imread('shape.gif');
mask = imresize(mask, [height width]);
mask = mask > 100; % make binary

% Compute signed distance transform (negative inside, positive outside)
%A = zeros(size(mask));
%dist = bwdist(mask);
%A(~mask) = dist(~mask);
%dist = -bwdist(~mask);
%A(mask) = dist(mask);
A = 10 - 20*mask; % alternative initialization only around zero-crossing

% Create levelset2D and plot
%LS = levelset2D(A);
%LS = levelset2D(A,10);
%LS = levelset2D(A,10, 'Euler', 'WENO');
LS = levelset2D(A,10, 'Euler', 'FirstOrder', 'FastMarching');
figure;
plot(LS, 'contour', 'phi', 'gradient 5');

% Check largest gradient  magnitude
[Dx,Dy] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))

% Reinitialize
tic; LS = rebuild_narrowband(LS); toc
figure;
plot(LS, 'contour', 'phi', 'gradient 5');

% Check largest gradient  magnitude
[Dx,Dy] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))
