
close all;
clear;
clc;
addpath('..');

% Grid size
width = 64;
height = 64;
depth = 64;

% Create test shape
mask = zeros(width,height,depth);
mask(width/2-width/6:width/2+width/6, ...
     height/2-height/6:height/2+height/6, ...
     depth/2-depth/6:depth/2+depth/6) = 1;
mask = mask==1;

narrowband = 5;
A = narrowband - 2*narrowband*mask; % create zero crossing to initialize level set function

% Create levelset3D and plot
%LS = levelset3D(A,narrowband);
%LS = levelset3D(A,narrowband, 'Euler', 'FirstOrder', 'PDE');
LS = levelset3D(A,narrowband, 'Euler', 'FirstOrder', 'FastMarching');
figure;
plot(LS, 'contour');

% Reinitialize
tic; LS = reinitialize(LS); toc
figure;
plot(LS, 'contour', 'narrowband 5');

% Check largest gradient  magnitude
[Dx,Dy,Dz] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2 + Dz.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))
