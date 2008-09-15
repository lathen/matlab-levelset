
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
%A = zeros(size(mask));
%dist = bwdist(mask);
%A(~mask) = dist(~mask);
%dist = -bwdist(~mask);
%A(mask) = dist(mask);
A = 5 - 10*mask; % alternative initialization only around zero-crossing

% Create levelset3D and plot
%LS = levelset3D(A,5);
%LS = levelset3D(A,5, 'Euler', 'FirstOrder', 'PDE');
LS = levelset3D(A,5, 'Euler', 'FirstOrder', 'FastMarching');
figure;
plot(LS, 'contour');

% Check largest gradient  magnitude
[Dx,Dy,Dz] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2 + Dz.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))

% Reinitialize
tic; LS = rebuild_narrowband(LS); toc
figure;
plot(LS, 'contour', 'narrowband 5');

% Check largest gradient  magnitude
[Dx,Dy,Dz] = diff_central(LS);
normgrad = sqrt(Dx.^2 + Dy.^2 + Dz.^2);
min_gradient = min(normgrad(:))
max_gradient = max(normgrad(:))
