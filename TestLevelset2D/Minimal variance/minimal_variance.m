
addpath('..');
init;

% Read input image for segmentation
A = imread('image.jpg');
A = imresize(A, size(LS));
figure;
imagesc(A);
colormap(gray);

% Propagate (until convergence)
time = 0.002;
hold on;
for i = 1:20
    tic; [LS,iter] = propagate(LS,time,'minimal_variance_operator',A); toc
    tic; LS = reinitialize(LS); toc

    clf;
    imagesc(A);
    plot(LS, 'contour');
    drawnow;
end

rmpath('..');
