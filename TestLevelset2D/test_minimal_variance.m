
test_common;

% Read input image for segmentation
A = imread('image.jpg');
A = imresize(A, size(LS));
figure;
imagesc(A);
colormap(gray);

% Propagate (until convergence)
time = 0.005;
hold on;
for i = 1:20
    tic; [LS,iter] = propagate(LS,time,'minimal_variance',A); toc
    tic; LS = rebuild_narrowband(LS); toc

    imagesc(A);
    plot(LS, 'contour', 'narrowband');
    drawnow;
end

imagesc(A);
plot(LS, 'contour');

rmpath('..');
