
addpath('..')
init;

% Create vector field (note that origin is in upper left corner)
s = size(LS);
rows = s(1); cols = s(2); slices = s(3);
Fx = ones(rows,cols,slices);
Fy = -ones(rows,cols,slices);
Fz = ones(rows,cols,slices);

figure;

% Propagate
time = 1;
for i = 1:20
    tic; [LS,iter] = propagate(LS,time,'advect_operator', Fx,Fy,Fz); toc
    tic; LS = reinitialize(LS); toc
    clf;
    plot(LS, 'contour', 'gradient 5');
    view(28,30);
    drawnow;
end

%d = 10;
%quiver([1:d:cols], [1:d:rows], Fx(1:d:end,1:d:end), Fy(1:d:end,1:d:end), 'g');

rmpath('..');
