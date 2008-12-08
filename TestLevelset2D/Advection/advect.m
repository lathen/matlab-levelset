
addpath('..');
init;

% Create vector field (note that origin is in upper left corner)
s = size(LS);
rows = s(1); cols = s(2);
Fx = ones(rows,cols);
Fy = -ones(rows,cols);

figure;

% Propagate
time = 1;
for i = 1:30
    tic; [LS,iter] = propagate(LS,time,'advect_operator', Fx,Fy); toc
    tic; LS = reinitialize(LS); toc
    clf;
    plot(LS, 'contour', 'phi', 'gradient 5');
    drawnow;
end

d = 10;
quiver([1:d:cols], [1:d:rows], Fx(1:d:end,1:d:end), Fy(1:d:end,1:d:end), 'g');

rmpath('..');
