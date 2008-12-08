
addpath('..')
init;

% Propagate
time = 1;
for i = 1:30
    tic; [LS,iter] = propagate(LS,time,'speed_normal',1); toc
    tic; LS = reinitialize(LS); toc
    figure(100);
    clf;
    plot(LS, 'contour', 'narrowband 5');
    view(28,30);
    drawnow;
end

rmpath('..');
