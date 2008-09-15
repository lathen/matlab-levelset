
test_common;

% Propagate
time = 5;
for i = 1:20
    tic; [LS,iter] = propagate(LS,time,'curvature_flow',1); toc
    tic; LS = rebuild_narrowband(LS); toc
    figure(100);
    plot(LS, 'contour', 'phi', 'gradient 5', 'narrowband');
    drawnow;
end

rmpath('..');
