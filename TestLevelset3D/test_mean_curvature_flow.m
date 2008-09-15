
test_common;

% Propagate
time = 1;
for i = 1:20
    tic; [LS,iter] = propagate(LS,time,'mean_curvature_flow',1); toc
    tic; LS = rebuild_narrowband(LS); toc
    figure;
    plot(LS, 'contour', 'gradient 5');
    drawnow;
end

rmpath('..');
