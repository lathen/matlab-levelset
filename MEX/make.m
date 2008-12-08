
% Build levelset2D
disp('Building levelset2D');

disp('Compiling @levelset2D/private/diff_central_order2.cpp');
mex -I./include -outdir ../@levelset2D/private/ ../@levelset2D/private/diff_central_order2.cpp 

disp('Compiling @levelset2D/private/diff_upwind_order1.cpp');
mex -I./include -outdir ../@levelset2D/private/ ../@levelset2D/private/diff_upwind_order1.cpp 

disp('Compiling @levelset2D/private/diff_upwind_WENO.cpp');
mex -I./include -outdir ../@levelset2D/private/ ../@levelset2D/private/diff_upwind_WENO.cpp 

disp('Compiling @levelset2D/private/diff2_order2.cpp');
mex -I./include -outdir ../@levelset2D/private/ ../@levelset2D/private/diff2_order2.cpp

disp('Compiling @levelset2D/private/reinitialize_fastmarching.cpp');
mex -I./include -outdir ../@levelset2D/private/ ../@levelset2D/private/reinitialize_fastmarching.cpp 

disp('Compiling @levelset2D/private/reinitialize_fastsweeping.cpp');
mex -I./include -outdir ../@levelset2D/private/ ../@levelset2D/private/reinitialize_fastsweeping.cpp 


% Build levelset3D
disp('Building levelset3D');

disp('Compiling @levelset3D/private/diff_central_order2.cpp');
mex -I./include -outdir ../@levelset3D/private/ ../@levelset3D/private/diff_central_order2.cpp 

disp('Compiling @levelset3D/private/diff_upwind_order1.cpp');
mex -I./include -outdir ../@levelset3D/private/ ../@levelset3D/private/diff_upwind_order1.cpp 

disp('Compiling @levelset3D/private/diff2_order2.cpp');
mex -I./include -outdir ../@levelset3D/private/ ../@levelset3D/private/diff2_order2.cpp

disp('Compiling @levelset3D/private/min_curvature.cpp');
mex -I./include -outdir ../@levelset3D/private/ ../@levelset3D/private/min_curvature.cpp 

disp('Compiling @levelset3D/private/reinitialize_fastmarching.cpp');
mex -I./include -outdir ../@levelset3D/private/ ../@levelset3D/private/reinitialize_fastmarching.cpp 

disp('Compiling @levelset3D/private/triangulate.cpp');
mex -I./include -outdir ../@levelset3D/private/ ../@levelset3D/private/triangulate.cpp 
