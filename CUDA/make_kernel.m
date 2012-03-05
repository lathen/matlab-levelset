function k = make_kernel(kernel, vol, blockDim, padding)

path = fileparts(mfilename('fullpath'));

delete([kernel,'.ptx']);

% Compile kernel
nvcc('-ptx', ...
     ['-D"BLOCKDIM_X=',int2str(blockDim(1)),'"'], ...
     ['-D"BLOCKDIM_Y=',int2str(blockDim(2)),'"'], ...
     ['-D"BLOCKDIM_Z=',int2str(blockDim(3)),'"'], ...
     ['-D"PADDING=',int2str(padding),'"'], ...
     ['"',path,'/',kernel,'.cu"']);

% Create kernel object and set thread and grid sizes
k = parallel.gpu.CUDAKernel([kernel,'.ptx'], [kernel,'.cu']);
k.ThreadBlockSize = blockDim+2*padding;

[rows,cols,slices] = size(vol);
k.GridSize = [rows/blockDim(1) cols*slices/(blockDim(2)*blockDim(3))];
