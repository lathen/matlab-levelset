
#ifndef _LS_PREPAREDATA
#define _LS_PREPAREDATA

#include "ls_common.cu"

__device__ bool prepareData(
	const float * data,
	SHARED_DATA_DEF,
	const int * dim,
	int & ind
	)
{
    // Separate blockIdx.yz
    const int blockIdx_y = blockIdx.y % (dim[1] / BLOCKDIM_Y);
    const int blockIdx_z = blockIdx.y / (dim[1] / BLOCKDIM_Y);

	// Compute grid coordinates
    const int baseX = blockIdx.x * BLOCKDIM_X - PADDING + threadIdx.x;
    const int baseY = blockIdx_y * BLOCKDIM_Y - PADDING + threadIdx.y;
    const int baseZ = blockIdx_z * BLOCKDIM_Z - PADDING + threadIdx.z;

	// Compute linear index to grid point
	ind = (baseZ*dim[1] + baseY)*dim[0] + baseX;

	// Load the data to the shared memory
	if (baseX >= 0 && baseX < dim[0] &&
		baseY >= 0 && baseY < dim[1] &&
		baseZ >= 0 && baseZ < dim[2])
		DATA = data[ind];
	else
		DATA = 0;

	// Return false if the thread is in the padded region
	if (threadIdx.x < PADDING || threadIdx.x >= BLOCKDIM_X+PADDING ||
		threadIdx.y < PADDING || threadIdx.y >= BLOCKDIM_Y+PADDING ||
		threadIdx.z < PADDING || threadIdx.z >= BLOCKDIM_Z+PADDING)
		return false;

	return true;
}

#endif
