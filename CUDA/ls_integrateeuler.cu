
#include <cuda.h>

__global__ void integrateEuler(
	float * phi,
	const float * ddt,
	const float dt,
	const int * dim
	)
{
    // Separate blockIdx.yz
    const int blockIdx_y = blockIdx.y % (dim[1] / BLOCKDIM_Y);
    const int blockIdx_z = blockIdx.y / (dim[1] / BLOCKDIM_Y);

	// Compute grid coordinates
    const int baseX = blockIdx.x * BLOCKDIM_X + threadIdx.x;
    const int baseY = blockIdx_y * BLOCKDIM_Y + threadIdx.y;
    const int baseZ = blockIdx_z * BLOCKDIM_Z + threadIdx.z;

	// Compute linear index to grid point
	int ind = (baseZ*dim[1] + baseY)*dim[0] + baseX;

	// Propagate
    phi[ind] += dt*ddt[ind];
}
