
#include <cuda.h>
#include "ls_common.cu"
#include "ls_preparedata.cu"
#include "ls_sign.cu"
#include "ls_godunov.cu"

__global__ void reinitialize(
	const float * phi,
	float * out,
	const int * dim
	)
{
	__shared__ SHARED_DATA_DEF;

	int ind = 0;
	if (!prepareData(phi, SHARED_DATA, dim, ind))
		return;

	__syncthreads();

	float S = sign(SHARED_DATA);

	float dx2, dy2, dz2;
	godunov(SHARED_DATA, dx2, dy2, dz2, S);

	out[ind] = S *(1.0 - sqrtf(dx2 + dy2 + dz2));
}
