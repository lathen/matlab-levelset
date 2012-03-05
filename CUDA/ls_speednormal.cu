
#include <cuda.h>
#include "ls_common.cu"
#include "ls_preparedata.cu"
#include "ls_godunov.cu"
#include "ls_meancurvature.cu"

__global__ void speedNormal(
	const float * phi,
	float * out,
	const int * dim,
	const float * F,
	const float a
	)
{
	__shared__ SHARED_DATA_DEF;

	int ind = 0;
	if (!prepareData(phi, SHARED_DATA, dim, ind))
		return;

	__syncthreads();

	float dx2, dy2, dz2;
	godunov(SHARED_DATA, dx2, dy2, dz2, F[ind]);

	float dx, dy, dz;
	diffCentral(SHARED_DATA, dx, dy, dz);

	float c;
	meanCurvature(SHARED_DATA, c);

	out[ind] = -F[ind]*sqrtf(dx2 + dy2 + dz2) + a*c*sqrtf(dx*dx + dy*dy + dz*dz);
}
