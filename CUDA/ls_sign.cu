
#ifndef _LS_SIGN
#define _LS_SIGN

#include "ls_common.cu"
#include "ls_diffcentral.cu"

__device__ float sign(SHARED_DATA_DEF)
{
	float dx, dy, dz;
	diffCentral(SHARED_DATA, dx, dy, dz);

	return DATA / sqrtf(DATA*DATA + dx*dx + dy*dy + dz*dz);
}

#endif
