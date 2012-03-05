
#ifndef _LS_DIFF2
#define _LS_DIFF2

#include "ls_common.cu"

__device__ void diff2(
	SHARED_DATA_DEF,
	float & d2x, float & d2y, float & d2z, float & dxdy, float & dxdz, float & dydz
	)
{
	d2x = NHOOD(1,0,0) - 2.0f*DATA + NHOOD(-1,0,0);
	d2y = NHOOD(0,1,0) - 2.0f*DATA + NHOOD(0,-1,0);
	d2z = NHOOD(0,0,1) - 2.0f*DATA + NHOOD(0,0,-1);

	dxdy = (-NHOOD(1,1,0) + NHOOD(1,-1,0) - NHOOD(-1,-1,0) + NHOOD(-1,1,0)) * 0.25f;
	dxdz = (-NHOOD(1,0,1) + NHOOD(1,0,-1) - NHOOD(-1,0,-1) + NHOOD(-1,0,1)) * 0.25f;
	dydz = (-NHOOD(0,1,1) + NHOOD(0,1,-1) - NHOOD(0,-1,-1) + NHOOD(0,-1,1)) * 0.25f;
}

#endif
