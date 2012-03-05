
#ifndef _LS_DIFFCENTRAL
#define _LS_DIFFCENTRAL

#include "ls_common.cu"

__device__ void diffCentral(
	SHARED_DATA_DEF,
	float & dx, float & dy, float & dz
	)
{
	dx = (NHOOD(1,0,0) - NHOOD(-1,0,0))*0.5;
	dy = (NHOOD(0,1,0) - NHOOD(0,-1,0))*0.5;
	dz = (NHOOD(0,0,1) - NHOOD(0,0,-1))*0.5;
}

#endif
