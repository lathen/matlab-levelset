
#ifndef _LS_DIFFONESIDED
#define _LS_DIFFONESIDED

#include "ls_common.cu"

__device__ void diffOnesided(
	SHARED_DATA_DEF,
	float & dxm, float & dxp,
	float & dym, float & dyp,
	float & dzm, float & dzp
	)
{
	dxm = DATA - NHOOD(-1,0,0);
	dxp = NHOOD(1,0,0) - DATA;

	dym = DATA - NHOOD(0,-1,0);
	dyp = NHOOD(0,1,0) - DATA;

	dzm = DATA - NHOOD(0,0,-1);
	dzp = NHOOD(0,0,1) - DATA;
}

#endif
