
#ifndef _LS_GODONUV
#define _LS_GODONUV

#include "ls_common.cu"
#include "ls_diffonesided.cu"

__device__ void godunov(
	SHARED_DATA_DEF,
	float & dx2, float & dy2, float & dz2, float a
	)
{
	float dxm, dxp, dym, dyp, dzm, dzp;
	diffOnesided(SHARED_DATA, dxm, dxp, dym, dyp, dzm, dzp);

	if (a > 0) {
		dx2 = max( powf(max(dxm,0.0f),2), powf(min(dxp,0.0f),2) );
		dy2 = max( powf(max(dym,0.0f),2), powf(min(dyp,0.0f),2) );
		dz2 = max( powf(max(dzm,0.0f),2), powf(min(dzp,0.0f),2) );
	}
	else {
		dx2 = max( powf(min(dxm,0.0f),2), powf(max(dxp,0.0f),2) );
		dy2 = max( powf(min(dym,0.0f),2), powf(max(dyp,0.0f),2) );
		dz2 = max( powf(min(dzm,0.0f),2), powf(max(dzp,0.0f),2) );
	}
}

#endif
