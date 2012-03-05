
#ifndef _LS_MEANCURVATURE
#define _LS_MEANCURVATURE

#include "ls_common.cu"
#include "ls_diffcentral.cu"
#include "ls_diff2.cu"

__device__ void meanCurvature(
	SHARED_DATA_DEF,
	float & c
	)
{
	float dx, dy, dz;
	diffCentral(SHARED_DATA, dx, dy, dz);

	float dx2 = dx*dx;
	float dy2 = dy*dy;
	float dz2 = dz*dz;
	float denominator = powf(2.0f*(dx2 + dy2 + dz2), 1.5f);

	if (denominator < 0.00035) { // Roughly sqrt(float epsilon)
		c = 0;
		return;
	}

	float d2x, d2y, d2z, dxdy, dxdz, dydz;
	diff2(SHARED_DATA, d2x, d2y, d2z, dxdy, dxdz, dydz);


	float nominator = dx2*(d2y + d2z) - 2.0f*dy*dz*dydz +
					  dy2*(d2x + d2z) - 2.0f*dx*dz*dxdz +
					  dz2*(d2x + d2y) - 2.0f*dx*dy*dxdy;
	c = nominator / denominator;
}

#endif
