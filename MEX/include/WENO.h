/*************************************************************************************************
 *
 * Copyright (c) 2008
 * Gunnar Läthén (gunnar.lathen@itn.liu.se)
 * Linköping University, Sweden.
 *
 *************************************************************************************************
 * Contributors:  
 *                  1) Gunnar Läthén (gunnar.lathen@itn.liu.se)
 *************************************************************************************************
 *
 * This class implements the WENO upwind scheme as presented in
 *  "S. Osher and R. Fedkiw. Level Set and Dynamic Implicit Surfaces.
 *   Springer-Verlag New York Inc., 2003."
 *
 *************************************************************************************************/

#ifndef _WENO_H
#define _WENO_H

#include <cmath>
#include <limits>

template <typename DataType>
DataType WENO(DataType v1,
              DataType v2,
              DataType v3,
              DataType v4,
              DataType v5)
{
    DataType forth = 1.0/4.0;
    DataType sixth = 1.0/6.0;
    DataType twelvth = 1.0/12.0;
	
	// Define the three potential HJ ENO approximations to df
    DataType df_1 = (2.0*v1 - 7.0*v2 + 11.0*v3) * sixth;
    DataType df_2 = (   -v2 + 5.0*v3 +  2.0*v4) * sixth;
    DataType df_3 = (2.0*v3 + 5.0*v4 -      v5) * sixth;
	
	// Estimate the smoothness of the stencils for the
	// three potential approximations
	DataType S1, S2, S3;
	DataType S_1, S_2;

    S_1 = v1 - 2.0*v2 + v3;	S_2 = v1 - 4.0*v2 + 3.0*v3;
    S1 = 13.0*twelvth * S_1*S_1 + forth * S_2*S_2;

	S_1 = v2 - 2.0*v3 + v4;	S_2 = v2 - v4;
	S2 = 13.0*twelvth * S_1*S_1 + forth * S_2*S_2;

	S_1 = v3 - 2.0*v4 + v5; S_2 = 3.0*v3 - 4.0*v4 + v5;
    S3 = 13.0*twelvth * S_1*S_1 + forth * S_2*S_2;
	
	DataType eps = 0.000001;

	// Calculate blending weights
    DataType a1 = 0.1 / ((S1 + eps)*(S1 + eps));
    DataType a2 = 0.6 / ((S2 + eps)*(S2 + eps));
    DataType a3 = 0.3 / ((S3 + eps)*(S3 + eps));

	// Return the polynomial approximation
	DataType norm = a1 + a2 + a3;
	return (a1*df_1 + a2*df_2 + a3*df_3) / norm;
}

#endif
