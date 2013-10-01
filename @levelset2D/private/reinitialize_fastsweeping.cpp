/*************************************************************************************************
 *
 * Copyright (c) 2007
 * Ola Nilsson (ola.nilsson@itn.liu.se)
 * Linköping University, Sweden
 * All rights reserved.
 *
 *************************************************************************************************
 * Contributors:
 *                  1) Ola Nilsson (ola.nilsson@itn.liu.se)
 *                  2) Gunnar Läthén (gunnar.lathen@itn.liu.se)
 *************************************************************************************************/

#include "mex.h"
#include "Array2D.h"
#include "IndexArray.h"
#include <limits>
#include <algorithm>

template <typename DataType>
DataType minAbs(DataType a, DataType b) {
    return abs(a) < abs(b) ? a : b;
}


template <typename DataType>
inline DataType sign(DataType s)
{
	return s <= 0 ? -1 : 1;
}


template <typename DataType, typename IndexType>
void solve(Array2D<DataType, IndexType> & M, int i, int j)
{
    DataType uxmin, uymin, value, S;
    uxmin = minAbs(M(i-1, j), M(i+1, j));
    uymin = minAbs(M(i, j-1), M(i, j+1));
    S = sign(M(i,j));
    
    if(std::abs(uxmin-uymin) >= 1)
        value = minAbs(uxmin, uymin) + S;
    else
        value = (uxmin + uymin + S*std::sqrt(2 - (uxmin-uymin)*(uxmin-uymin)) )*0.5;
                 
    M(i,j) = minAbs(M(i,j), value);
}


/*
 * Interpolates new distance values for center point given
 * the four neighbors. Assumes center > 0.
 */
template <typename DataType>
DataType interpolateDistance(const DataType & center,
                             const DataType & im,
                             const DataType & ip,
                             const DataType & jm,
                             const DataType & jp)
{
  /* Approximate distance by linear interpolation along each
   grid axis */
  
   DataType b, val = std::numeric_limits<DataType>::max();
   if (im <= 0) {
   b = center / (center - im);
   if (val > b) val = b;
   }
   if (ip <= 0) {
   b = center / (center - ip);
   if (val > b) val = b;
   }
   if (jm <= 0) {
   b = center / (center - jm);
   if (val > b) val = b;
   }
   if (jp <= 0) {
   b = center / (center - jp);
   if (val > b) val = b;
   }
   return val;    
  
  /* Improved distance estimate by Adalsteinsson and Sethian in
   "The Fast Construction of Extension Velocities in Level Set Methods"
   (see fig. 4) */
  /*
  DataType s1 = std::numeric_limits<DataType>::max();
  DataType s2 = std::numeric_limits<DataType>::max();
  DataType t1 = std::numeric_limits<DataType>::max();
  DataType t2 = std::numeric_limits<DataType>::max();
  if (jm <= 0) s1 = center / (center - jm);
  if (jp <= 0) s2 = center / (center - jp);
  if (im <= 0) t1 = center / (center - im);
  if (ip <= 0) t2 = center / (center - ip);
  
  DataType s = std::min(s1,s2);
  DataType t = std::min(t1,t2);
  
  if (s == std::numeric_limits<DataType>::max())
    return t;
  
  if (t == std::numeric_limits<DataType>::max())
    return s;
  
  return s*t/std::sqrt(s*s + t*t); 
   */
}


/**
 * Doesn't handle boundary conditions.
 */
template <typename DataType, typename IndexType>
void compute(mxArray * plhs[], const mxArray * prhs[])
{
    typedef Array2D<DataType, IndexType> MyArray2D;
    typedef IndexArray<IndexType> MyIndexArray;
    
    MyArray2D A(prhs[0]);
    MyIndexArray indices(prhs[1]);
    
    
    DataType bandwidth = *(static_cast<DataType *>(mxGetData(prhs[2])));
    
    mwSize i, j;
    
    mxArray * array = mxCreateNumericMatrix(A.getRows(), A.getCols(), mxGetClassID(prhs[0]), mxREAL);
    MyArray2D B(array);
    
    
    // init with large values
    for (i = 0; i < A.getRows(); i++) {
        for (j = 0; j < A.getCols(); j++) {
            
            DataType im = A(i-1,j);
            DataType ip = A(i+1,j);
            DataType jm = A(i,j-1);
            DataType jp = A(i,j+1);
            DataType b, val;
            
            // If we have a zero-crossing, add the point to the accepted values
            // Flip the sign if we are inside
            if (A(i,j) > 0 && (im <= 0 || ip <= 0 || jm <= 0 || jp <= 0)) {
                B(i,j) = interpolateDistance(A(i,j), im, ip, jm, jp);
                //mexPrintf("Interpolated distance at (%i,%i) to (%f)\n", i,j, B(i,j));
            }
            else if (A(i,j) <= 0 && (im > 0 || ip > 0 || jm > 0 || jp > 0)) {
                B(i,j) = -interpolateDistance(-A(i,j), -im, -ip, -jm, -jp);
                //mexPrintf("Interpolated distance at (%i,%i) to (%f)\n", i,j, B(i,j));
            }
            else
                B(i,j) = std::numeric_limits<DataType>::max()*sign(A(i,j));
        }
    }
    
    // Sweep.
    // 1st sweep. i+j+
    for (i = 0; i < B.getRows(); i++)
        for (j = 0; j < B.getCols(); j++)
            if (std::abs(B(i, j)) > 1)
                solve(B, i, j);
    
    // 2nd sweep. i-j+
    for (i = B.getRows()-1; i >= 0; i--)
    	for (j = 0; j < B.getCols(); j++)
            if (std::abs(B(i, j)) > 1)
            	solve(B, i, j);
                
    // 3rd sweep. i-j-
    for (i = B.getRows()-1; i >= 0; i--)
    	for (j = B.getCols()-1; j >= 0; j--)
            if (std::abs(B(i, j)) > 1)
            	solve(B, i, j);
                
    // 4th sweep. i+j-
    for (i = 0; i < B.getRows(); i++)
    	for (j = B.getCols()-1; j >= 0; j--)
            if (std::abs(B(i, j)) > 1)
            	solve(B, i, j);


    IndexType count = 0;
    for (i = 0; i < B.getRows(); i++)
        for (j = 0; j < B.getCols(); j++)
            if (std::abs(B(i,j)) <= bandwidth)
                count++;

    // Create arrays for output
    plhs[0] = mxCreateNumericMatrix(1, count, mxGetClassID(prhs[0]), mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, count, mxGetClassID(prhs[1]), mxREAL);

    DataType * phi = static_cast<DataType *>(mxGetData(plhs[0]));
    IndexType * band = static_cast<IndexType *>(mxGetData(plhs[1]));
    
    count = 0;
    for (i = 0; i < B.getRows(); i++) {
        for (j = 0; j < B.getCols(); j++) {
            if (std::abs(B(i,j)) <= bandwidth) {
                phi[count] = B(i,j);
                band[count] = B.sub2ind(i,j) + 1;  // Matlab uses 1-based indexing
                count++;
            }
        }
    }
    
     mxDestroyArray(array);
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check for proper number of arguments
    if (nrhs != 3)
        mexErrMsgTxt("Three input arguments required.");
    if (nlhs != 2)
        mexErrMsgTxt("Two output arguments required.");
    
    // The input must be noncomplex
    if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) || mxIsComplex(prhs[2]) ||
    !mxIsNumeric(prhs[0]) || !mxIsNumeric(prhs[1]) || !mxIsNumeric(prhs[2]))
        mexErrMsgTxt("The input must be noncomplex (and numeric).");
    
    // The second input must be an unsigned integer
    if (mxIsSingle(prhs[1]) || mxIsDouble(prhs[1]))
        mexErrMsgTxt("The second input must be an integer.");
    
    if (mxIsSingle(prhs[0])) {
        if (mxIsInt8(prhs[1]))
            compute<float, char>(plhs, prhs);
        else if (mxIsUint8(prhs[1]))
            compute<float, unsigned char>(plhs, prhs);
        else if (mxIsInt16(prhs[1]))
            compute<float, short>(plhs, prhs);
        else if (mxIsUint16(prhs[1]))
            compute<float, unsigned short>(plhs, prhs);
        else if (mxIsInt32(prhs[1]))
            compute<float, int>(plhs, prhs);
        else if (mxIsUint32(prhs[1]))
            compute<float, unsigned int>(plhs, prhs);
        else if (mxIsInt64(prhs[1]))
            compute<float, long long int>(plhs, prhs);
        else if (mxIsUint64(prhs[1]))
            compute<float, unsigned long long int>(plhs, prhs);
    }
    else if (mxIsDouble(prhs[0])) {
        if (mxIsInt8(prhs[1]))
            compute<double, char>(plhs, prhs);
        else if (mxIsUint8(prhs[1]))
            compute<double, unsigned char>(plhs, prhs);
        else if (mxIsInt16(prhs[1]))
            compute<double, short>(plhs, prhs);
        else if (mxIsUint16(prhs[1]))
            compute<double, unsigned short>(plhs, prhs);
        else if (mxIsInt32(prhs[1]))
            compute<double, int>(plhs, prhs);
        else if (mxIsUint32(prhs[1]))
            compute<double, unsigned int>(plhs, prhs);
        else if (mxIsInt64(prhs[1]))
            compute<double, long long int>(plhs, prhs);
        else if (mxIsUint64(prhs[1]))
            compute<double, unsigned long long int>(plhs, prhs);
    }
    else
        mexErrMsgTxt("Input types not supported.");
}

