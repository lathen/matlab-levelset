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
 * 3rd-5th order finite differences according to WENO scheme.
 *
 *************************************************************************************************/

#include "mex.h"
#include "Array2D.h"
#include "IndexArray.h"
#include "WENO.h"

template <typename DataType, typename IndexType>
void compute(mxArray * plhs[], const mxArray * prhs[])
{
    typedef Array2D<DataType, IndexType> MyArray2D;
    typedef IndexArray<IndexType> MyIndexArray;

    MyArray2D A(prhs[0]);
    MyIndexArray indices(prhs[1]);

    typename MyArray2D::IndexedIterator iter = A.begin(indices);
    typename MyArray2D::IndexedIterator iend = A.end(indices);

    plhs[0] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);
    plhs[3] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);

    DataType * Dxm = static_cast<DataType *>(mxGetData(plhs[0]));
    DataType * Dxp = static_cast<DataType *>(mxGetData(plhs[1]));
    DataType * Dym = static_cast<DataType *>(mxGetData(plhs[2]));
    DataType * Dyp = static_cast<DataType *>(mxGetData(plhs[3]));

    IndexType ind = 0;
    DataType v1, v2, v3, v4, v5;
    while (iter != iend) {

        //int i = iter.getI();
        //int j = iter.getJ();
        //mexPrintf("Visiting (%i,%i)\n", i, j);
        
        v1 = iter.get(0,-2) - iter.get(0,-3);
        v2 = iter.get(0,-1) - iter.get(0,-2);
        v3 = iter.get(0, 0) - iter.get(0,-1);
        v4 = iter.get(0, 1) - iter.get(0, 0);
        v5 = iter.get(0, 2) - iter.get(0, 1);
        Dxm[ind] = WENO(v1,v2,v3,v4,v5);
        
        v1 = iter.get(0, 3) - iter.get(0, 2);
        v2 = iter.get(0, 2) - iter.get(0, 1);
        v3 = iter.get(0, 1) - iter.get(0, 0);
        v4 = iter.get(0, 0) - iter.get(0,-1);
        v5 = iter.get(0,-1) - iter.get(0,-2);
        Dxp[ind] = WENO(v1,v2,v3,v4,v5);

        v1 = iter.get(-2,0) - iter.get(-3,0);
        v2 = iter.get(-1,0) - iter.get(-2,0);
        v3 = iter.get( 0,0) - iter.get(-1,0);
        v4 = iter.get( 1,0) - iter.get( 0,0);
        v5 = iter.get( 2,0) - iter.get( 1,0);
        Dym[ind] = WENO(v1,v2,v3,v4,v5);
        
        v1 = iter.get( 3,0) - iter.get( 2,0);
        v2 = iter.get( 2,0) - iter.get( 1,0);
        v3 = iter.get( 1,0) - iter.get( 0,0);
        v4 = iter.get( 0,0) - iter.get(-1,0);
        v5 = iter.get(-1,0) - iter.get(-2,0);
        Dyp[ind] = WENO(v1,v2,v3,v4,v5);

        ind++;
        iter++;
    }    

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check for proper number of arguments
    if (nrhs != 2)
        mexErrMsgTxt("Two input arguments required.");
    if (nlhs != 4)
        mexErrMsgTxt("Four output arguments required.");
    
    // The input must be noncomplex
    if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) ||
        !mxIsNumeric(prhs[0]) || !mxIsNumeric(prhs[1]))
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

