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
 * 2nd order second differences.
 *
 *************************************************************************************************/

#include "mex.h"
#include "Array3D.h"
#include "IndexArray.h"

template <typename DataType, typename IndexType>
void compute(mxArray * plhs[], const mxArray * prhs[])
{
    typedef Array3D<DataType, IndexType> MyArray3D;
    typedef IndexArray<IndexType> MyIndexArray;

    MyArray3D A(prhs[0]);
    MyIndexArray indices(prhs[1]);

    typename MyArray3D::IndexedIterator iter = A.begin(indices);
    typename MyArray3D::IndexedIterator iend = A.end(indices);

    plhs[0] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);
    plhs[3] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);
    plhs[4] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);
    plhs[5] = mxCreateNumericMatrix(1, indices.length(), mxGetClassID(prhs[0]), mxREAL);

    DataType * D2x = static_cast<DataType *>(mxGetData(plhs[0]));
    DataType * D2y = static_cast<DataType *>(mxGetData(plhs[1]));
    DataType * D2z = static_cast<DataType *>(mxGetData(plhs[2]));
    DataType * DxDy = static_cast<DataType *>(mxGetData(plhs[3]));
    DataType * DxDz = static_cast<DataType *>(mxGetData(plhs[4]));
    DataType * DyDz = static_cast<DataType *>(mxGetData(plhs[5]));

    IndexType ind = 0;
    while (iter != iend) {

        D2x[ind] = iter.get(0,1,0) - 2*iter.get(0,0,0) + iter.get(0,-1,0);
        D2y[ind] = iter.get(1,0,0) - 2*iter.get(0,0,0) + iter.get(-1,0,0);
        D2z[ind] = iter.get(0,0,1) - 2*iter.get(0,0,0) + iter.get(0,0,-1);
        
        DxDy[ind] = (-iter.get(1,1,0) + iter.get(-1,1,0) - iter.get(-1,-1,0) + iter.get(1,-1,0)) * 0.25;
        DxDz[ind] = (-iter.get(0,1,1) + iter.get(0,1,-1) - iter.get(0,-1,-1) + iter.get(0,-1,1)) * 0.25;
        DyDz[ind] = (-iter.get(1,0,1) + iter.get(1,0,-1) - iter.get(-1,0,-1) + iter.get(-1,0,1)) * 0.25;

        ind++;
        iter++;
    }    

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check for proper number of arguments
    if (nrhs != 2)
        mexErrMsgTxt("Two input arguments required.");
    if (nlhs != 6)
        mexErrMsgTxt("Six output arguments required.");

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

