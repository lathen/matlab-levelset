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
 * Computes the two principal curvatures of a surface.
 *
 *************************************************************************************************/

#include "mex.h"
#include "Array3D.h"
#include "IndexArray.h"
#include "Matrix3x3.h"
#include <limits>
#include <algorithm>

template <typename DataType>
class Vector
  {
    public :
    Vector(DataType x = 0, DataType y = 0, DataType z = 0) {
      mData[0] = x;
      mData[1] = y;
      mData[2] = z;
    }
    
    inline const DataType & operator [] (unsigned int ind) const {
      mxAssert(ind < 3, "Vector out-of-bounds access");
      return mData[ind];
    }
    
    inline DataType & operator [] (unsigned int ind) {
      mxAssert(ind < 3, "Vector out-of-bounds access");
      return mData[ind];
    }
    
    inline DataType length() const { return std::sqrt(mData[0]*mData[0] + mData[1]*mData[1] + mData[2]*mData[2]); }
    
    void normalize() {
      DataType len = length();
      if (len < std::sqrt(std::numeric_limits<DataType>::epsilon())) return;
      mData[0] /= len;
      mData[1] /= len;
      mData[2] /= len;
    }
    
    Vector operator * (const DataType & val) const {
      return Vector(mData[0]*val, mData[1]*val, mData[2]*val);
    }
    
    Vector operator + (const Vector & vec) const {
      return Vector(mData[0]+vec[0], mData[1]+vec[1], mData[2]+vec[2]);
    }
    
    protected :
    DataType mData[3];
  };



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
  
  DataType * c1 = static_cast<DataType *>(mxGetData(plhs[0]));
    DataType * c2 = static_cast<DataType *>(mxGetData(plhs[1]));
  
  IndexType ind = 0;
  while (iter != iend) {
    
    //int i = iter.getI();
    //int j = iter.getJ();
    //int k = iter.getK();
    //mexPrintf("Visiting (%i,%i,%i): %f\n", i, j, k, *iter);
    
    // Compute normals in neighborhood
    
    // (x+0.5,y,z)
    Vector<DataType> Nxp;
    Nxp[0] = iter.get(1,0,0) - *iter;
    Nxp[1] = (iter.get(1,1,0) - iter.get(1,-1,0) + iter.get(0,1,0) - iter.get(0,-1,0)) * 0.25;
    Nxp[2] = (iter.get(1,0,1) - iter.get(1,0,-1) + iter.get(0,0,1) - iter.get(0,0,-1)) * 0.25;
    Nxp.normalize();
    
    // (x-0.5,y,z)
    Vector<DataType> Nxm;
    Nxm[0] = *iter - iter.get(-1,0,0);
    Nxm[1] = (iter.get(-1,1,0) - iter.get(-1,-1,0) + iter.get(0,1,0) - iter.get(0,-1,0)) * 0.25;
    Nxm[2] = (iter.get(-1,0,1) - iter.get(-1,0,-1) + iter.get(0,0,1) - iter.get(0,0,-1)) * 0.25;
    Nxm.normalize();
    
    // (x,y+0.5,z)
    Vector<DataType> Nyp;
    Nyp[0] = (iter.get(1,1,0) - iter.get(-1,1,0) + iter.get(1,0,0) - iter.get(-1,0,0)) * 0.25;
    Nyp[1] = iter.get(0,1,0) - *iter;
    Nyp[2] = (iter.get(0,1,1) - iter.get(0,1,-1) + iter.get(0,0,1) - iter.get(0,0,-1)) * 0.25;
    Nyp.normalize();
    
    // (x,y-0.5,z)
    Vector<DataType> Nym;
    Nym[0] = (iter.get(1,-1,0) - iter.get(-1,-1,0) + iter.get(1,0,0) - iter.get(-1,0,0)) * 0.25;
    Nym[1] = *iter - iter.get(0,-1,0);
    Nym[2] = (iter.get(0,-1,1) - iter.get(0,-1,-1) + iter.get(0,0,1) - iter.get(0,0,-1)) * 0.25;
    Nym.normalize();
    
    // (x,y,z+0.5)
    Vector<DataType> Nzp;
    Nzp[0] = (iter.get(1,0,1) - iter.get(-1,0,1) + iter.get(1,0,0) - iter.get(-1,0,0)) * 0.25;
    Nzp[1] = (iter.get(0,1,1) - iter.get(0,-1,1) + iter.get(0,1,0) - iter.get(0,-1,0)) * 0.25;
    Nzp[2] = iter.get(0,0,1) - *iter;
    Nzp.normalize();
    
    // (x,y,z-0.5)
    Vector<DataType> Nzm;
    Nzm[0] = (iter.get(1,0,-1) - iter.get(-1,0,-1) + iter.get(1,0,0) - iter.get(-1,0,0)) * 0.25;
    Nzm[1] = (iter.get(0,1,-1) - iter.get(0,-1,-1) + iter.get(0,1,0) - iter.get(0,-1,0)) * 0.25;
    Nzm[2] = *iter - iter.get(0,0,1);
    Nzm.normalize();
    
    // Compute normal derivative matrix
    Matrix3x3<DataType> N;
    
    // along x
    N(0,0) = Nxp[0] - Nxm[0];
    N(0,1) = Nxp[1] - Nxm[1];
    N(0,2) = Nxp[2] - Nxm[2];
    
    // along y
    N(1,0) = Nyp[0] - Nym[0];
    N(1,1) = Nyp[1] - Nym[1];
    N(1,2) = Nyp[2] - Nym[2];
    
    // along z
    N(2,0) = Nzp[0] - Nzm[0];
    N(2,1) = Nzp[1] - Nzm[1];
    N(2,2) = Nzp[2] - Nzm[2];
    
    // Compute normal at (x,y,z)
    Vector<DataType> n;
    n[0] = (iter.get(1,0,0) - iter.get(-1,0,0)) * 0.5;
    n[1] = (iter.get(0,1,0) - iter.get(0,-1,0)) * 0.5;
    n[2] = (iter.get(0,0,1) - iter.get(0,0,-1)) * 0.5;
    n.normalize();
    
    // Compute tensor (outer) product of n
    Matrix3x3<DataType> nxn;
    nxn(0,0) = n[0]*n[0];  nxn(0,1) = n[1]*n[0];  nxn(0,2) = n[2]*n[0];
    nxn(1,0) = n[0]*n[1];  nxn(1,1) = n[1]*n[1];  nxn(1,2) = n[2]*n[1];
    nxn(2,0) = n[0]*n[2];  nxn(2,1) = n[1]*n[2];  nxn(2,2) = n[2]*n[2];
    
    // Compute shape matrix
    Matrix3x3<DataType> B = N*(Matrix3x3<DataType>::identity() - nxn);
    
    // Compute the eigenvalues of B
    // x^3 + ax^2 + bx + c = 0
    DataType a,b;
    a = -(B(0,0) + B(1,1) + B(2,2));
    b = -(B(1,0)*B(0,1) + B(2,0)*B(0,2) + B(1,2)*B(2,1) - B(0,0)*B(1,1) - B(0,0)*B(2,2) - B(1,1)*B(2,2));
    
    DataType k1 = 0, k2 = 0, root;
    root = a*a - 4.0*b;
    if (root < 0)
      root = 0;
    k1 = 0.5*(-a + std::sqrt(root));
    k2 = 0.5*(-a - std::sqrt(root));
    
    // Mean curvature is 1/2 * Tr(B)
    //DataType K = 0.5*(B(0,0) + B(1,1) + B(2,2));
    //mexPrintf("Mean curvature: %f\n", K);

    //mexPrintf("Principal curvatures: %f, %f\n", k1, k2);
    c1[ind] = std::min(k1,k2);
    c2[ind] = std::max(k1,k2);      
    
    ind++;
    iter++;
  }    
  
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check for proper number of arguments
  if (nrhs != 2)
    mexErrMsgTxt("Two input arguments required.");
  if (nlhs != 2)
    mexErrMsgTxt("Two output arguments required.");
  
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

