/*************************************************************************************************
 *
 * Copyright (c) 2005
 * Michael Bang Nielsen (bang@daimi.au.dk)
 * University of Aarhus, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************
 * Contributors:  
 *                  1) Michael Bang Nielsen (bang@daimi.au.dk)
 *                  2) Ola Nilsson (ola.nilsson@itn.liu.se)
 *                  3) Gunnar Läthén (gunnar.lathe@itn.liu.se)
 *************************************************************************************************/

#ifndef _matrix3x3_h
#define _matrix3x3_h

template <typename DataType>
class Matrix3x3
  {
  public:
    Matrix3x3();
    
    // Input: i (row), j (column)
    const DataType & operator()(unsigned int i, unsigned int j) const;
    DataType & operator()(unsigned int i, unsigned int j);
    Matrix3x3 operator*(const Matrix3x3& m2) const;
    Matrix3x3 operator-(const Matrix3x3& m2) const;
    Matrix3x3 operator+(const Matrix3x3& m2) const;
    
    static Matrix3x3 identity();
    
  protected:
    DataType m[3][3];
  };

#include "Matrix3x3_Impl.h"

#endif
