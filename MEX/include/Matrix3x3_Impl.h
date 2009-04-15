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

#include <cmath>

template <typename DataType>
Matrix3x3<DataType>::Matrix3x3()
{
	m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; 
	m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; 
	m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; 
}

template <typename DataType>
Matrix3x3<DataType> Matrix3x3<DataType>::identity()
{
	return Matrix3x3();
}


template <typename DataType>
Matrix3x3<DataType> Matrix3x3<DataType>::operator*(const Matrix3x3& m2) const
{
	Matrix3x3 res;
	unsigned int i, j, k;
  
	for (i=0; i<3; i++)
	{
    for (j=0; j<3; j++)
    {
      res.m[i][j] = DataType(0.0);
      for (k=0; k<3; k++)
      {
		    res.m[i][j] += m[i][k] * m2.m[k][j];
      }
    }
	}
  
	return res;
}


template <typename DataType>
Matrix3x3<DataType> Matrix3x3<DataType>::operator-(const Matrix3x3& m2) const
{
	Matrix3x3 res;
	unsigned int i, j;
  
	for (i=0; i<3; i++)
	{
    for (j=0; j<3; j++)
    {
      res.m[i][j] = m[i][j] - m2.m[i][j];
    }
	}
  
	return res;
}


template <typename DataType>
Matrix3x3<DataType> Matrix3x3<DataType>::operator+(const Matrix3x3& m2) const
{
	Matrix3x3 res;
	unsigned int i, j;
  
	for (i=0; i<3; i++)
	{
    for (j=0; j<3; j++)
    {
      res.m[i][j] = m[i][j] + m2.m[i][j];
    }
	}
  
	return res;
}


// Input: i (row), j (column)
template <typename DataType>
DataType &Matrix3x3<DataType>::operator()(unsigned int i, unsigned int j)
{
	return m[i][j];
}

// Input: i (row), j (column)
template <typename DataType>
const DataType &Matrix3x3<DataType>::operator()(unsigned int i, unsigned int j) const
{
	return m[i][j];
}

