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
 * This class is a wrapper around MATLABs 3D array. Given a pointer, it gives
 * convenient access to elements using the (row,col) operator. Furthermore, you
 * can use an IndexArray to efficiently iterate only a set of elements.
 *
 *************************************************************************************************/

#ifndef _ARRAY_3D_H
#define _ARRAY_3D_H

#include "mex.h"
#include "IndexArray.h"

template <typename DataType, typename IndexType>
class Array3D
{
public :
    
    Array3D(const mxArray * A) {
        mxAssert(mxGetNumberOfDimensions(A) == 3, "Matrix of dimension 3 required");
        mA = static_cast<DataType *>(mxGetData(A));
        const mwSize * dimA = mxGetDimensions(A);
        mRows = dimA[0];
        mCols = dimA[1];
        mSlices = dimA[2];
        mElements = mRows*mCols*mSlices;
        //mexPrintf("Creating 3D-array of size %ix%ix%i\n", mRows, mCols, mSlices);
    }
    
    
    class IndexedIterator
    {
    public :

        IndexedIterator(Array3D * A, IndexArray<IndexType> * indices, IndexType index)
          : mA(A), mIndices(indices), mIndex(index), mI(0), mJ(0), mK(0) {
              if (mIndices->length() > 0)
                  updateIndices();
        }

        inline bool operator != (const IndexedIterator & iter) const {
            return mIndex != iter.mIndex;
        }

        inline bool operator == (const IndexedIterator & iter) const {
            return mIndex != iter.mIndex;
        }

        inline DataType & operator*() const {
            return (*mA)[(*mIndices)[mIndex]];
        }

        inline void operator++(int) {
            mIndex++;
            updateIndices();
        }
        
        inline int getI() const { return mI; }
        inline int getJ() const { return mJ; }
        inline int getK() const { return mK; }

        DataType & get(int i, int j, int k) const {
            return (*mA)(mI+i, mJ+j, mK+k);
        }

    protected :
        Array3D<DataType, IndexType> * mA;
        IndexArray<IndexType> * mIndices;
        IndexType mIndex;
        int mI, mJ, mK;
        
        inline void updateIndices() {
            mA->ind2sub((*mIndices)[mIndex % mIndices->length()], mI, mJ, mK);
        }

    };
    
    
    inline mwSize getRows() const { return mRows; }
    inline mwSize getCols() const { return mCols; }
    inline mwSize getSlices() const { return mSlices; }
    
    inline DataType & operator[](IndexType ind) const {
        mxAssert(ind < mElements, "Out of bounds in 3D-array");
        return mA[ind];
    }
    
    DataType & operator()(int i, int j, int k) const {
        // Clamp indices to boundaries
        if (i < 0)              i = 0;
        else if (i > mRows-1)   i = mRows-1;
        
        if (j < 0)              j = 0;
        else if (j > mCols-1)   j = mCols-1;
        
        if (k < 0)              k = 0;
        else if (k > mSlices-1) k = mSlices-1;

        return (*this)[sub2ind(i,j,k)];
    }
    
    IndexedIterator begin(IndexArray<IndexType> & indices) {
        return IndexedIterator(this, &indices, 0);
    }
    IndexedIterator end(IndexArray<IndexType> & indices) {
        return IndexedIterator(this, &indices, indices.length());
    }
    
    inline void ind2sub(IndexType ind, int & i, int & j, int & k) const {
        i = ind % mRows;  ind /= mRows;
        j = ind % mCols;
        k = ind / mCols;
    }

    inline IndexType sub2ind(int i, int j, int k) const {
        return k*mCols*mRows + j*mRows + i;
    }


protected :
    DataType * mA;
    mwSize mRows, mCols, mSlices;
    IndexType mElements;

};

#endif
