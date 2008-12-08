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
 * This class is a wrapper around MATLABs 2D array. Given a pointer, it gives
 * convenient access to elements using the (row,col) operator. Furthermore, you
 * can use an IndexArray to efficiently iterate only a set of elements.
 *
 *************************************************************************************************/

#ifndef _ARRAY_2D_H
#define _ARRAY_2D_H

#include <cmath>
#include "mex.h"
#include "IndexArray.h"

template <typename DataType, typename IndexType>
class Array2D
{
public :
    
    Array2D(const mxArray * A) {
        mxAssert(mxGetNumberOfDimensions(A) == 2, "Matrix of dimension 2 required");
        mA = static_cast<DataType *>(mxGetData(A));
        const mwSize * dimA = mxGetDimensions(A);
        mRows = dimA[0];
        mCols = dimA[1];
        mElements = mRows*mCols;
        //mexPrintf("Creating 2D-array of size %ix%i\n", mRows, mCols);
    }
    
    
    class IndexedIterator
    {
    public :

        IndexedIterator(Array2D * A, IndexArray<IndexType> * indices, IndexType index)
          : mA(A), mIndices(indices), mIndex(index), mI(0), mJ(0) {
              if (mIndices->length() > 0)
                  updateIndices();
        }

        inline bool operator != (const IndexedIterator & iter) const {
            return mIndex != iter.mIndex;
        }

        inline bool operator == (const IndexedIterator & iter) const {
            return mIndex == iter.mIndex;
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

        DataType & get(int i, int j) const {
            return (*mA)(mI+i, mJ+j);
        }

    protected :
        Array2D<DataType, IndexType> * mA;
        IndexArray<IndexType> * mIndices;
        IndexType mIndex;
        int mI, mJ;
        
        inline void updateIndices() {
            mA->ind2sub((*mIndices)[mIndex % mIndices->length()], mI, mJ);
        }

    };
    
    
    inline mwSize getRows() const { return mRows; }
    inline mwSize getCols() const { return mCols; }
    
    inline DataType & operator[](IndexType ind) const {
        mxAssert(ind < mElements, "Out of bounds in 2D-array");
        return mA[ind];
    }
    
    DataType & operator()(int i, int j) const {
        // Clamp indices to boundaries
        if (i < 0)              i = 0;
        else if (i > mRows-1)   i = mRows-1;
        
        if (j < 0)              j = 0;
        else if (j > mCols-1)   j = mCols-1;
        
        return (*this)[sub2ind(i,j)];
    }
    
    IndexedIterator begin(IndexArray<IndexType> & indices) {
        return IndexedIterator(this, &indices, 0);
    }
    IndexedIterator end(IndexArray<IndexType> & indices) {
        return IndexedIterator(this, &indices, indices.length());
    }
    
    inline void ind2sub(IndexType ind, int & i, int & j) const {
        i = ind % mRows;
        j = ind / mRows;
    }

    inline IndexType sub2ind(int i, int j) const {
        return j*mRows + i;
    }


protected :
    DataType * mA;
    mwSize mRows, mCols;
    IndexType mElements;

};

#endif
