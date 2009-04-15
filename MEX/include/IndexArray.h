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
 * This class is used in conjunction with Array2D and Array3D to efficiently
 * iterate a set of elements in a MATLAB array.
 *
 *************************************************************************************************/

#ifndef _INDEX_ARRAY_H
#define _INDEX_ARRAY_H

#include "mex.h"

template <typename IndexType>
class IndexArray
{
public :
    
    IndexArray(const mxArray * indices) {
        mIndices = static_cast<IndexType *>(mxGetData(indices));
        mElements = mxGetNumberOfElements(indices);
        //mexPrintf("Creating index array of length %i\n", mElements);
    }
    
    inline IndexType operator[](IndexType ind) {
        mxAssert(ind < mElements, "Out of bounds in index array");
        return mIndices[ind]-1;
    }
    
    inline IndexType length() const { return mElements; }
    
protected :
    IndexType * mIndices;
    IndexType mElements;
};

#endif
