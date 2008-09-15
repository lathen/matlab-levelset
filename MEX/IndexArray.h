

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
