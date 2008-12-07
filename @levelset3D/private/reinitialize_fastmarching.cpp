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
 * Fast marching implementation for reinitialization.
 *
 *************************************************************************************************/

#include "mex.h"
#include "Array3D.h"
#include "IndexArray.h"
#include "Heap.h"
#include <map>
#include <cmath>

struct Coordinate
{
    Coordinate(int I, int J, int K) : i(I), j(J), k(K), inside(false) { }
    Coordinate(int I, int J, int K, bool Inside) : i(I), j(J), k(K), inside(Inside) { }
    int i, j, k;
    bool inside;
    
    bool operator < (const Coordinate & coord) const {
        return                  i < coord.i  ||
               (i == coord.i && j < coord.j) ||
               (i == coord.i && j == coord.j && k < coord.k);
    }
};

template <typename DataType>
class Node : public Heap<DataType>::Heapable
{
public :
	Node(Coordinate coord) : coordinate(coord) { }
    Node(int i, int j, int k) : coordinate(i,j,k) { }
	Coordinate coordinate;
};


template <typename DataType>
DataType computeTentativeValue(int i, int j, int k,
                               std::map<Coordinate, DataType> & accepted)
{
    const DataType inf = std::numeric_limits<DataType>::max();
    typename std::map<Coordinate, DataType>::const_iterator iter;

    // Fetch (i-1,j,k)   (returns infinity if not accepted)
    iter = accepted.find(Coordinate(i-1,j,k));
    DataType im = (iter == accepted.end() ? inf : (*iter).second);

    // Fetch (i+1,j,k)
    iter = accepted.find(Coordinate(i+1,j,k));
    DataType ip = (iter == accepted.end() ? inf : (*iter).second);

    // Fetch (i,j-1,k)
    iter = accepted.find(Coordinate(i,j-1,k));
    DataType jm = (iter == accepted.end() ? inf : (*iter).second);

    // Fetch (i,j+1,k)
    iter = accepted.find(Coordinate(i,j+1,k));
    DataType jp = (iter == accepted.end() ? inf : (*iter).second);
    
    // Fetch (i,j,k-1)
    iter = accepted.find(Coordinate(i,j,k-1));
    DataType km = (iter == accepted.end() ? inf : (*iter).second);

    // Fetch (i,j,k+1)
    iter = accepted.find(Coordinate(i,j,k+1));
    DataType kp = (iter == accepted.end() ? inf : (*iter).second);

    // Get the minimum along each dimension
    DataType f1 = std::min(im, ip);
    DataType f2 = std::min(jm, jp);
    DataType f3 = std::min(km, kp);

    // Sort the values in ascending order
    if (f1 > f2)  std::swap(f1,f2);
    if (f2 > f3)  std::swap(f2,f3);
    if (f1 > f2)  std::swap(f1,f2);
    
    mxAssert(f1 != inf, "Invalid tentative point - no accepted neighbors");
    
    // Case 1: We have accepted values only along one dimension
    if (f2 == inf)  return f1 + 1;
    
    // Case 2: We have accepted values along two dimensions
    if (f3 == inf) {
        DataType b = f1 + f2;
        DataType root = b*b - 2*(f1*f1 + f2*f2 - 1);
        if (root < 0)  root = 0;
        return (b + std::sqrt(root))*0.5;
    }

    // Case 3: We have accepted values along all three dimensions
    DataType b = f1 + f2 + f3;
    DataType root = b*b - 3*(f1*f1 + f2*f2 + f3*f3 - 1);
    if (root < 0)  root = 0;
    return (b + std::sqrt(root)) / 3.0;
}


template <typename DataType, typename IndexType>
void updateTentative(int i, int j, int k,
                     std::map<Coordinate, Node<DataType> *> & tentative,
                     std::map<Coordinate, DataType> & accepted,
                     Heap<DataType> & heap,
                     Array3D<DataType, IndexType> & A,
                     const DataType & bandwidth,
                     bool inside)
{
    // Check bounds
    if (i < 0 || i > A.getRows()-1 ||
        j < 0 || j > A.getCols()-1 ||
        k < 0 || k > A.getSlices()-1)
        return;
    
    Coordinate coord(i,j,k,inside);
    
    // Exit if already accepted
    if (accepted.find(coord) != accepted.end())  return;

    // Find the tentative value and see if it's been deleted
    // (it was outside the narrow band)
    typename std::map<Coordinate, Node<DataType> *>::iterator iter = tentative.find(coord);
    if (iter != tentative.end() && (*iter).second == NULL) return;

    // Otherwise, compute tentative value
    DataType value = computeTentativeValue(i,j,k, accepted);
    //mexPrintf("Tentative value for point (%i,%i,%i) computed: %f (bandwidth: %f)\n", i,j,k, value, bandwidth);
    
    // The point is new and within the narrowband
    if (iter == tentative.end() && value <= bandwidth) {
        Node<DataType> * node = new Node<DataType>(coord);
        node->cost = value;
        
        heap.push(node);
        tentative[coord] = node;
        //mexPrintf("Adding (%i,%i,%i) as new tentative point with value %f\n", i,j,k, value);
    }
    // The point is not new and within the narrowband
    else if (value <= bandwidth) {
        Node<DataType> * node = (*iter).second;
        node->cost = value;
        
        heap.update(node);
        //mexPrintf("Updating tentative point (%i,%i,%i) with value %f\n", i,j,k, value);
    }
    // The point is outside the narrowband and should be deleted
    // if possible
    else if (iter != tentative.end()) {
        //mexPrintf("Outside narrowband at (%i,%i,%i)\n", i,j,k);
        Node<DataType> * node = (*iter).second;
        heap.remove(node);
        delete node;
        (*iter).second = NULL;
    }
    //mexPrintf(" with value %f\n", value);
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
                             const DataType & jp,
                             const DataType & km,
                             const DataType & kp)
{
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
    if (km <= 0) {
        b = center / (center - km);
        if (val > b) val = b;
    }
    if (kp <= 0) {
        b = center / (center - kp);
        if (val > b) val = b;
    }
    
    return val;    
}


template <typename DataType, typename IndexType>
void compute(mxArray * plhs[], const mxArray * prhs[])
{
    typedef Array3D<DataType, IndexType> MyArray3D;
    typedef IndexArray<IndexType> MyIndexArray;
    typedef Node<DataType> MyNode;

    MyArray3D A(prhs[0]);
    MyIndexArray indices(prhs[1]);
    
    DataType bandwidth = *(static_cast<DataType *>(mxGetData(prhs[2])));
    if (bandwidth == 0)  mexPrintf("Warning: narrowband width is 0\n");

    // Create map of accepted grid points
    std::map<Coordinate, DataType> accepted;
    
    // Create map of tentative grid points
    std::map<Coordinate, MyNode *> tentative;

    // Heap used to store and fetch smallest tentative grid points
    Heap<DataType> heap;

    // Iterate the narrow band and find all points adjacent to
    // the zero-crossing
    typename MyArray3D::IndexedIterator iter = A.begin(indices);
    typename MyArray3D::IndexedIterator iend = A.end(indices);
    while (iter != iend) {
        int i = iter.getI();
        int j = iter.getJ();
        int k = iter.getK();
        
        DataType im = iter.get(-1, 0, 0);
        DataType ip = iter.get( 1, 0, 0);
        DataType jm = iter.get( 0,-1, 0);
        DataType jp = iter.get( 0, 1, 0);
        DataType km = iter.get( 0, 0,-1);
        DataType kp = iter.get( 0, 0, 1);
        DataType b, val;

        // If we have a zero-crossing, add the point to the accepted values
        // Flip the sign if we are inside
        if (*iter > 0 && (im <= 0 || ip <= 0 || jm <= 0 || jp <= 0 || km <= 0 || kp <= 0)) {
            accepted[Coordinate(i,j,k)] = interpolateDistance(*iter, im, ip, jm, jp, km, kp);
            //mexPrintf("Adding outside (%i,%i,%i) to accepted values (%f)\n", i,j,k, interpolateDistance(*iter, im, ip, jm, jp, km, kp));
        }
        else if (*iter <= 0 && (im > 0 || ip > 0 || jm > 0 || jp > 0 || km > 0 || kp > 0)) {
            accepted[Coordinate(i,j,k,true)] = interpolateDistance(-*iter, -im, -ip, -jm, -jp, -km, -kp);
            //mexPrintf("Adding inside (%i,%i,%i) to accepted values (%f)\n", i,j,k, -interpolateDistance(-*iter, -im, -ip, -jm, -jp, -km, -kp));
        }
        
        iter++;
    }
    
    // Iterate the accepted points and compute tentative distances
    // for all neighbors
    typename std::map<Coordinate, DataType>::const_iterator iterAccepted = accepted.begin();
    typename std::map<Coordinate, DataType>::const_iterator iendAccepted = accepted.end();
    while (iterAccepted != iendAccepted) {
        
        int i = (*iterAccepted).first.i;
        int j = (*iterAccepted).first.j;
        int k = (*iterAccepted).first.k;
        //mexPrintf("Visiting (%i,%i,%i)\n", i,j,k);

        updateTentative(i-1,j  ,k  , tentative, accepted, heap, A, bandwidth, (*iterAccepted).first.inside);
        updateTentative(i+1,j  ,k  , tentative, accepted, heap, A, bandwidth, (*iterAccepted).first.inside);
        updateTentative(i  ,j-1,k  , tentative, accepted, heap, A, bandwidth, (*iterAccepted).first.inside);
        updateTentative(i  ,j+1,k  , tentative, accepted, heap, A, bandwidth, (*iterAccepted).first.inside);
        updateTentative(i  ,j  ,k-1, tentative, accepted, heap, A, bandwidth, (*iterAccepted).first.inside);
        updateTentative(i  ,j  ,k+1, tentative, accepted, heap, A, bandwidth, (*iterAccepted).first.inside);
     
        iterAccepted++;
    }
    
    
    // Accept the smallest tentative value, update its neighbors
    // and continue until we're done
    while (!heap.isEmpty()) {
        MyNode * node = static_cast<MyNode *>(heap.pop());
        int i = node->coordinate.i;
        int j = node->coordinate.j;
        int k = node->coordinate.k;
        //mexPrintf("Accepting (%i,%i) with value %f\n", i, j, node->cost);

        accepted[node->coordinate] = node->cost;
        
        updateTentative(i-1,j  ,k  , tentative, accepted, heap, A, bandwidth, node->coordinate.inside);
        updateTentative(i+1,j  ,k  , tentative, accepted, heap, A, bandwidth, node->coordinate.inside);
        updateTentative(i  ,j-1,k  , tentative, accepted, heap, A, bandwidth, node->coordinate.inside);
        updateTentative(i  ,j+1,k  , tentative, accepted, heap, A, bandwidth, node->coordinate.inside);
        updateTentative(i  ,j  ,k-1, tentative, accepted, heap, A, bandwidth, node->coordinate.inside);
        updateTentative(i  ,j  ,k+1, tentative, accepted, heap, A, bandwidth, node->coordinate.inside);

        delete node;
    }
    
    
    // Create arrays for output
    plhs[0] = mxCreateNumericMatrix(1, accepted.size(), mxGetClassID(prhs[0]), mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, accepted.size(), mxGetClassID(prhs[1]), mxREAL);

    DataType * phi = static_cast<DataType *>(mxGetData(plhs[0]));
    IndexType * band = static_cast<IndexType *>(mxGetData(plhs[1]));
    
    // Loop through all accepted narrowband points
    iterAccepted = accepted.begin();
    IndexType ind = 0;
    while (iterAccepted != iendAccepted) {
        
        int i = (*iterAccepted).first.i;
        int j = (*iterAccepted).first.j;
        int k = (*iterAccepted).first.k;

        // Flip sign if we're on the inside
        if ((*iterAccepted).first.inside)
            phi[ind] = -(*iterAccepted).second;
        else
            phi[ind] = (*iterAccepted).second;
        
        band[ind] = A.sub2ind(i,j,k) + 1;  // Matlab uses 1-based indexing
     
        ind++;
        iterAccepted++;
    }
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

