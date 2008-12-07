/*************************************************************************************************
 *
 * Copyright (c) 2008
 * Gunnar Läthén (gunnar.lathen@itn.liu.se)
 * Linköping University, Sweden.
 *
 *************************************************************************************************
 * Contributors:  
 *                  1) Ken Museth (ken.museth@itn.liu.se)
 *                  2) Gunnar Läthén (gunnar.lathen@itn.liu.se)
 *************************************************************************************************
 *
 * Triangulates the zero isosurface of a 3D scalar field using Marching cubes.
 *
 *************************************************************************************************/

#include "mex.h"
#include "Array3D.h"
#include "IndexArray.h"
#include "MC_Table.h"
#include <cmath>
#include <algorithm>
#include <map>
#include <vector>
#include <limits>

struct Edge
{
    Edge(unsigned int V0, unsigned int V1) : v0(V0), v1(V1) {
        if (v0 > v1) std::swap(v0, v1);
    }
    
    bool operator < (const Edge & edge) const {
        return  v0 < edge.v0 ||
                (v0 == edge.v0 && v1 < edge.v1);
    }
    
    unsigned int v0, v1;
};


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
Vector<DataType> computeNormal(unsigned int ind1, unsigned int ind2, DataType blend, const Array3D<DataType, IndexType> & A)
{
    DataType dx,dy,dz;
    int i,j,k;
    
    A.ind2sub(ind1, i,j,k);
    dx = (A(i+1,j,k) - A(i-1,j,k)) * 0.5;
    dy = (A(i,j+1,k) - A(i,j-1,k)) * 0.5;
    dz = (A(i,j,k+1) - A(i,j,k-1)) * 0.5;
    
    Vector<DataType> n1(dx,dy,dz);
    n1.normalize();
    
    A.ind2sub(ind2, i,j,k);
    dx = (A(i+1,j,k) - A(i-1,j,k)) * 0.5;
    dy = (A(i,j+1,k) - A(i,j-1,k)) * 0.5;
    dz = (A(i,j,k+1) - A(i,j,k-1)) * 0.5;
    
    Vector<DataType> n2(dx,dy,dz);
    n2.normalize();
    
    Vector<DataType> n = n1*(1-blend) + n2*blend;
    n.normalize();
        
    return n;
}


template <typename DataType, typename IndexType>
void compute(mxArray * plhs[], const mxArray * prhs[])
{
    typedef Array3D<DataType, IndexType> MyArray3D;
    typedef IndexArray<IndexType> MyIndexArray;

    MyArray3D A(prhs[0]);
    MyIndexArray indices(prhs[1]);

    std::map<Edge, unsigned int> edgeMap;
    std::vector<Vector<DataType> > vertices;
    std::vector<Vector<DataType> > normals;
    std::vector<Vector<unsigned int> > triangles;
    
    DataType voxelValues[8];
    unsigned int voxelIndices[8];

	typename MyArray3D::IndexedIterator iter = A.begin(indices);
    typename MyArray3D::IndexedIterator iend = A.end(indices);
    while (iter != iend) {

        int i = iter.getI();
        int j = iter.getJ();
        int k = iter.getK();
        
        voxelValues[0] = iter.get(0,0,0);
        voxelValues[1] = iter.get(0,1,0);
        voxelValues[2] = iter.get(1,1,0);
        voxelValues[3] = iter.get(1,0,0);
        voxelValues[4] = iter.get(0,0,1);
        voxelValues[5] = iter.get(0,1,1);
        voxelValues[6] = iter.get(1,1,1);
        voxelValues[7] = iter.get(1,0,1);
        
        // Determine the MC case
        short bitIndex = 0;
        for (int ind = 0; ind < 8; ind++)
            if (voxelValues[ind] > 0)
                bitIndex = bitIndex | MC::bit(ind);

        if (bitIndex > 0 && bitIndex < 255) {
            //mexPrintf("Visiting (%i,%i,%i): %f ", i, j, k, *iter);
            //mexPrintf("(Case %i)\n", bitIndex);
            voxelIndices[0] = A.sub2ind(i,  j,  k  );
            voxelIndices[1] = A.sub2ind(i,  j+1,k  );
            voxelIndices[2] = A.sub2ind(i+1,j+1,k  );
            voxelIndices[3] = A.sub2ind(i+1,j,  k  );
            voxelIndices[4] = A.sub2ind(i,  j,  k+1);
            voxelIndices[5] = A.sub2ind(i,  j+1,k+1);
            voxelIndices[6] = A.sub2ind(i+1,j+1,k+1);
            voxelIndices[7] = A.sub2ind(i+1,j,  k+1);

            for (int iTri = MC::first(bitIndex); iTri < MC::last(bitIndex); iTri++) {
                int * iEdge = MC::Vtx(iTri); // three edges crossed by triangle iTri;
                Vector<unsigned int> tri;// holds 3 pointer into mVtxList
                for (int ind = 0; ind < 3; ind++) {//loop over 3 vertices located on voxel edges
                    int vtx1 = MC::EdgeVtx1(iEdge[ind]);//first vertex on voxel edge
                    int vtx2 = MC::EdgeVtx2(iEdge[ind]);//second vertex on voxel edge

                    Edge edge(voxelIndices[vtx1], voxelIndices[vtx2]);
                    //mexPrintf("Looking up edge (%i,%i)\n", edge.v0, edge.v1);
                    std::map<Edge, unsigned int>::iterator iterEdge = edgeMap.find(edge);
                    if (iterEdge == edgeMap.end()) {
                        //mexPrintf("Edge not found - creating new\n");
                        DataType blend = voxelValues[vtx1] / (voxelValues[vtx1] - voxelValues[vtx2]);//0<= between <=1 
                        Vector<DataType> vertex(MC::X(vtx1)+i, MC::Y(vtx1)+j, MC::Z(vtx1)+k);
                        vertex[MC::EdgeAxis(iEdge[ind])] += blend;
                        Vector<DataType> normal = computeNormal(voxelIndices[vtx1], voxelIndices[vtx2], blend, A);

                        edgeMap[edge] = vertices.size();
                        tri[ind] = vertices.size();

                        vertices.push_back(vertex);
                        normals.push_back(normal);
                    }
                    else {
                        //mexPrintf("Edge found - fetching\n");
                        tri[ind] = (*iterEdge).second;
                    }

                }
                triangles.push_back(tri);
            }
        }

        iter++;
    }
    
    //mexPrintf("Isosurface extraction complete:\n");
    //mexPrintf("  Vertices/Normals: %i\n", vertices.size());
    //mexPrintf("  Triangles: %i\n", triangles.size());

    plhs[0] = mxCreateNumericMatrix(triangles.size(), 3, mxUINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(vertices.size(), 3, mxGetClassID(prhs[0]), mxREAL);
    plhs[2] = mxCreateNumericMatrix(normals.size(), 3, mxGetClassID(prhs[0]), mxREAL);

    unsigned int * tri = static_cast<unsigned int *>(mxGetData(plhs[0]));
    DataType * vert = static_cast<DataType *>(mxGetData(plhs[1]));
    DataType * norm = static_cast<DataType *>(mxGetData(plhs[2]));

    IndexType ind = 0;
    std::vector<Vector<unsigned int> >::const_iterator triIter = triangles.begin();
    std::vector<Vector<unsigned int> >::const_iterator triIend = triangles.end();
    while (triIter != triIend) {
        tri[ind] = (*triIter)[0]+1;
        tri[ind+triangles.size()] = (*triIter)[1]+1;
        tri[ind+triangles.size()*2] = (*triIter)[2]+1;
        ind++;
        triIter++;
    }

    ind = 0;
    typename std::vector<Vector<DataType> >::const_iterator vertIter = vertices.begin();
    typename std::vector<Vector<DataType> >::const_iterator vertIend = vertices.end();
    while (vertIter != vertIend) {
        vert[ind] = (*vertIter)[1]+1;
        vert[ind+vertices.size()] = (*vertIter)[0]+1;
        vert[ind+vertices.size()*2] = (*vertIter)[2]+1;
        ind++;
        vertIter++;
    }

    ind = 0;
    typename std::vector<Vector<DataType> >::const_iterator normIter = normals.begin();
    typename std::vector<Vector<DataType> >::const_iterator normIend = normals.end();
    while (normIter != normIend) {
        norm[ind] = (*normIter)[1];
        norm[ind+normals.size()] = (*normIter)[0];
        norm[ind+normals.size()*2] = (*normIter)[2];
        ind++;
        normIter++;
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check for proper number of arguments
    if (nrhs != 2)
        mexErrMsgTxt("Two input arguments required.");
    if (nlhs != 3)
        mexErrMsgTxt("Three output arguments required.");

    // The input must be noncomplex
    if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) ||
        !mxIsNumeric(prhs[0]) || !mxIsNumeric(prhs[1]))
        mexErrMsgTxt("The input must be noncomplex (and numeric).");

    // The second input must be an unsigned integer
    if (mxIsSingle(prhs[1]) || mxIsDouble(prhs[1]))
        mexErrMsgTxt("The second input must be integer.");

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

