
#ifndef _LS_COMMON
#define _LS_COMMON

#define SHARED_DATA s_Data
#define SHARED_DATA_DEF float SHARED_DATA[BLOCKDIM_Z + 2*PADDING][BLOCKDIM_Y + 2*PADDING][BLOCKDIM_X + 2*PADDING]

#define DATA SHARED_DATA[threadIdx.z][threadIdx.y][threadIdx.x]
#define NHOOD(i,j,k) SHARED_DATA[threadIdx.z+k][threadIdx.y+j][threadIdx.x+i]

#endif
