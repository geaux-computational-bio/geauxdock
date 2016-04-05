/*
#include <cstdlib>
#include <cstdio>

#include "geauxdock.h"
#include "gpu.cuh"

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
 */




__forceinline__ __device__
float static
MyRand_d (curandState * curandstate_d)
{
    const int gidx = blockDim.x * blockIdx.x + threadIdx.x;
    curandState myseed = curandstate_d[gidx];
    float randdd = curand_uniform (&myseed);
    curandstate_d[gidx] = myseed;

    return randdd;
}




/*
minimizeing:
ceil ((float) data_x / x) * ceil ((float) data_y / y)
ST:
x = 32, 64, 128, 256, ... blockDim.x
y = blockDim.x / x

 */

__device__
void
CalcThreadBlockShape (const int data_x, const int data_y, int & bdx, int & bdy)
{
    int var[11]; // at most 10 candidates, log2 (1024) = 10
    int imax = (int) log2 ((float) blockDim.x);

    for (int i = WARPshift; i <= imax; ++i) { // 32, 64, 128, 256, ... blockDim.x
        const int x = 1 << i;
        const int y = blockDim.x / x;
        var[i] = ceil ((float) data_x / x) * ceil ((float) data_y / y);
    }

    int iselect = WARPshift;
    for (int i = WARPshift; i <= imax; ++i) {
        if (var[iselect] >= var[i]) { // find the smallest (tend to adapt a larger i)
            iselect = i;
            //printf ("%03dx%03d: %3d     ", 1 << i, blockDim.x / (1 << i), var[i]);
        }
    }
    //printf ("\n");

    bdx = 1 << iselect;
    bdy = blockDim.x / bdx;

    //if (blockIdx.x == 0 && threadIdx.x == 0)
    //printf ("data: %3d x %3d     TB: %03dx%03d\n", data_x, data_y, bdx, bdy);
}

