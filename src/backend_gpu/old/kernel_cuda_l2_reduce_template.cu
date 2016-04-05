/*
#include <cstdlib>
#include <cstdio>

#include "geauxdock.h"
#include "gpu.cuh"

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
 */



// warp_id == bidx / warpSize == bidx >> WARPshift
// warp_land == bidx % warpSize == bidx & WARPmask
#ifndef WARPshift
#define WARPshift 5
#endif
#ifndef WARPmask
#define WARPmask 0b11111
#endif







template <typename T0, typename T1, typename T2, typename T3, typename T4>
__inline__ __device__
static void
BlockReduceSum_5_d (T0 &a0, T1 &a1, T2 &a2, T3 &a3, T4 &a4)
{

    const int bidx = threadIdx.x;

    __shared__ T0 a0s[WARPperB];
    __shared__ T1 a1s[WARPperB];
    __shared__ T2 a2s[WARPperB];
    __shared__ T3 a3s[WARPperB];
    __shared__ T4 a4s[WARPperB];

    // 1st level warp reduction
#pragma unroll
    for (int stride = 16; stride > 0; stride >>= 1) {
        a0 += __shfl_xor (a0, stride);
        a1 += __shfl_xor (a1, stride);
        a2 += __shfl_xor (a2, stride);
        a3 += __shfl_xor (a3, stride);
        a4 += __shfl_xor (a4, stride);
    }

    {
        const int warp_lane = bidx & WARPmask;
        const int warp_id = bidx >> WARPshift;
        if (warp_lane == 0) {
            a0s[warp_id] = a0;
            a1s[warp_id] = a1;
            a2s[warp_id] = a2;
            a3s[warp_id] = a3;
            a4s[warp_id] = a4;
        }
    }

    __syncthreads ();
    a0 = (bidx < WARPperB) ? a0s[bidx] : 0;
    a1 = (bidx < WARPperB) ? a1s[bidx] : 0;
    a2 = (bidx < WARPperB) ? a2s[bidx] : 0;
    a3 = (bidx < WARPperB) ? a3s[bidx] : 0;
    a4 = (bidx < WARPperB) ? a4s[bidx] : 0;

    // 2nd level warp reduction
#pragma unroll
    for (int stride = 16; stride > 0; stride >>= 1) {
        a0 += __shfl_xor (a0, stride);
        a1 += __shfl_xor (a1, stride);
        a2 += __shfl_xor (a2, stride);
        a3 += __shfl_xor (a3, stride);
        a4 += __shfl_xor (a4, stride);
    }
}









template <typename T0>
__inline__ __device__
static void
WarpReduceSum_1_d (T0 &a0)
{
#pragma unroll
    for (int stride = 16; stride > 0; stride >>= 1) {
        a0 += __shfl_xor (a0, stride);
    }
}






#if 0

template <typename T0, typename T1>
__inline__ __device__
static void
BlockReduceSum_2D_2_d (const int bdy, const int bdx, T0 &a0, T1 &a1)
{
    const int bidx = threadIdx.x;

    __shared__ T0 a0s[WARPperB];
    __shared__ T1 a1s[WARPperB];
    const int warp_lane = bidx & WARPmask;


    // warp reduction
#pragma unroll
    for (int stride = 16; stride > 0; stride >>= 1) {
        a0 += __shfl_xor (a0, stride);
        a1 += __shfl_xor (a1, stride);
    }

    {
        const int warp_id = bidx >> WARPshift;
        if (warp_lane == 0) {
            a0s[warp_id] = a0;
            a1s[warp_id] = a1;
        }
    }

    __syncthreads ();

    const int warp_x_per_b = bdx >> WARPshift;

    if (bidx < bdy) {
        a0 = 0;
        a1 = 0;
#pragma unroll
        for (int s = 0; s < warp_x_per_b; ++s) {
            a0 += a0s[warp_x_per_b * bidx + s];
            a1 += a1s[warp_x_per_b * bidx + s];
        }
    }


}

#endif







template <typename T0>
__inline__ __device__
static void
BlockReduceSum_2D_1_d (const int bdy, const int bdx, T0 &a0)
{
    const int bidx = threadIdx.x;

    __shared__ T0 a0s[WARPperB];
    const int warp_lane = bidx & WARPmask;


    // warp reduction
#pragma unroll
    for (int stride = 16; stride > 0; stride >>= 1) {
        a0 += __shfl_xor (a0, stride);
    }

    {
        const int warp_id = bidx >> WARPshift;
        if (warp_lane == 0) {
            a0s[warp_id] = a0;
        }
    }

    __syncthreads ();

    const int warp_x_per_b = bdx >> WARPshift;

    if (bidx < bdy) {
        a0 = 0;
#pragma unroll
        for (int s = 0; s < warp_x_per_b; ++s) {
            a0 += a0s[warp_x_per_b * bidx + s];
        }
    }

}
