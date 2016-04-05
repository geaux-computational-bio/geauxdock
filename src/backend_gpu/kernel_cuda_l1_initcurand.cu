// initialize curand status
// CURAND_Library.pdf, pp21

__global__ void
InitCurand_d (curandState* curandstate_d, const int seed)
{
    const int gidx = blockDim.x * blockIdx.x + threadIdx.x;

    curand_init (seed, gidx, 0, &curandstate_d[gidx]);

    // seed, subsequence, offset, gpuseed
    // skipahead(100000, &curandstate_d[gidx]);
}

