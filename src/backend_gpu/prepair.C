#include <cstdio>

#include <geauxdock.h>

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <yeah/cuda/runtime/wrapper.h>


void
SetDevice ()
{
    for (int g = 0; g < NGPU; ++g) {
        cudaSetDevice (g);
        //cudaDeviceSetCacheConfig (cudaFuncCachePreferEqual);
        //cudaDeviceSetCacheConfig (cudaFuncCachePreferL1);
        cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
    }
}


static int
minimal_int (const int a, const int b)
{
    return a < b ? a : b;
}





void
DeviceAlloc (Complex ** c, Record **r, curandState **s)
{
    for (int g = 0; g < NGPU; g++) {
        cudaSetDevice (g);
        CUDA_ERR (cudaMalloc ((void **) &c[g], sizeof (Complex)));
        CUDA_ERR (cudaMalloc ((void **) &r[g], sizeof (Record) * MAX_REP));
        CUDA_ERR (cudaMalloc ((void **) &s[g], sizeof (curandState) * BD * GD));
    }
}



void
DeviceFree (Complex ** c, Record **r, curandState **s)
{
    for (int g = 0; g < NGPU; g++) {
        cudaSetDevice (g);
        CUDA_ERR (cudaFree (c[g]));
        CUDA_ERR (cudaFree (r[g]));
        CUDA_ERR (cudaFree (s[g]));
    }
}






// calculate the upper and lower bound
// for multi-GPU data decomposition
void
SetParaT (Complex *ch, ParaT **pt)
{
    const int n_rep = ch->size.n_rep;


    // communication optimization split
#if 0
    const int n_lig = ch->size.n_lig;
    const int n_prt = ch->size.n_prt;
    const int n_tmp = ch->size.n_tmp;

    const int n_rep_per_gpu_max =
        (int) ceilf ((float) (n_prt * n_tmp) / NGPU) * n_lig;
    const int n_active_gpu =
        (int) ceilf ((float) n_rep / n_rep_per_gpu_max);
    if (n_active_gpu < NGPU) {
        printf ("error: n_active_gpu < NGPU\n");
        exit (9999);
    }
#endif


    // maximizing balance split
#if 1
    const int n_rep_per_gpu_max =
        (int) ceilf ((float) n_rep / NGPU);
    const int n_active_gpu =
        (int) ceilf ((float) n_rep / n_rep_per_gpu_max);
    if (n_active_gpu < NGPU) {
        printf ("error: n_active_gpu < NGPU\n");
        exit (9999);
    }
#endif

    //printf ("NGPU = %d\n", NGPU);
    //printf ("n_rep_per_gpu_max = %d\n", n_rep_per_gpu_max);
    //printf ("n_avtive_gpu = %d\n", n_active_gpu);




    for (int g = 0; g < NGPU; ++g) {
        pt[g]->rep_begin = n_rep_per_gpu_max * g;
        pt[g]->rep_end = minimal_int (pt[g]->rep_begin + n_rep_per_gpu_max, n_rep) - 1;
        pt[g]->n_rep = pt[g]->rep_end - pt[g]->rep_begin + 1;
        pt[g]->record_sz = sizeof (Record) * pt[g]->n_rep;

        //printf ("%d replicas in GPU %d: %d - %d\n",
        //pt[g]->n_rep, g, pt[g]->rep_begin, pt[g]->rep_end);
    }
}



static void
CopyH2D_1device (Complex * ch, Complex * cd, ParaT * pt)
{
    ch->rep_begin = pt->rep_begin;
    ch->rep_end = pt->rep_end;
    CUDA_ERR (cudaMemcpyAsync (cd, ch, sizeof (Complex), cudaMemcpyHostToDevice));
}


void
CopyH2D (Complex * ch, Complex ** cd, ParaT ** pt)
{
    for (int g = 0; g < NGPU; g++) {
        cudaSetDevice (g);
        CopyH2D_1device (ch, cd[g], pt[g]);
    }
}


