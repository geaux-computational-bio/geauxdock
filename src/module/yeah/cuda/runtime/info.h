#ifndef _YEAH_CUDA_RUNTIME_INFO_H_
#define _YEAH_CUDA_RUNTIME_INFO_H_


#include <cuda_runtime.h>
#include <yeah/cuda/runtime/wrapper.h>


// /cuda/7.5/include/cuda_occupancy.h

struct cudaExeConfig {
    int blockDim;
    int gridDim;

    int gridDim_x;
    int gridDim_y;
    int gridDim_z;
    int blockDim_x;
    int blockDim_y;
    int blockDim_z;
};





// reference:
// nvidia/cuda7example/helper_cuda.h _ConvertSMVer2Cores_2

int Fp32UnitsPerMp (const cudaDeviceProp * const prop)
{
    const int major = prop->major;
    const int minor = prop->minor;


    int fp32units_per_mp = 0;

    typedef struct {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] =
    {
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
        { 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
        { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
        { 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
        {   -1, -1 }
    };

    int i = 0;
    while (nGpuArchCoresPerSM[i].SM != -1) {
        if (nGpuArchCoresPerSM[i].SM == ((major << 4) + minor)) {
            fp32units_per_mp = nGpuArchCoresPerSM[i].Cores;
            break;
        }
        i++;
    }

    // an alternative implementation
    /*
       const int fp32units_per_mp =
       major == 1 ? 8 :
       major == 2 ? (minor == 0 ? 32 : 48) :
       major == 3 ? 192 :
       major == 5 ? 128 :
       0;
     */


    if (fp32units_per_mp == 0)
        printf ("Fp32UnitsPerMp for SM %d.%d is undefined.\n", major, minor);
    return fp32units_per_mp;
}


int Fp64UnitsPerMp (const cudaDeviceProp * const prop)
{
    const int major = prop->major;
    const int minor = prop->minor;
    const int is_geforce = strncmp("GeForce", prop->name, 7) == 0;

    const int fp64units_per_mp =
        major == 1 ? 1 :
        major == 2 ? ( minor == 0 ? 16 : 4 ) :
        major == 3 ? ( minor < 3 || is_geforce ? 8 : 64 ) :
        major == 5 ? 4 :  // all geforce
        0;

    if (fp64units_per_mp == 0)
        printf ("Fp64UnitsPerMp for SM %d.%d is undefined.\n", major, minor);
    return fp64units_per_mp;
};




void
PrintCudaDeviceProp (const cudaDeviceProp * const prop, const int d)
{

    if (prop->major == 9999 && prop->minor == 9999) {
        printf ("device %d is invalid\n", d);
    }

    const char *NOorYES[] = {"NO", "YES"}; // 0-NO, 1-YES

    const int fp32units_per_mp = Fp32UnitsPerMp (prop);
    const int fp64units_per_mp = Fp64UnitsPerMp (prop);
    const double chip_bw_Bps = 2 * prop->memoryClockRate * 1000.0 * (prop->memoryBusWidth >> 3);
    const double chip_fp32_flops = 1000.0 * fp32units_per_mp * prop->clockRate * prop->multiProcessorCount;
    const double chip_fp64_flops = 1000.0 * fp64units_per_mp * prop->clockRate * prop->multiProcessorCount;

    printf ("GPU %d: %s (CC %d.%d)\n", d, prop->name, prop->major, prop->minor);
    //printf ("GPU %d:\n", d);
    //printf ("%s (CC %d.%d)\n", prop->name, prop->major, prop->minor);
    printf ("---------------------------------------\n");
    // core
    printf ("core clock rate \t\t%.2f GHz\n", prop->clockRate / 1e6);
    printf ("MP/GPU \t\t\t\t%d\n", prop->multiProcessorCount);
    printf ("FP32/MP \t\t\t%d\n", fp32units_per_mp);
    printf ("FP64/MP \t\t\t%d\n", fp64units_per_mp);
    printf ("FP32 FLOPS/GPU \t\t\t%.2f GFLOPS/s\n", chip_fp32_flops * 1e-9);
    printf ("FP64 FLOPS/GPU \t\t\t%.2f GFLOPS/s\n", chip_fp64_flops * 1e-9);
    putchar ('\n');
    // memory
    printf ("memory clock rate \t\t%.2f GHz\n", (float) prop->memoryClockRate / 1e6);
    printf ("memory bus width \t\t%d b\n", prop->memoryBusWidth);
    printf ("memory bandwidth \t\t%.2f GB/s\n", chip_bw_Bps * 1e-9);
    printf ("size reg/block \t\t\t%.3f KiB\n", (float) prop->regsPerBlock / 1024.0f );
    printf ("size scratchpad/block \t\t%.3f KiB\n", (float) prop->sharedMemPerBlock / 1024.0f);
    printf ("size L2 cache/GPU \t\t%.2f KiB\n", (float) prop->l2CacheSize / 1024.0f);
    printf ("size constant/GPU \t\t%.2f KiB\n", (float) prop->totalConstMem / 1024.0f);
    printf ("size global/GPU \t\t%.2f MiB\n", (float) prop->totalGlobalMem / (1024.0f * 1024.0f));
    printf ("\n");
    // core / memory
    printf ("FP32 FLOPS/BW ratio \t\t%.2f\n", 4 * chip_fp32_flops / chip_bw_Bps);
    printf ("FP64 FLOPS/BW ratio \t\t%.2f\n", 8 * chip_fp64_flops / chip_bw_Bps);
    putchar ('\n');
    // feature
    printf ("unifiedAddressing: \t\t%s\n", NOorYES[prop->unifiedAddressing]);
    printf ("canMapHostMemory: \t\t%s\n", NOorYES[prop->canMapHostMemory]); // support cudaHostAlloc(), cudaHostGetDevicePointer()
    printf ("concurrentKernels: \t\t%s\n", NOorYES[prop->concurrentKernels]); // simultaneously executing multiple kernels
    printf ("deviceOverlap: \t\t\t%s\n", NOorYES[prop->deviceOverlap]); // overlaping kernel and memcpy
    printf ("asyncEngineCount: \t\t%d\n", prop->asyncEngineCount); // overlaping kernel and memcpy (single/bi-direction)
    printf ("ECCEnabled:\t\t\t%s\n", NOorYES [prop->ECCEnabled]);
    putchar ('\n');
    // CUDA restrictions
    printf ("maxWarps: \t\t\t%d\n", prop->warpSize);
    printf ("maxThreadsPerBlock: \t\t%d\n", prop->maxThreadsPerBlock);
    printf ("maxThreadsPerMultiProcessor: \t%d\n", prop->maxThreadsPerMultiProcessor);
    printf ("kernelExecTimeoutEnabled: \t%s\n", NOorYES[prop->kernelExecTimeoutEnabled]);
    putchar ('\n');
}



void GetPrintCudaDeviceProp (const int d)
{
    cudaDeviceProp prop;
    CUDA_ERR (cudaGetDeviceProperties (&prop, d));
    PrintCudaDeviceProp (&prop, d);
}

void GetPrintCurrentCudaDeviceProp ()
{
    int d;
    cudaDeviceProp prop;
    CUDA_ERR (cudaGetDevice (&d));
    GetPrintCudaDeviceProp (d);
}






void PrintCudaFuncArributes (const cudaFuncAttributes * fa, const char * s)
{
    printf ("static cuda kernel information for \"%s\"\n", s);
    printf ("ptxVersion: \t\t\t%d\n", fa->ptxVersion);
    printf ("binaryVersion: \t\t\t%d\n", fa->binaryVersion);
    printf ("reg usage/T: \t\t\t%d\n", fa->numRegs);
    printf ("lmem usage/T: \t\t\t%.3f KiB\n", (float) fa->localSizeBytes / 1024.0f);
    printf ("smem usage/TB: \t\t\t%.3f KiB\n", (float) fa->sharedSizeBytes / 1024.0f);
    printf ("cmem usage/grid: \t\t%.3f KiB\n", (float) fa->constSizeBytes / 1024.0f);
    printf ("max T per TB: \t\t\t%d\n", fa->maxThreadsPerBlock);
    putchar ('\n');
}




// eg:
// GetPrintCudaFuncArributes ((void (*)) MyKernel, "MyKernel");

void GetPrintCudaFuncArributes (void (*func), const char * s)
{
    cudaFuncAttributes fa;
    CUDA_ERR (cudaFuncGetAttributes (&fa, func));
    PrintCudaFuncArributes (&fa, s);
}




// eg:
// GetPrint_potencial_occupancy ((void (*)) MyKernel, "MyKernel", 32, 512, 0);

void GetPrint_potencial_occupancy (void (*func), const char * s,
    const int bperg,
    const int tperb,
    const int dynamic_smem_sz)
{
    int d;
    cudaDeviceProp prop;
    CUDA_ERR (cudaGetDevice (&d));
    CUDA_ERR (cudaGetDeviceProperties (&prop, d));


    // report the potential maxium occupancy, given a execution configuration
    int threads_per_block = tperb;
    int active_blocks_per_mp;
    CUDA_ERR (cudaOccupancyMaxActiveBlocksPerMultiprocessor
        (&active_blocks_per_mp, func, threads_per_block, dynamic_smem_sz));
    printf ("warp occupancy for exec config <<<..., %d>>:\n", threads_per_block);
    printf ("active_blocks_per_mp: \t\t%d\n", active_blocks_per_mp);

    int active_warps_per_mp = threads_per_block * active_blocks_per_mp / prop.warpSize;
    int max_warps_per_mp = prop.maxThreadsPerMultiProcessor / prop.warpSize;
    float warp_occupancy = (float) active_warps_per_mp / (float) max_warps_per_mp;
    printf ("warp_occupancy: \t\t%d/%d, %.3f\n", active_warps_per_mp, max_warps_per_mp, warp_occupancy);
    printf ("\n");
}


// g++ found the fucntion prototype of "cudaOccupancyMaxActiveBlocksPerMultiprocessor"
// g++ can't find the fucntion prototype of "cudaOccupancyMaxPotentialBlockSize"
#if defined(__CUDACC__) // compiled by NVCC

// eg:
// GetPrint_execconfig_max_occupancy ((void (*)) MyKernel, "MyKernel", 0);

// suggest a execution configuration that maximize the warps occupancy
void GetPrint_execconfig_max_occupancy (void (*func), const char * s, const int dynamic_smem_sz)
{
    int minimal_blocks_per_grid, threads_per_block;
    const int blocks_per_grid_upper_bound = 0;
    CUDA_ERR (cudaOccupancyMaxPotentialBlockSize
        (&minimal_blocks_per_grid, &threads_per_block, func, dynamic_smem_sz, blocks_per_grid_upper_bound));
    printf ("suggested exec config:\t\t<<<%d(+), %d>>>\n", minimal_blocks_per_grid, threads_per_block);
    printf ("\n");

    GetPrint_potencial_occupancy (func, s, minimal_blocks_per_grid, threads_per_block, 0);
}



void GetPrintCudaFuncArributes2 (void (*func), const char * s,
    const int bperg,
    const int tperb,
    const int dynamic_smem_sz)
{
    GetPrint_execconfig_max_occupancy (func, s, 0);
    //GetPrint_potencial_occupancy (func, s, 32, 512, 0);
    //GetPrint_potencial_occupancy (func, s, 32, 768, 0);
    //GetPrint_potencial_occupancy (func, s, 64, 512, 0);
    GetPrint_potencial_occupancy (func, s, bperg, tperb, dynamic_smem_sz);

    //GetPrintCurrentCudaDeviceProp ();
}

#endif


#endif
