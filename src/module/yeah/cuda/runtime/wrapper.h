#ifndef _YEAH_CUDA_RUNTIME_WRAPPER_H_
#define _YEAH_CUDA_RUNTIME_WRAPPER_H_


#include <cuda_runtime.h>
#include <stdlib.h>







// CUDA API Error-Checking Wrappers
// This will output the proper CUDA error strings in the event that a CUDA host call returns an error

#ifdef USE_SIMPLE_WRAPPER

// reference:
// nvidia/cuda4example/cutil_inline_runtime.h: cudaSafeCall
#define CUDA_ERR(err) __GetError1 (err, __FILE__, __LINE__)

inline void
__GetError1 (const cudaError_t err, const char * const file, const int line)
{
  if (err != cudaSuccess) {
    printf ("%s(%d) : CUDA error : (%d) : %s.\n",
            file, line, (int) err, cudaGetErrorString (err));
    exit (EXIT_FAILURE);
  }
}

#else

// reference:
// nvidia/cuda7example/helper_cuda.h: checkCudaErrors
#include <nvidia/cuda7example/helper_cuda.h>
#ifdef __DRIVER_TYPES_H__
#define CUDA_ERR(err) check ( (err), #err, __FILE__, __LINE__ )
#endif

#endif







// This will output the proper error string when calling cudaGetLastError
//
// reference:
// nvidia/cuda4example/cutil_inline_runtime.h: cutilCheckMsg
// nvidia/cuda7example/helper_cuda.h: getLastCudaError


#define CUDA_LAST_ERR() __GetLastError1 (__FILE__, __LINE__)

inline void
__GetLastError1 (const char * const file, const int line)
{
  cudaError_t err = cudaGetLastError ();
  if (err != cudaSuccess) {
    printf ("%s(%d) : CUDA error : (%d) : %s.\n",
            file, line, (int) err, cudaGetErrorString (err));
    exit (EXIT_FAILURE);
  }
}






#define CUDAMALLOCHOST(ptr, sz, type) \
  CUDA_ERR (cudaMallocHost ((void **) &ptr, sz))

#define CUDAMALLOC(ptr, sz, type) \
  CUDA_ERR (cudaMalloc ((void **) &ptr, sz))

#define CUDAMEMCPY_SYMBOL(dst, src, type) \
  CUDA_ERR (cudaMemcpyToSymbol (dst, src, sizeof (type), 0, cudaMemcpyHostToDevice))

#define CUDAKERNEL_ASYNC(func, dim_grid, dim_block, ...) \
  func <<< dim_grid, dim_block >>> (__VA_ARGS__); \
  CUDA_LAST_ERR ()

#define CUDAKERNEL_STREAM_ASYNC(func, dim_grid, dim_block, n, stream, ...) \
  func <<< dim_grid, dim_block, n, stream >>> (__VA_ARGS__); \
  CUDA_LAST_ERR ()

#endif
