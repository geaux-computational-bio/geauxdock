#ifndef _YEAH_CUDA_DRIVER_WRAPPER_H_
#define _YEAH_CUDA_DRIVER_WRAPPER_H_


#include <cuda.h>


// CUDA API Error-Checking Wrappers
// This will output the proper CUDA error strings in the event that a CUDA host call returns an error

// reference:
// nvidia/cuda7example/helper_cuda.h: checkCudaErrors
#include <nvidia/cuda7example/helper_cuda.h>
//#ifdef __DRIVER_TYPES_H__
#define CU_ERR(err) check ( (err), #err, __FILE__, __LINE__ )
//#endif


#endif
