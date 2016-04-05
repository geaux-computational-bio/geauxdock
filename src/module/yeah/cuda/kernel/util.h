#ifndef _YEAH_CUDA_KERNEL_UTIL_H_
#define _YEAH_CUDA_KERNEL_UTIL_H_

// supported by CC > 3.5
// force the load to go through the read-only data cache
// generally
// significantly improve performance on Kepler GPU
// slightly improves performance on Maxwell GPU


//#define ENABLE_CUDA_LDG


#ifdef ENABLE_CUDA_LDG
#define CUDA_LDG_P(p) __ldg (p)
#define CUDA_LDG_D(d) __ldg (&d)
#else
#define CUDA_LDG_P(p) *p
#define CUDA_LDG_D(d) d
#endif


#endif
