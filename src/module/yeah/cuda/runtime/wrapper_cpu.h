/*
 * implement some CUDA APIs in C stdlib
 * to help write universal code for both CPU/GPU
 */

#ifndef _YEAH_CUDA_RUNTIME_WRAPPER_CPU_H_
#define _YEAH_CUDA_RUNTIME_WRAPPER_CPU_H_



#define cudaFuncSetCacheConfig(...) \
    ;

#define cudaSetDevice(...) \
    ;

#define cudaGetLastError() \
    ;




#if 0
// emulate CUDA malloc/free via CPU malloc/free
#define CUDAMALLOC(ptr, sz, type) \
    ptr = (type) malloc (sz)

#define CUDAFREE(...) \
    free (__VA_ARGS__);

#define CUDAMEMCPY(dst, src, sz, direction) \
    memcpy (dst, src, sz);
#endif

#if 1

// emulate CUDA malloc/free via pointer assignment
#define CUDAMALLOC(ptr, sz, type) \
    ;

#define CUDAFREE(...) \
    ;

#define CUDAMEMCPY(dst, src, sz, direction) \
    dst = src;

#endif



#define CUDAMEMCPYTOSYMBOL(dst, src, type)\
    dst = *src;

#define CUDAKERNEL_SYNC(funcname, dim_grid, dim_block, ...) \
    funcname (__VA_ARGS__);


#endif



