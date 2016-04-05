
#ifndef _YEAH_C_UTIL_H_
#define _YEAH_C_UTIL_H_

#include <stdio.h>
#include <stdlib.h>


// usage:
// NULLPTR (a = malloc (xxx));
#define NULLPTR(prt) __GetNullPointerError1 (prt, __FILE__, __LINE__)


inline void
__GetNullPointerError1 (void * prt, char *file, int line)
{
    if (prt == NULL) {
        printf ("%s: line %d: null ptr failure\n", file, line);
        exit (EXIT_FAILURE);
    }
}







#if defined(__CUDACC__) // NVCC
#define ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
#define ALIGN(n) __attribute__((aligned(n)))
#elif defined(__INTEL_COMPILER) // INTEL
#define ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
#define ALIGN(n) __declspec(align(n))
#else
  #error "Missing the definition of "ALIGN macro" for this compiler."
#endif




#define RESTRICT __restrict__



#endif
