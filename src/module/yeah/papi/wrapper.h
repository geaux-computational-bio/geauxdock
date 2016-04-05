
#ifndef _YEAH_PAPI_WRAPPER_H_
#define _YEAH_PAPI_WRAPPER_H_

#include <memory.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <papi.h>



// can't use wrapper in the followsing APIs
//PAPI_ERR (PAPI_shutdown ()); // error: invalid use of void expression





#define PAPI_ERR(err) __GetPapiError1 (err, __FILE__, __LINE__)

#define PAPI_ERR_INIT(err) __GetPapiError2 (err, __FILE__, __LINE__)



#ifdef __cplusplus
extern "C" {
#endif

void
__GetPapiError1 (int err, const char * const file, int line)
{
    if (err < PAPI_OK) {
        printf ("%s\tFAILED\nLine # %d\n", file, line);
        if (err == PAPI_ESYS) {
            char buf[128];
            memset (buf, '\0', sizeof (buf));
            sprintf (buf, "System error:");
            perror (buf);
        }
        else if (err > 0) {
            printf ("Error calculating: \n");
        }
        else {
            printf ("Error: %s\n", PAPI_strerror (err));
        }

        exit (EXIT_FAILURE);
    }
}



void
__GetPapiError2 (int err, const char * const file, int line)
{
    if (err != PAPI_VER_CURRENT) {
        printf ("%s\tFAILED\nLine # %d\n", file, line);
        printf ("Error: %s\n", PAPI_strerror (err));
        exit (EXIT_FAILURE);
    }
}


#ifdef __cplusplus
}
#endif




#endif
