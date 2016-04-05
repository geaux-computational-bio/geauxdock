#ifndef _YEAH_OPENCL_WRAPPER_H_
#define _YEAH_OPENCL_WRAPPER_H_


#include <assert.h>
#include <yeah/opencl/cl.h>

//#ifdef __cplusplus
//extern "C" {
//#endif


/* Return the OpenCL error string for a given error number.
 */
const char * opencl_error_string (cl_int error);

const char * opencl_error_string_2 (cl_int error);




#define CL_ERR(err) __GetErrCL (err, __FILE__, __LINE__)


inline void
__GetErrCL (const cl_int err, const char * const file, const int line)
{
    if (err != CL_SUCCESS) {
        printf ("%s(%d) : OpenCL error : (%d) : %s.\n",
            file, line, (int) err,opencl_error_string_2 (err));
        fflush (stderr);
        assert (err == CL_SUCCESS);
    }
}


//#ifdef __cplusplus
//}
//#endif



#endif

