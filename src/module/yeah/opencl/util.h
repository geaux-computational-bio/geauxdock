// from ppcg git 201506

#ifndef _YEAH_OPENCL_UTIL_H_
#define _YEAH_OPENCL_UTIL_H_

#include <yeah/opencl/cl.h>


#ifdef __cplusplus
extern "C" {
#endif


/* Find a GPU or a CPU associated with the first available platform.
 * If use_gpu is set, then this function first tries to look for a GPU
 * in the first available platform.
 * If this fails or if use_gpu is not set, then it tries to use the CPU.
 */
cl_device_id opencl_create_device (int use_gpu);


cl_device_id opencl_create_device_2 (cl_device_type dev_type);


cl_device_id opencl_choose_device (cl_device_type dev_type);



/* Create an OpenCL program from a string and compile it.
 */
cl_program opencl_build_program_from_string (cl_context ctx,
                        cl_device_id dev,
                        const char *program_source,
                        size_t program_size,
                        const char *opencl_options);




/* Create an OpenCL program from a source file and compile it.
 */
cl_program opencl_build_program_from_file (cl_context ctx,
                        cl_device_id dev,
                        const char *filename,
                        const char *opencl_options);


#ifdef __cplusplus
}
#endif


#endif
