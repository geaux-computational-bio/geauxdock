// from ppcg git 201506


#include <stdio.h>
#include <stdlib.h>
#include <yeah/opencl/cl.h>
#include <yeah/opencl/util.h>
#include <yeah/opencl/wrapper.h>
#include <yeah/opencl/info.h>



/* Find a GPU or a CPU associated with the first available platform.
 * If use_gpu is set, then this function first tries to look for a GPU
 * in the first available platform.
 * If this fails or if use_gpu is not set, then it tries to use the CPU.
 */
cl_device_id
opencl_create_device (int use_gpu)
{
    cl_platform_id platform;
    cl_device_id dev;
    int err;

    CL_ERR (clGetPlatformIDs (1, &platform, NULL));
    err = CL_DEVICE_NOT_FOUND;
    if (use_gpu)
        err = clGetDeviceIDs (platform, CL_DEVICE_TYPE_GPU, 1, &dev, NULL);
    if (err == CL_DEVICE_NOT_FOUND)
        CL_ERR (clGetDeviceIDs (platform, CL_DEVICE_TYPE_CPU, 1, &dev, NULL));
    return dev;
}



cl_device_id
opencl_create_device_2 (cl_device_type dev_type)
{
    cl_platform_id platform;
    cl_device_id dev;

    CL_ERR (clGetPlatformIDs (1, &platform, NULL));
    CL_ERR (clGetDeviceIDs (platform, dev_type, 1, &dev, NULL));
    return dev;
}





/*
 * Enumerate all platform & devices, and choose the first match
 * CL_DEVICE_TYPE_CPU
 * CL_DEVICE_TYPE_GPU
 * CL_DEVICE_TYPE_ACCELERATOR
 */ 

cl_device_id
opencl_choose_device (cl_device_type dev_type)
{
#define MAX_PLATFORMS 10
#define MAX_DEVICES 10

    cl_platform_id platforms[MAX_PLATFORMS];
    cl_uint platforms_n = 0;
    cl_device_id devices[MAX_DEVICES];
    cl_uint devices_n = 0;  // devices per platform


    CL_ERR (clGetPlatformIDs (MAX_PLATFORMS, platforms, &platforms_n));
    if (platforms_n < 1) {
        printf ("Error: no platforms found while looking for a platform.\n");
        exit (EXIT_FAILURE);
    }

    cl_uint p;
    for (p = 0; p < platforms_n; ++p) {
        print_platform_info (platforms[p]);
        CL_ERR (clGetDeviceIDs (platforms[p], CL_DEVICE_TYPE_ALL, MAX_DEVICES, devices, &devices_n));
        //CL_ERR (clGetDeviceIDs (platforms[p], dev_type, MAX_DEVICES, devices, &devices_n));

        cl_uint d;
        for (d = 0; d < devices_n; ++d) {
            print_device_info (devices[d]);
            if (dev_type == opencl_device_type (devices[d])) {
                printf ("Choose device %i\n", devices[d]);
                return devices[d];
            }
        }
    }

    printf ("No requested type of device found in the system\n");
    exit (EXIT_FAILURE);
}









/* Create an OpenCL program from a string and compile it.
 */
cl_program
opencl_build_program_from_string (cl_context ctx, cl_device_id dev,
    const char *program_source,
    size_t program_size,
    const char *opencl_options)
{
    int err;
    cl_program program;
    char *program_log;
    size_t log_size;

    program = clCreateProgramWithSource (ctx, 1, &program_source, &program_size, &err);
    if (err < 0) {
        fprintf (stderr, "Could not create the program\n");
        exit (EXIT_FAILURE);
    }
    err = clBuildProgram (program, 0, NULL, opencl_options, NULL, NULL);
    if (err < 0) {
        fprintf (stderr, "Could not build the program.\n");
        CL_ERR (clGetProgramBuildInfo (program, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size));
        program_log = (char *) malloc (log_size + 1);
        program_log[log_size] = '\0';
        CL_ERR (clGetProgramBuildInfo (program, dev, CL_PROGRAM_BUILD_LOG, log_size + 1, program_log, NULL));
        fprintf (stderr, "%s\n", program_log);
        free (program_log);
        exit (EXIT_FAILURE);
    }
    return program;
}



/* Create an OpenCL program from a source file and compile it.
 */
cl_program
opencl_build_program_from_file (cl_context ctx, cl_device_id dev,
    const char *filename,
    const char *opencl_options)
{
    cl_program program;
    FILE *program_file;
    char *program_source;
    size_t program_size, read;

    program_file = fopen (filename, "r");
    if (program_file == NULL) {
        fprintf (stderr, "Could not find the source file.\n");
        exit (EXIT_FAILURE);
    }
    fseek (program_file, 0, SEEK_END);
    program_size = ftell (program_file);
    rewind (program_file);
    program_source = (char *) malloc (program_size + 1);
    program_source[program_size] = '\0';
    read = fread (program_source, sizeof (char), program_size, program_file);
    if (read != program_size) {
        fprintf (stderr, "Error while reading the kernel.\n");
        exit (EXIT_FAILURE);
    }
    fclose (program_file);

    program = opencl_build_program_from_string (ctx, dev, program_source, program_size, opencl_options);
    free (program_source);

    return program;
}



