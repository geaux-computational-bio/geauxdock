#ifndef _YEAH_OPENCL_INFO_H_
#define _YEAH_OPENCL_INFO_H_


#include <yeah/opencl/cl.h>


#ifdef __cplusplus
extern "C" {
#endif



typedef struct {
    const cl_device_info info;
    const char * str;
} Device_info;


typedef struct {
    const cl_platform_info info;
    const char * str;
} Platform_info;


typedef struct {
    const cl_program_info info;
    const char * str;
} Program_info;




void
opencl_device_type_to_string (const cl_device_type type, char* desc);

cl_device_type
opencl_device_type (const cl_device_id device);




void
list_devices (cl_platform_id platform, const int is_print);

void
list_platforms (const int is_print);



void
print_device_info (cl_device_id device);

void
print_platform_info (cl_platform_id platform);




void
print_program_info (cl_program program);

void
print_program_info_binary (cl_program program);



#ifdef __cplusplus
}
#endif


#endif
