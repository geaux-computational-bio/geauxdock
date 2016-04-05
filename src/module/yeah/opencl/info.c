#include <stdio.h>
#include <string.h>
#include <yeah/opencl/cl.h>
#include <yeah/opencl/info.h>
#include <yeah/opencl/wrapper.h>


#define STR_MAX_LEN 256


void
opencl_device_type_to_string (const cl_device_type type, char* desc)
{
    strcpy (desc, "");

    if (CL_DEVICE_TYPE_CPU == (type & CL_DEVICE_TYPE_CPU))
        strcat (desc,"CPU ");
    if (CL_DEVICE_TYPE_GPU == (type & CL_DEVICE_TYPE_GPU))
        strcat (desc,"GPU ");
    if (CL_DEVICE_TYPE_ACCELERATOR == (type & CL_DEVICE_TYPE_ACCELERATOR))
        strcat (desc,"ACCELERATOR ");
    if (CL_DEVICE_TYPE_DEFAULT == (type & CL_DEVICE_TYPE_DEFAULT))
        strcat (desc,"DEFAULT");
}



cl_device_type
opencl_device_type (const cl_device_id device)
{
    cl_device_type info;
    size_t sz = sizeof (cl_device_type);
    CL_ERR (clGetDeviceInfo (device, CL_DEVICE_TYPE, sz, &info, NULL));
    return info;
}









void
print_program_info (cl_program program)
{

    const Program_info program_info[] = {
        {CL_PROGRAM_SOURCE,       "CL_PROGRAM_SOURCE"}
    };

    int j;
    for (j = 0; j < sizeof (program_info) / sizeof (Program_info); j++) {
        size_t sz;
        CL_ERR (clGetProgramInfo (program, program_info[j].info, 0, NULL, &sz));
        char * info = (char*) malloc (sz);
        CL_ERR (clGetProgramInfo (program, program_info[j].info, sz, info, NULL));
        printf("%s:\n%s\n", program_info[j].str, info);
        free(info);
    }

    printf ("\n");
}




// ??????????????????????
// http://stackoverflow.com/questions/12868889/clgetprograminfo-cl-program-binary-sizes-incorrect-results

void
print_program_info_binary (cl_program program)
{

#if 1
    // GOOD, arrayfire
    size_t sz;
    CL_ERR (clGetProgramInfo (program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &sz, NULL));
    unsigned char *info = (unsigned char *) malloc (sz);
    CL_ERR (clGetProgramInfo (program, CL_PROGRAM_BINARIES, sizeof (unsigned char *), &info, NULL));
    printf("%s:\n%s\n", "Kenerl Target Code", info); 
#endif

#if 0
    // bad
    size_t sz;
    CL_ERR (clGetProgramInfo (program, CL_PROGRAM_BINARIES, 0, NULL, &sz));
    unsigned char *info = (unsigned char *) malloc (sz);
    CL_ERR (clGetProgramInfo (program, CL_PROGRAM_BINARIES, sizeof (unsigned char *), &info, NULL));
    printf("%s:\n%s\n", "Kenerl Target Code", info); 
#endif


#if 1
    FILE *fp = fopen("a.ptx", "wb");
    fwrite (info, sizeof(char), sz, fp);
    fclose (fp);
#endif
    free (info);

    printf ("\n");
}









void
print_device_info (cl_device_id device)
{

    const Device_info device_info_str[] = {
        {CL_DEVICE_NAME,             "CL_DEVICE_NAME"},
        {CL_DEVICE_VERSION,          "CL_DEVICE_VERSION"},
        {CL_DRIVER_VERSION,          "CL_DRIVER_VERSION"},
        {CL_DEVICE_OPENCL_C_VERSION, "CL_DEVICE_OPENCL_C_VERSION"}
    };

    const Device_info device_info_uint[] = {
        {CL_DEVICE_MAX_CLOCK_FREQUENCY,           "CL_DEVICE_MAX_CLOCK_FREQUENCY"},
        {CL_DEVICE_ADDRESS_BITS,                  "CL_DEVICE_ADDRESS_BITS"},
        {CL_DEVICE_MAX_COMPUTE_UNITS,             "CL_DEVICE_MAX_COMPUTE_UNITS"},
        {CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,   "CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR"},
        {CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,  "CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT"},
        {CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,    "CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT"},
        {CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,   "CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG"},
        {CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF,   "CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF"},
        {CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,  "CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT"},
        {CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE"}
    };


    int j;
    for (j = 0; j < sizeof (device_info_str) / sizeof (Device_info); ++j) {
        size_t sz;
        CL_ERR (clGetDeviceInfo (device, device_info_str[j].info, 0, NULL, &sz));
        char *info = (char *) malloc (sz);
        CL_ERR (clGetDeviceInfo (device, device_info_str[j].info, sz, info, NULL));
        printf ("          %-40s: %s\n", device_info_str[j].str, info);
        free (info);
    }


    for (j = 0; j < sizeof (device_info_uint) / sizeof (Device_info); ++j) {
        cl_uint info;
        size_t sz = sizeof (info);
        CL_ERR (clGetDeviceInfo (device, device_info_uint[j].info, sz, &info, NULL));
        printf ("          %-40s: %d\n", device_info_uint[j].str, info);
    }


    cl_device_type info;
    size_t sz = sizeof (cl_device_type);
    char str[STR_MAX_LEN];
    CL_ERR (clGetDeviceInfo (device, CL_DEVICE_TYPE, sz, &info, NULL));
    opencl_device_type_to_string (info, str);
    printf ("          %-40s: %s\n", "CL_DEVICE_TYPE", str);

    printf ("\n");

}




void
print_platform_info (cl_platform_id platform)
{
    const Platform_info platform_info[] = {
        {CL_PLATFORM_NAME,            "CL_PLATFORM_NAME"},
        {CL_PLATFORM_VENDOR,          "CL_PLATFORM_VENDOR"},
        {CL_PLATFORM_VERSION,         "CL_PLATFORM_VERSION"},
        {CL_PLATFORM_PROFILE,         "CL_PLATFORM_PROFILE"},
        {CL_PLATFORM_EXTENSIONS,      "CL_PLATFORM_EXTENSIONS"}
    };

    int j;
    for (j = 0; j < sizeof (platform_info) / sizeof (Platform_info); j++) {
        size_t sz;
        CL_ERR (clGetPlatformInfo (platform, platform_info[j].info, 0, NULL, &sz));
        char *info = (char *) malloc (sz);
        CL_ERR (clGetPlatformInfo (platform, platform_info[j].info, sz, info, NULL));
        printf ("      %-30s: %s\n", platform_info[j].str, info);
        free (info);
    }
    printf ("\n");

}









void
list_devices (cl_platform_id platform, const int is_print)
{
    cl_uint device_count;
    CL_ERR (clGetDeviceIDs (platform, CL_DEVICE_TYPE_ALL, 0, NULL, &device_count));
    cl_device_id *device = (cl_device_id *) malloc (sizeof (cl_device_id) * device_count);
    CL_ERR (clGetDeviceIDs (platform, CL_DEVICE_TYPE_ALL, device_count, device, NULL));

    int i;
    for (i = 0; i < device_count; i++) {
        printf ("      .%d %-11s\n", i, "Device");
        if (is_print)
            print_device_info (device[i]);
    }

    free (device);
}




void
list_platforms (const int is_print)
{
    cl_uint platform_count;
    CL_ERR (clGetPlatformIDs (1, NULL, &platform_count));
    cl_platform_id *platform = (cl_platform_id *) malloc (sizeof (cl_platform_id) * platform_count);
    CL_ERR (clGetPlatformIDs (platform_count, platform, NULL));

    int i;
    for (i = 0; i < platform_count; i++) {
        printf ("  .%d %-11s\n", i, "Platform");
        if (is_print)
            print_platform_info (platform[i]);
        list_devices (platform[i], is_print);
    }

    free (platform);
}


