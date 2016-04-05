
#ifndef _YEAH_CUDA_NVML_INFO_H_
#define _YEAH_CUDA_NVML_INFO_H_


#include <stdio.h>
#include <stdlib.h>
#include <nvidia/nvml.h>
#include <yeah/cuda/nvml/wrapper.h>


#define NVML_STR_MAX_LENG 64

typedef struct
{
    nvmlDevice_t device;                            // important
    nvmlPciInfo_t pciinfo;                          // important
    nvmlEccErrorCounts_t eccerrorcounts;            // important
    nvmlUtilization_t utilization;                  // important
    nvmlMemory_t memory;                            // important
    nvmlBAR1Memory_t bar1memory;
    nvmlProcessInfo_t processinfo;
    nvmlBridgeChipType_t bridgechiptype;
    nvmlBridgeChipInfo_t bridgechipinfo;
    nvmlBridgeChipHierarchy_t bridgechiphierachy;
    nvmlSample_t sample;
    nvmlViolationTime_t violationtime;
    nvmlUnit_t unit;                                // important
    nvmlHwbcEntry_t hwbcentry;
    nvmlLedState_t ledstate;
    nvmlUnitInfo_t unitinfo;                        // important, name ...
    nvmlPSUInfo_t psuinfo;                          // important, may calcuate the power
    nvmlUnitFanInfo_t faninfo;
    nvmlUnitFanSpeeds_t fanspeeds;
    nvmlEventData_t eventdata;
    nvmlAccountingStats_t accountingstats;

    char name[NVML_STR_MAX_LENG];
} nvmlDeviceProp;





int
MatchCudaNvmlDevices (const cudaDeviceProp * const cuda_device_prop,
               const nvmlDeviceProp * const nvml_device_prop)
{
    int rv =
        (cuda_device_prop->pciDomainID == nvml_device_prop->pciinfo.domain) &&
        (cuda_device_prop->pciBusID == nvml_device_prop->pciinfo.bus) &&
        (cuda_device_prop->pciDeviceID == nvml_device_prop->pciinfo.device);

    return rv;
}






/*
 * https://github.com/al42and/cuda-smi
 *
 * For a number of reasons nVidia uses different device enumeration
 * in `nvidia-smi` monitoring utility and in their CUDA API,
 * making it extremely frustrating to choose GPU on multi-GPU machine.
 *
 * With the release of CUDA 7.0, it became possible to use `nvidia-smi` device order
 * in CUDA applications by setting environment variable `CUDA_DEVICE_ORDER=PCI_BUS_ID`
 * */

void
GetNvmlDeviceHandleByCudaID (nvmlDeviceProp *nvml_device_prop, const int id_cuda, int * id_nvml, int nd_nvml)
{
    cudaDeviceProp cuda_device_prop;
    CUDA_ERR (cudaGetDeviceProperties (&cuda_device_prop, id_cuda));

    int i_match = -1;
    for (int i = 0; i < nd_nvml; ++i) {
        NVML_ERR (nvmlDeviceGetHandleByIndex (i, &nvml_device_prop->device));
        NVML_ERR (nvmlDeviceGetPciInfo (nvml_device_prop->device, &nvml_device_prop->pciinfo));

        if (MatchCudaNvmlDevices (&cuda_device_prop, nvml_device_prop)) {
            i_match = i;
            *id_nvml = i;
            break;
        }
    }

    if (i_match == -1) {
        printf ("No NVML device matches CUDA ID\n");
        NVML_ERR (nvmlShutdown ());
        exit (EXIT_FAILURE);
    }

    printf ("CUDA Device %d = NVML device %d\n", id_cuda, *id_nvml);
}





void
GetNvmlDeviceProp (nvmlDeviceProp * const prop, const int id_nvml, const int id_cuda)
{
    const int one_MiB = 1 << 20;
    unsigned int temp;
    unsigned int enforcedpowerlimit;
    unsigned int powerlimit; // If the card's total power draw reaches this limit the power management algorithm kicks in
    unsigned int powerdefaultlimit;
    unsigned int powerconstaints_min, powerconstaints_max;
    unsigned int powerusage; // On Fermi and Kepler GPUs the reading is accurate to within +/- 5% of current power draw

    //NVML_ERR (nvmlUnitGetUnitInfo (prop->unit, &prop->unitinfo));
    NVML_ERR (nvmlDeviceGetName (prop->device, prop->name, NVML_STR_MAX_LENG));
    NVML_ERR (nvmlDeviceGetMemoryInfo (prop->device, &prop->memory));
    NVML_ERR (nvmlDeviceGetTemperature (prop->device, NVML_TEMPERATURE_GPU, &temp));
    NVML_ERR (nvmlDeviceGetEnforcedPowerLimit (prop->device, &enforcedpowerlimit));
    NVML_ERR (nvmlDeviceGetPowerManagementLimit (prop->device, &powerlimit));
    NVML_ERR (nvmlDeviceGetPowerManagementDefaultLimit (prop->device, &powerdefaultlimit));
    NVML_ERR (nvmlDeviceGetPowerManagementLimitConstraints (prop->device, &powerconstaints_min, &powerconstaints_max));
    NVML_ERR (nvmlDeviceGetPowerUsage (prop->device, &powerusage));


    printf ("%d: %s\n", id_nvml, prop->name);
    printf ("---------------------------------------\n");
    printf ("CUDA Device id \t\t\t\t\t%d\n", id_cuda);
    printf ("NVML Device id \t\t\t\t\t%d\n", id_nvml);
    printf ("PCIe \t\t\t\t\t\t%04x:%02x:%02x.0\n", prop->pciinfo.domain, prop->pciinfo.bus, prop->pciinfo.device);
    printf ("Memory usage (MiB) \t\t\t\t%zu / %zu\n", prop->memory.used / one_MiB, prop->memory.total / one_MiB);
    printf ("Temperature (C) \t\t\t\t%zu\n", temp);
    printf ("Enforced power limit (W) \t\t\t%.2f\n", (float) enforcedpowerlimit / 1000.0f);
    printf ("Power limit (W) \t\t\t\t%.2f\n", (float) powerlimit / 1000.0f);
    printf ("Power default limit (W) \t\t\t%.2f\n", (float) powerdefaultlimit / 1000.0f);
    printf ("Power management limit constraints (W) \t\t%.2f - %.2f\n",
            (float) powerconstaints_min / 1000.0f, (float) powerconstaints_max / 1000.0f);
    printf ("Power usage (W) \t\t\t\t%.2f\n", (float) powerusage / 1000.0f);
}



void
GetNvmlDeviceProp_2 (nvmlDeviceProp * const prop)
{
    const int one_MiB = 1 << 20;
    unsigned int temp;
    unsigned int powerusage; // On Fermi and Kepler GPUs the reading is accurate to within +/- 5% of current power draw

    NVML_ERR (nvmlDeviceGetMemoryInfo (prop->device, &prop->memory));
    NVML_ERR (nvmlDeviceGetTemperature (prop->device, NVML_TEMPERATURE_GPU, &temp));
    NVML_ERR (nvmlDeviceGetPowerUsage (prop->device, &powerusage));


    printf ("Memory usage (MiB) \t\t\t\t%zu / %zu\n", prop->memory.used / one_MiB, prop->memory.total / one_MiB);
    printf ("Temperature (C) \t\t\t\t%zu\n", temp);
    printf ("Power usage (W) \t\t\t\t%.2f\n", (float) powerusage / 1000.0f);
}






#endif
