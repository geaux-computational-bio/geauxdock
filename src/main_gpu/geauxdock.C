//#include <mpi.h>

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstdio>


#include <size.h>
#include <geauxdock.h>
#include <util.h>
#include <util_optimize.h>
#include <util_print.h>
#include <toggle.h>
#include <load.h>
//#include <sys/resource.h>

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>


#include "../backend_gpu/prepair.h"
#include "../backend_gpu/backend_gpu.h"
#include <yeah/c/timing.h>
#include <yeah/cpp/timer.hpp>
#include <yeah/cuda/runtime/wrapper.h>





int main (int argc, char **argv)
{
/*
  // increase the stack size
  const rlim_t kStackSize = 16 * 1024 * 1024;   // min stack size = 16 MB
  struct rlimit rl;
  int result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0) {
      if (rl.rlim_cur < kStackSize) {
          rl.rlim_cur = kStackSize;
          result = setrlimit(RLIMIT_STACK, &rl);
          if (result != 0)
              printf ("setrlimit returned result = %d\n", result);
      }
  }
*/



#include "../frontend/server1.C"

  yeah::Timer tc[8];

  tc[0].Start ();

  printf ("Initilizing GPU\n");

  // host data
  Complex *ch;
  Record *rh;
  CUDA_ERR (cudaMallocHost ((void **) &ch, sizeof (Complex)));
  CUDA_ERR (cudaMallocHost ((void **) &rh, sizeof (Record) * MAX_REP));
  //ch = (Complex *) malloc (sizeof (Complex));
  //rh = (Record *) malloc (sizeof (Record) * MAX_REP);

  printf ("sizeof record = %f MB\n", (float) sizeof (Record) * MAX_REP / 1024 / 1024);


  // device data and trandfer parameters
  Complex *cd[NGPU];
  Record *rd[NGPU];
  curandState *sd[NGPU];
  ParaT parat[NGPU];
  ParaT *pt[NGPU];
  for (int g = 0; g < NGPU; g++)
    pt[g] = &parat[g];


  SetDevice ();
  DeviceAlloc (cd, rd, sd);
  InitCurand (sd);

  tc[0].Stop ();

  tc[1].Start ();

  //int i = 2;
  for (int i = 0; i < ncomplex; ++i) {
    t[2].Start ();
    *ch = complex[i];
    printf ("complex %d\n", i);


    printf ("%s start docking\n", argv[0]);
    SetParaT (ch, pt);
    CopyH2D (ch, cd, pt);
    t[2].Stop ();

    t[3].Reset ();
    t[3].Start ();
    Dock (ch, rh, cd, rd, pt, sd);
    t[3].Stop ();
    printf ("copy H2D time\t\t\t%8.3f\n", t[2].Span ());
    printf ("Dock time\t\t\t%8.3f\n", t[15].Span ());
  }


  tc[1].Stop ();


  tc[3].Start ();

  DeviceFree (cd, rd, sd);
  CUDA_ERR (cudaFreeHost (rh));
  CUDA_ERR (cudaFreeHost (ch));
  //free (rh);
  //free (ch);

  tc[3].Stop ();



  for (int g = 0; g < NGPU; ++g) {
    cudaSetDevice (g);
    cudaDeviceReset ();
  }



#include "../frontend/server2.C"


  printf ("client mem alloc time\t\t%8.3f\n", tc[0].Span ());
  printf ("client compute time\t\t%8.3f\n", tc[1].Span ());
  printf ("client mem free time\t\t%8.3f\n", tc[3].Span ());


  return 0;
}


