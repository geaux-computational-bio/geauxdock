#include <mpi.h>
#include <cstdio>

#include <geauxdock.h>

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

#include "../backend_gpu/prepair.h"
#include "../backend_gpu/backend_gpu.h"
#include <yeah/cpp/timer.hpp>
#include <yeah/cuda/runtime/wrapper.h>


MPI_Datatype MPI_JOB;




void Client (int argc, char **argv, const int id)
{
  yeah::Timer t[2];
  t[0].Start ();

  // host data
  Complex *ch;
  Record *rh;
  CUDA_ERR (cudaMallocHost ((void **) &ch, sizeof (Complex)));
  CUDA_ERR (cudaMallocHost ((void **) &rh, sizeof (Record) * MAX_REP));


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


  int send_msg = id;
  const int dst = 0;
  MPI_Status status;
  const int mytag = 0;

  while (1) {
    MPI_Send (&send_msg, 1, MPI_INT, dst, mytag, MPI_COMM_WORLD);
    printf ("%-20s \t\t\t\t\t client %02d request\n", argv[0], id);
    MPI_Recv (ch, 1, MPI_JOB, dst, mytag, MPI_COMM_WORLD, &status);

    if (ch->signal == FINISH_SIGNAL) {
      printf ("%-20s \t\t\t\t\t\t\t\t\t\t client %02d retired\n", argv[0], id);
      break;
    }
    else {
      printf ("%s start docking\n", argv[0]);
      SetParaT (ch, pt);
      CopyH2D (ch, cd, pt);

      t[1].Start ();
      Dock (ch, rh, cd, rd, pt, sd);
      t[1].Stop ();
    }
  }

  DeviceFree (cd, rd, sd);
  CUDA_ERR (cudaFreeHost (rh));
  CUDA_ERR (cudaFreeHost (ch));


  for (int g = 0; g < NGPU; ++g) {
    cudaSetDevice (g);
    cudaDeviceReset ();
  }




  t[0].Stop ();
  printf ("client:\n");
  printf ("client wall time\t\t%8.3f\n", t[0].Span ());
  printf ("Run MC time time\t\t%8.3f\n", t[1].Span ());
}




int
main (int argc, char **argv)
{
  int id, nprocs;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  if (nprocs < 2) {
    printf ("nprocs < 2. exit\n");
    exit (-1);
  }

  MPI_Type_contiguous (sizeof (Complex), MPI_BYTE, &MPI_JOB);
  MPI_Type_commit (&MPI_JOB);


  if (id != 0)
    Client (argc, argv, id);
  else
    printf ("client is assigned with a wrong MPI rank\n");


  return 0;
}
 
