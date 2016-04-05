#include <mpi.h>

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

#include <yeah/c/timing.h>
#include <yeah/cpp/timer.hpp>



MPI_Datatype MPI_JOB;



void Server (int argc, char **argv, const int id, const int nclient)
{
#include "../frontend/server1.C"

  // send jobs

  MPI_Status status;
  MPI_Request req[ncomplex];
  MPI_Status st[ncomplex];

  // dispatch jobs on requests
  for (int i = 0; i < ncomplex; ++i) {
    Complex *job = complex + i;
    //Complex *job = complex + 4;
    int dst;
    const int mytag = 0;
    MPI_Recv (&dst, 1, MPI_INT, MPI_ANY_SOURCE, mytag, MPI_COMM_WORLD, &status);
    MPI_Isend (job, 1, MPI_JOB, dst, mytag, MPI_COMM_WORLD, &req[i]);
    printf ("%-20s server sent job %02d to client %02d\n", argv[0], job->signal, dst);
  }
  for (int i = 0; i < ncomplex; ++i) {
    MPI_Wait (&req[i], &st[i]);
  }

  // send finish signals
  Complex * done = complex + ncomplex - 1; // send the last complex in the queue
  done->signal = FINISH_SIGNAL;
  for (int i = 0; i < nclient; ++i) {
    int dst;
    const int mytag = 0;
    MPI_Recv (&dst, 1, MPI_INT, MPI_ANY_SOURCE, mytag, MPI_COMM_WORLD, &status);
    MPI_Isend (done, 1, MPI_JOB, dst, mytag, MPI_COMM_WORLD, &req[i]);
    printf ("%-20s server kill client %02d\n", argv[0], dst);
  }
  for (int i = 0; i < ncomplex; ++i) {
    MPI_Wait (&req[i], &st[i]);
  }

#include "../frontend/server2.C"
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


  if (id == 0)
    Server (argc, argv, id, nprocs - 1);
  else
    printf ("server is assigned with a wrong MPI rank\n");


  return 0;
}
