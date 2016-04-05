#include <mpi.h>
#include <cstdio>
#include <cstdlib>

#include <geauxdock.h>

#include "../backend_cpu_mic/backend_cpu_mic.h"
#include <yeah/cpp/timer.hpp>


MPI_Datatype MPI_JOB;



void Client (int argc, char **argv, const int id)
{
  yeah::Timer t[2];
  t[0].Start ();

  // data
  Complex *recv_msg = (Complex *) malloc (sizeof (Complex));
  Record *record = (Record *) malloc (sizeof (Record) * MAX_REP);


  int send_msg = id;
  const int dst = 0;
  MPI_Status status;
  const int mytag = 0;

  while (1) {
    MPI_Send (&send_msg, 1, MPI_INT, dst, mytag, MPI_COMM_WORLD);
    printf ("%-20s \t\t\t\t\t client %02d request\n", argv[0], id);
    MPI_Recv (recv_msg, 1, MPI_JOB, dst, mytag, MPI_COMM_WORLD, &status);

    if (recv_msg->signal == FINISH_SIGNAL) {
      printf ("%-20s \t\t\t\t\t\t\t\t\t\t client %02d retired\n", argv[0], id);
      break;
    }
    else {
      printf ("%s start docking\n", argv[0]);
      t[1].Start ();
      Dock (recv_msg, record);
      t[1].Stop ();
    }
  }

  free (record);
  free (recv_msg);


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
 
