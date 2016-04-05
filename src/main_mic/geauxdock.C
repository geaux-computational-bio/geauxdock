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

#include "../backend_cpu_mic/backend_cpu_mic.h"
#include <yeah/c/timing.h>
#include <yeah/cpp/timer.hpp>






int main (int argc, char **argv)
{

#include "../frontend/server1.C"


  Record *record = (Record *) malloc (sizeof (Record) * MAX_REP);
  for (int i = 0; i < ncomplex; ++i) {
    Complex *job = complex + i;
    printf ("%s start docking\n", argv[0]);
    Dock (job, record);
  }
  free (record);


#include "../frontend/server2.C"

  return 0;
}


