#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>


// return current time (number of mili-seconds since the Epoch).
double
HostTimeNow2 ()
{
  struct timeval t;
  gettimeofday (&t, NULL);
  double mytime_second = (double) t.tv_sec + (double) t.tv_usec * 1e-6;
  return mytime_second * 1000;
}



// link with -lrt
double
HostTimeNow ()
{
  struct timespec t; 
  clock_gettime (CLOCK_REALTIME, &t);
  double mytime_second = (double) t.tv_sec + (double) t.tv_nsec * 1e-9;
  return mytime_second * 1000;
}

