// do not use flag -std=c99

#include <stdio.h>
#include <time.h>

// milli m 1e-3
// micro u 1e-6
// nano  n 1e-9


void
Sleep_milliseconds (double ms)
{
  const long us = (double) (ms * 1e3);
  const long sec_ = (us / 1000000);
  const long nsec_ = (us % 1000000) * 1000;

  struct timespec t1, t2;
  t1.tv_sec = sec_;
  t1.tv_nsec = nsec_;
  if (nanosleep (&t1, &t2) < 0)
    printf("nanosleep () failed\n");
  //else
  //  printf("sleep %ds %dns\n", sec_, nsec_);
}


// an example of random sleep
//Sleep_milliseconds (rand () % 2000));


