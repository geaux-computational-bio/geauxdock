#include <papi.h>


/* Return current time (number of milli-seconds since the Epoch). */
double
HostTimeNow ()
{
    return (double) PAPI_get_virt_usec () / 1000.0f;
}


