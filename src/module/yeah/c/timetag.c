#include <time.h>


void GetTimeTag (char * s, const int len)
{
    time_t t1 = time (NULL);
    struct tm *t2 = localtime (&t1);
    strftime (s, len, "%Y%m%d_%H%M%S", t2);
}



/*
void xxxxx (char * s, const int len)
{
    time_t t1 = time (NULL);
    struct tm *t2 = localtime (&t1);
    printf ("%s", asctime (t2));
}
*/
