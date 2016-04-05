
#ifndef _YEAH_PAPI_STRUCT_H_
#define _YEAH_PAPI_STRUCT_H_

typedef struct {
  const int event;
  const char * str;
} Papi_event_struct;


#define PAPI_NULL_YEAH 0x00000000



int Papi_event_struct_n (Papi_event_struct * tt)
{
    int n = 0;
    while (1) {
        if (tt[n].event != PAPI_NULL_YEAH)
            n++;
        else
            break;
    }
    return n;
}




#endif
