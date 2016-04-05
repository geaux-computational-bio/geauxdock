#ifndef _YEAH_CUDA_RUNTIME_EVENT_H_
#define _YEAH_CUDA_RUNTIME_EVENT_H_

#include <cuda_runtime.h>
#include <iostream>



namespace yeah {

namespace cuda {

  class Event
  {
  public:
    cudaEvent_t e1; // start
    cudaEvent_t e2; // stop
    cudaEvent_t e3; // set point
    double ts; // span (miliseconds)
    int n; // counter

    Event () { ts = 0; n = 0; }
    ~Event () {} 

    void Start (cudaStream_t s = 0) { cudaEventRecord (e1, s); }
    void Stop (cudaStream_t s = 0) { cudaEventRecord (e2, s); }
    void Calc () { float tt; cudaEventElapsedTime (&tt, e1, e2); ts += (double) tt; n++; }
    void SyncStop () { cudaEventSynchronize (e2); }

    void Set (cudaStream_t s = 0) { cudaEventRecord (e3, s); }
    void Sync () { cudaEventSynchronize (e3); }
    cudaError_t Query () { return cudaEventQuery (e3); }

    void Flush () { ts = 0; n = 0;}
    float Span () { return ts; }
    int Count () { return n; }
    void Print () { std::cout << "span (ms): " << ts << "\tcounter: " << n << std::endl; }
  };


  class EventSD: public Event
  {
  public:
    EventSD () { cudaEventCreate (&e1); cudaEventCreate (&e2); cudaEventCreate (&e3); }
    ~EventSD () { cudaEventDestroy (e1); cudaEventDestroy (e2); cudaEventDestroy (e3); } 
  };


  class EventMD: public Event
  {
  public:
    void Create () { cudaEventCreate (&e1); cudaEventCreate (&e2); cudaEventCreate (&e3); }
    void Destroy () { cudaEventDestroy (e1); cudaEventDestroy (e2); cudaEventDestroy (e3); } 
  };


}
}







#endif

