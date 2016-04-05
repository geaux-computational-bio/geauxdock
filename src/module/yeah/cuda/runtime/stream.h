#ifndef _YEAH_CUDA_RUNTIME_STREAM_H_
#define _YEAH_CUDA_RUNTIME_STREAM_H_

#include <cuda_runtime.h>
#include <iostream>


namespace yeah {

    namespace cuda {

        class Stream
        {
        public:
            cudaStream_t s;

            Stream () {}
            ~Stream () {}
            void Sync () { cudaStreamSynchronize (s); }
            void SyncEvent (cudaEvent_t e) { cudaStreamWaitEvent (s, e, 0); }
        };


        class StreamSD: public Stream
        {
        public:
            StreamSD () { cudaStreamCreate (&s); }
            ~StreamSD () { cudaStreamDestroy (s); }
        };


        class StreamMD: public Stream
        {
        public:
            void Create () { cudaStreamCreate (&s); }
            void Destroy () { cudaStreamDestroy (s); }
        };

    }
}

#endif
