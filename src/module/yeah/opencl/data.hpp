#ifndef _YEAH_OPENCL_DATA_HPP_
#define _YEAH_OPENCL_DATA_HPP_


#include <assert.h>
#include <yeah/opencl/cl.h>


namespace yeah {

    namespace opencl {

        template <class T0>
            class Data0
            {
            public:
                T0 *dh;
                cl_mem dd;
                size_t sz;
                int n;

                //Data0 (int num = 1) { this->Init(num); }
                //~Data0 () {}
                void Init (int num) { n = num; sz = sizeof (T0) * num; }

                void AllocH () { this->dh = (T0 *) malloc (this->sz); assert (this->dh != NULL); }
                void FreeH () { free (this->dh); }
                //GE (cudaMalloc ((void **) &this->dd, this->sz))
                //clReleaseMemObject (this->dd);
            };



        template <class T0>
            class Data: public Data0 <T0>
        {
        public:
            Data (int num = 1) { this->Init(num); }
            ~Data() {}
        };


    }
}

#endif
