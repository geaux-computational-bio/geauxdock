#ifndef _YEAH_PAPI_MEASURER_H_
#define _YEAH_PAPI_MEASURER_H_


#include <cstdio>
#include <cstdlib>
#include <papi.h>
#include <yeah/papi/wrapper.h>
#include <yeah/papi/struct.h>



/*
// usage:
Papi_event_struct papi_event_struct1[] = {
    //{PAPI_L1_TCM,       "PAPI_L1_TCM"},
    //{PAPI_L2_TCM,       "PAPI_L2_TCM"},
    {PAPI_L1_DCM,       "PAPI_L1_DCM"},
    //{PAPI_L2_DCM,       "PAPI_L2_DCM"},
    //{PAPI_L1_ICM,       "PAPI_L1_ICM"},
    //{PAPI_L2_ICM,       "PAPI_L2_ICM"},

    //{PAPI_LD_INS,       "PAPI_LD_INS"},
    //{PAPI_SR_INS,       "PAPI_SR_INS"},


    //{PAPI_BR_MSP,       "PAPI_BR_MSP"},
    //{PAPI_BR_PRC,       "PAPI_BR_PRC"},
    //{PAPI_BR_INS,       "PAPI_BR_INS"},


    //{PAPI_SP_OPS,       "PAPI_SP_OPS"},
    //{PAPI_DP_OPS,       "PAPI_DP_OPS"},
    //{PAPI_FP_OPS,       "PAPI_FP_OPS"},
    //{PAPI_VEC_SP,       "PAPI_VEC_SP"},
    //{PAPI_VEC_DP,       "PAPI_VEC_DP"},


    {PAPI_TOT_INS,        "PAPI_TOT_INS"},
    {PAPI_TOT_CYC,        "PAPI_TOT_CYC"},
    {PAPI_NULL_YEAH,      "PAPI_NULL_YEAH"}
};
*/





namespace yeah {
    namespace papi {

        class Measurer0
        {
        public:
            Papi_event_struct * papi_event_struct;
            int papi_event_n;
            int *papi_event;
            long long *papi_event_val;
            float t0, t1; // time

            void Init (Papi_event_struct * tt)
            {
                papi_event_struct = tt;
                papi_event_n = Papi_event_struct_n (papi_event_struct);
                papi_event_val = new long long[papi_event_n];
                papi_event = new int[papi_event_n];
                for (int i = 0; i < papi_event_n; ++i) {
                    papi_event[i] = papi_event_struct[i].event;
                }
            }

            void Shutdown ()
            {
                PAPI_shutdown ();
                delete[] papi_event_val;
                delete[] papi_event;
            }

            void StartCounters ()
            {
                t0 = (float) PAPI_get_virt_usec () / 1000000.0f;
                PAPI_ERR (PAPI_start_counters (papi_event, papi_event_n));
            }

            void ReadCounters ()
            {
                t1 = (float) PAPI_get_virt_usec () / 1000000.0f;
                PAPI_ERR (PAPI_read_counters (papi_event_val, papi_event_n));
            }

            void Print ()
            {
                printf ("time\t\t%10.3f s\n", t1 - t0);
                for (int i = 0; i < papi_event_n; ++i) {
                    printf ("%s\t%10.3f M\n",
                        papi_event_struct[i].str,
                        (float) papi_event_val[i] / 1000000.0f);
                }
                printf ("\n");
            }

            void Start ()
            {
                StartCounters ();
            }

            void Stop ()
            {
                ReadCounters ();
            }

        };


        // usage:
        // yeah::papi::Measurer m0 (papi_event_struct1);
        // m0.Start ();
        // funcx ();
        // m0.Stop ();
        // m0.Print ();
        // m0.Shutdown ();

        class Measurer: public Measurer0
        {
        public:
            Measurer (Papi_event_struct * papi_event_struct1)
            {
                Init (papi_event_struct1);
            }
            ~Measurer ()
            {
            }
        };


        // usage:
        // yeah::papi::MeasurerAuto m0 (papi_event_struct1);
        // funcx ();

        class MeasurerAuto: public Measurer0
        {
        public:
            MeasurerAuto (Papi_event_struct * papi_event_struct1)
            {
                Init (papi_event_struct1);
                Start ();
            }

            ~MeasurerAuto ()
            {
                Stop ();
                Print ();
                Shutdown ();
            }
        };

    }
}



#endif
