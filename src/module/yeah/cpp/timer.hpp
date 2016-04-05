#ifndef _YEAH_CPP_TIMER_HPP_
#define _YEAH_CPP_TIMER_HPP_

#include <iostream>
#include <yeah/c/timing.h>


// try switch to std::time

namespace yeah {

    class Timer
    {
    protected:
        double t1; // start
        double t2; // stop
        double ts; // span
        int n; // counter

    public:
        Timer () { t1 = t2 = ts = 0; n = 0; }
        void Start () { t1 = HostTimeNow (); }
        void Stop () { t2 = HostTimeNow (); ts += t2 - t1; n++; }
        void Reset () { t1 = t2 = ts = 0; n = 0; }
        double Span () { return ts; }
        int Count () { return n; }
        void Print () { std::cout << "span (ms): " << ts <<"\tcounter: " << n << std::endl; }
    };




    // usage:
    // yeah::TimerAuto t ("func ()");

    class TimerAuto: public Timer
    {
    protected:
        std::string name;

    public:
        TimerAuto (const std::string & name1)
        {
            name = name1;
            Start ();
        }

        ~TimerAuto ()
        {
            Stop ();
            std::cout <<
                "TimerAuto: " <<
                name << ": " <<
                Span () << " ms" << std::endl;
        }
    };

}


#endif
