#include <iostream>
#include <yeah/c/timing.h>
#include <yeah/cpp/timer.hpp>

#define N 100

void
func1 ()
{
    yeah::TimerAuto t_auto ("func1()");
    double a = 0.234;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                a += 0.0001234;
            }
        }
    }
    std::cout << a << std::endl;
}
int main ()
{
    yeah::Timer t0;
    yeah::TimerAuto t_auto ("main()");

    t0.Start ();
    func1 ();
    t0.Stop ();
    t0.Print ();
    //std::cout << HostTimeNow () << std::endl;

    return 0;
}

