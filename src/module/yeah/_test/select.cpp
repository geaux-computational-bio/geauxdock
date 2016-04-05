#include <iostream>
#include <yeah/cpp/select.hpp>

int main ()
{
    int n = 7;
    double a[] = {0, 1, 2, 3, 4, 5, 6, 7};
    double x;

    //x = yeah::Min_array <double> (&a[0], &a[n]);
    //x = yeah::Max_array <double> (&a[0], &a[n]);
    //x = yeah::Avg_array <double> (&a[0], &a[n]);
    //x = yeah::Avg_array <double> (&a[0], &a[0]);


    //x = yeah::Min_array_exclude <double> (&a[0], &a[n], 0);
    //x = yeah::Max_array_exclude <double> (&a[0], &a[n], 7);
    //x = yeah::Max_array_exclude <double> (&a[0], &a[n], 5);
    //x = yeah::Min_array_exclude <double> (&a[0], &a[n], 3);
    //x = yeah::Avg_array_exclude <double> (&a[0], &a[n], 7);
    //x = yeah::Avg_array_exclude <double> (&a[0], &a[n], 0);
    x = yeah::Avg_array_exclude <double> (&a[1], &a[1], 1);

    std::cout << x << std::endl;
    return 0;
}



