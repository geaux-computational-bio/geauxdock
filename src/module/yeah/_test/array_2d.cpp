#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <yeah/cpp/array_2d.hpp>

int
main (int argc, char ** argv)
{
#define SZ 8

    int **a = yeah::Alloc_2d <int> (SZ, SZ);
    a[SZ - 1][SZ - 1] = 9;
    yeah::Print_2d <int> (a, SZ, SZ);
    printf ("\n");
    yeah::Print_2d_v2 <int> (a[0], SZ, SZ);
    yeah::Free_2d <int> (a);

    return 0;
}

