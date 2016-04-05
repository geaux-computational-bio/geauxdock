
#ifndef _YEAH_CPP_ARRAY_2D_HPP_
#define _YEAH_CPP_ARRAY_2D_HPP_

#include <iostream>

namespace yeah {

// allocate a contiguous region of 2D array
    template <class T>
    T ** Alloc_2d (const int rows, const int cols)
    {
        //T *data = (T *) malloc (sizeof (T) * rows * cols);
        //T **a= (T **) malloc (sizeof (T *) * rows);
        T *data = new T[rows * cols];
        T **a= new T*[rows];

        for (int r = 0; r < rows; ++r)
            a[r] = &(data[cols * r]);

        return a;
    }


    template <class T>
    void Free_2d (T** a)
    {
        //free (a[0]);
        //free (a);
        delete[] a[0];
        delete[] a;
    }


// print_2d_v2 (a, rows, cols)
    template <class T>
    void Print_2d (T** a, const int rows, const int cols)
    {
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                std::cout << a[r][c] << " ";
            }
            std::cout << std::endl;
        }
    }


// print_2d_v2 (a[0], rows, cols)
    template <class T>
    void Print_2d_v2 (T* a, const int rows, const int cols)
    {
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                std::cout << a[cols * r + c] << " ";
            }
            std::cout << std::endl;
        }
    }



}




#endif
