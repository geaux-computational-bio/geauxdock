#!/bin/sh -x
mpirun \
-np 1 ../src/main_mpi/server -nt 20 : \
-np 1 ../src/main_mpi/client_gpu
#-np 1 ../src/main_mpi/client_mic
#-np 1 ../src/main_mpi/client_cpu
#-np 1 ../src/main_mpi/client_dummy : \


