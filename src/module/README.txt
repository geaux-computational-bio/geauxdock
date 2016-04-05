drop this whole directory ("module") under "src_root" of your project


code:
  #include <yeah/c/timing.h>
  ...

build:
  cc -I./module -c a.c
  cc a.o -L./module/lib -lyeahc



Please read the copyright of each components before distributing the code.
For what I understand, you should announce the following statement(s):

This software contains source code provided by NVIDIA Corporation.
This software contains source code (Random123) developed by D. E. Shaw Research.
This software contains source code Catch developed by Phil Nash


