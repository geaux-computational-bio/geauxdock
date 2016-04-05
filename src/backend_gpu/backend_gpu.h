#ifndef BACKEND_GPU_H
#define BACKEND_GPU_H

void InitCurand (curandState **s);

void Dock (Complex *ch,
           Record *rh,
           Complex **cd,
           Record **rd,
           ParaT **pt,
           curandState **curandstate_d);


#endif


