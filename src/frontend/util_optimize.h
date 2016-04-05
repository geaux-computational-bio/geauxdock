#ifndef  UTIL_OPTIMIZE_H
#define  UTIL_OPTIMIZE_H

void OptimizeLigand (const Ligand0 *, const Kde *, Ligand *, const int);
void OptimizeProtein (Protein0 *, Protein *, const EnePara0 *, const int);
void OptimizePsp (const Psp0 *, Psp *);
void OptimizeKde (const Kde0 *, Kde *);
void OptimizeMcs (const Mcs0 *, Mcs *, Mcs_R *, Mcs_ELL *, Mcs_CSR *, const int);
void OptimizeEnepara (const EnePara0 *, EnePara *);

#endif

