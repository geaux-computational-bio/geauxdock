#ifndef RECORD_H
#define RECORD_H


#include "size.h"




struct Energy
{
  float e[MAXWEI];
  /*
     0 - vdw
     1 - ele
     2 - pmf
     3 - psp
     4 - hdb

     5 - hpc
     6 - kde
     7 - lhm

     8 - dst

     9 - total
   */
  float cms;
  float rmsd;
};


struct Replica
{
  int idx_rep; // n_rep, replica
  int idx_prt; // n_prt, protein
  int idx_tmp; // n_tmp, temperature
  int idx_lig; // n_lig, ligand
};


struct LigRecordSingleStep
{
  Replica replica;
  Energy energy;
  float movematrix[6]; // // translation x y z, rotation x y z
  int step;
};


/*
struct Medoid
{
  LigRecordSingleStep step;
  int cluster_sz;
  float relative_cluster_sz;
  float z_score;
  float ave_cluster_ener;
  float config_integral;
};
*/



#endif /* RECORD_H */
