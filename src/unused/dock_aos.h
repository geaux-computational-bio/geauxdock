#ifndef  DOCK_AOS_H
#define  DOCK_AOS_H

// array of struct


struct LigandPoint
{
  float x, y, z;
  int t;
  float c;
  int n;
};


struct MYALIGN Ligand
{
  LigandPoint a[MAXLIG];

  // coord_center is under the lab system
  float center[3];		// ligand geometric center
  int lig_natom;			// number of ligand atoms
};



struct ProteinPoint
{
  float x, y, z;
  int t;
  int c;
  float ele;
  int seq3r;
  int c0_and_d12_or_c2;
  float hpp;
};



struct MYALIGN Protein
{
  ProteinPoint a[MAXPRO];
  float pocket_center[3];
  int prt_npoint;
};




struct MYALIGN Psp
{
  float psp[MAXLIG][MAXPRO];                                           // replaced
  //float sparse_psp[MAXLIG][MAXPRO];                                    // replaced
};


struct KdePoint
{
  float x, y, z;
  int t;
};

struct MYALIGN Kde
{
  KdePoint a[MAXKDE];
  int kde_npoint;			// number of kde points
};


struct McsPoint
{
  float x, y, z;
};

struct MYALIGN Mcs
{
  McsPoint a[MAX_MCS_COL];
  float tcc;
};



struct MYALIGN EnePara
{
  // L-J
  float p12[MAXTP2][MAXTP1][2];
  float lj0, lj1;

  // electrostatic
  float el1;
  float el0;
  float a1; // 4.0f - 3.0f * el0;
  float b1; // 2.0f * el0 - 3.0f;

  // contact
  float pmf[MAXTP2][MAXTP1][2];  // interchange to [lig][prt]
  float hdb[MAXTP2][MAXTP1][2];  // interchange to [lig][prt]

  // hydrophobic
  float hpp[MAXTP4];
  float hpl[MAXTP2][3];

  // kde
  float kde2; // -0.5f / (kde * kde)
  float kde3; // powf (kde * sqrtf (2.0f * PI), 3.0f)

  // weights for energy terms
  float ab_para[MAXWEI][2];         // the a,b parameter for normalization
  float w[MAXWEI];
};



#endif // DOCK_SOA_H

