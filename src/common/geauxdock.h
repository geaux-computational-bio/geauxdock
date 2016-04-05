#ifndef  DOCK_H
#define  DOCK_H


#include <string>

#include "size.h"
#include "toggle.h"



#if defined(__CUDACC__) // NVCC
   #define ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
  #define ALIGN(n) __attribute__((aligned(n)))
#elif defined(__INTEL_COMPILER) // INTEL
  #define ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
  #define ALIGN(n) __declspec(align(n))
#else
  #error "Please provide a definition for MY_ALIGN macro for your host compiler!"
#endif


//#define MYALIGN ALIGN(32)
//#define MYALIGN ALIGN(128)

// CPU AVX
#define MYALIGN ALIGN(64)

#define RESTRICT __restrict__



struct TraceFile
{
  std::string path;  // ligand trace path
};

struct LigandFile
{
  std::string path; // lig file path
  std::string conf_path;
  std::string id;	// ligand name
  std::string molid; // MOLID;

  int conf_total;	// number of conformation for one lig
  int raw_conf;		// raw total conf in the .sdf without small rmsd excluded

  int lig_natom;	// total effective pts number
  int lnb;	// total bonds number

  //int ensemble_total;	// ENSEMBLE_TOTAL in the .sdf file
};


struct ProteinFile
{
  std::string path;  // prt file path
  std::string id;    // prt name

  int conf_total;	// number of conformations for one lig
  int prt_npoint;	// total effective pts number
  int pnr;	// total bonds number
};


struct LhmFile
{
  std::string path;
  std::string ligand_id;
  int mcs_nrow; // number of mcs positions
};


struct EneParaFile
{
  std::string path;
};

struct NorParaFile
{
  std::string path_a;
  std::string path_b;
};

struct WeightFile
{
  std::string path;
};

struct InputFiles
{
  std::string lig_list;
  std::string prt_list;

  int lig_files_num;
  LigandFile * lig_files;
  LhmFile * lhm_files;

  ProteinFile prt_file;
  EneParaFile enepara_file;
  WeightFile weight_file;
  NorParaFile norpara_file;
  TraceFile trace_file;
};






struct ComplexSize
{
  // replica numbers
  int n_lig; // number of ligand conf
  int n_prt; // number of protein conf
  int n_tmp; // number of temperature
  int n_rep; // n_rep = n_lig * n_prt * n_tmp;

  // residue numbers (of per replica)
  int lig_natom; // number of ligand points
  int prt_npoint; // number of protein points
  int kde_npoint; // number of kde points
  int mcs_nrow; // number of mcs positions
};




// Energy
//float e[MAXWEI];
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




struct LigCoord
{
  float x[MAXLIG];		// Ligand x coords
  float y[MAXLIG];		// Ligand y coords
  float z[MAXLIG];		// Ligand z coords
  float center[3];		// ligand geometric center
};




//#if TARGET_DEVICE == TARGET_GPU
#include "geauxdock_soa.h"
//#elif TARGET_DEVICE == TARGET_CPU || TARGET_DEVICE == TARGET_MIC
//#include "dock_aos.h"
//#endif







struct Ligand0
{
  LigCoord coord_orig;           //                                      used

  int t[MAXLIG];		// atom type                            used
  float c[MAXLIG];		// atom charge                          used
  int n[MAXLIG];		// atom number                          after loading, lig0.n[i] == i + 1

  int lig_natom;			// number of ligand atoms               used
  int lnb;			// number of ligand bonds               NOT USED

  float pocket_center[3];	// pocket center                        used
                                // should belong to "Protein" structure


  float lens_rmsd;		// ensemble rmsd                        NOT USED

  float mw;			// ligand molecular weight              NOT USED
  float logp;			// ligand water/octanol partition coeff NOT USED
  float psa;			// ligand polar surface area            NOT USED
  float mr;			// ligand molar refractivity            NOT USED

  int hbd;			// ligand hydrogen bond donors          NOT USED
  int hba;			// ligand hydrogen bond acceptors       NOT USED

  std::string a[MAXLIG];	// atom name                            NOT USED
  std::string id;		// ligand id                            NOT USED
  std::string smiles;		// ligand smiles                        NOT USED
};



struct Protein0
{
  float x[MAXPRO];		// residue x coord                      used
  float y[MAXPRO];		// residue y coord                      used
  float z[MAXPRO];		// residue z coord                      used
  int n[MAXPRO];		// effective point number               NOT USED
  int t[MAXPRO];		// effective point type                 used
  int c[MAXPRO];		// effective point class                used
  int d[MAXPRO];		// redidue code                         used

  int prt_npoint;			// number of protein effective points   used
  int pnr;			// number of protein residues           NOT USED

  int r[MAXPRO];		// residue number                       replaced
  int seq3[MAXPRO];		// aa sequence numbering                replaced


  int cdc[MAXPRO]; // cdc = (prt_c == 0 && prt_d == 12) || (prt_c == 2)
  int seq3r[MAXPRO]; // seq3r[i] == prt->seq3[prt->r[i]]


  //std::string protein_seq1;  // aa sequence
  //char protein_seq2[MAXPRO]; // aa sequence
};




struct Psp0
{
  float psp[MAXPRO][MAXLIG];                                           // replaced
  int n;			// total number of PSP point           NOT USED
};



struct Kde0
{
  float x[MAXKDE];		// KDE x coord                          used
  float y[MAXKDE];		// KDE y coord                          used
  float z[MAXKDE];		// KDE z coord                          used
  int n[MAXKDE];		// KDE point number                     NOT USED
  int t[MAXKDE];		// KDE atom type                        used

  int kde_npoint;			// number of kde points                 used
  int pns[MAXTP2];		// number of specific kde points        NOT USED
};







struct Mcs0
{
  float x[MAX_MCS_COL];              //                                      used
  float y[MAX_MCS_COL];              //                                      used
  float z[MAX_MCS_COL];              //                                      used


  int   idx_col[MAX_MCS_COL];         // index for sparse matrix              used
  float x2[MAX_MCS_COL];              //                                      used
  float y2[MAX_MCS_COL];              //                                      used
  float z2[MAX_MCS_COL];              //                                      used
  int ncol;			// column number

  float tcc;                    //                                      used
};



struct EnePara0
{
  float vdw[MAXTP1][MAXTP2][2];	// L-J potential                        vdw[prt_t][lig_t][]
  float ele[MAXTP3];		// electrostatic potential              ele[prt_d + 30], ele[prt_t]
  float pmf[MAXTP1][MAXTP2][2];	// contact potential                    pmf[prt_t][lig_t][]
  float hpp[MAXTP4];		// protein hydrophobicity               hpp[prt_d]
  float hpl[MAXTP2][2];		// ligand hydrophobicity                hpl[lig_t][]
  float hdb[MAXTP1][MAXTP2][2];	// ligand hydrophobicity                hdb[prt_t][lig_t][]

  float lj[3];			// L-J params
  float el[2];			// electrostatic params
  float kde;			// kde bandwidth

  float w[MAXWEI];		// weights for energy terms
  float a_para[MAXWEI];         // the a parameter in normalization
  float b_para[MAXWEI];         // the b parameter in normalization
};













struct Temp
{
  float t;
  float minus_beta;
  int order;
};





// replica[n_rep]
// replica[n_prt][n_tmp][n_lig]

struct MYALIGN ReplicaMC
{
  int idx_rep; // n_rep, replica

  int idx_lig; // n_lig, ligand
  int idx_prt; // n_prt, protein
  int idx_tmp; // n_tmp, temperature

  float movematrix[6];       // translation x y z, rotation x y z
  float energy[MAXWEI];

  int step; // step counter
  int is_accept; // was the last purturb accepted? record data if is_accpet == 1
};



struct ExchgPara
{
  float floor_temp;   // lowest temperature in all replicas
  float ceiling_temp;   // highest temperature in all replicas
  int num_temp;      // number of temperatures for the same ligand and protein conformations
};


struct McPara
{
  int steps_total;
  int steps_per_dump;
  int steps_per_exchange;

  float move_scale[6]; // translation x y z, rotation x y z

  char outputdir[MAXSTRINGLENG];
  char outputfile[MAXSTRINGLENG];
};

struct McLog
{
  int ac_mc;  // MC acceptance counter
  int acs_mc[MAX_REP];   // MC acceptance counter for all replicas
  int ac_temp_exchg;
  int acs_temp_exchg[MAX_REP]; 
  
  // int ac_lig_exchg;
  // int acs_lig_exchg[MAX_REP];
};




struct Record
{
  ReplicaMC replica[STEPS_PER_DUMP];
  int next_entry;
};




struct MoveVector
{
  float ele[6]; // translation xyz, rotation xyz
};




struct Complex
{
  // GPU read only arrays
  Ligand lig[MAX_CONF_LIG];
  Protein prt[MAX_CONF_PRT];
  Psp psp;
  Kde kde;
  Mcs mcs[MAX_MCS_ROW];
  Mcs_R mcs_r;  
  Mcs_ELL mcs_ell;
  Mcs_CSR mcs_csr;
  EnePara enepara;
  Temp temp[MAX_TMP];

  // GPU writable arrays that duplicate on multiple GPUs
  ReplicaMC replica[MAX_REP];
  //float *etotal; // used by replica exchange
  //MoveVector *movevector; // used by replcia exchange

  // remove "record" pointers from Complex struct
  // due to the limitation of Xeon Phi offloading mode
  // GPU writable arrays that spreads over multiple GPUs
  // Record *record;


  // extra parameters
  ComplexSize size;
  McPara mcpara;


  // file names
  char lig_file[MAX_STR_LENG];
  char prt_file[MAX_STR_LENG];
  char lhm_file[MAX_STR_LENG];
  char enepara_file[MAX_STR_LENG];
  char weight_file[MAX_STR_LENG];


  // GPU read only scalars
  // sizes for multi-device decomposition
  int rep_begin;
  int rep_end;
  //int n_rep;
  //int record_sz;


  // MPI message signal
  int signal;
};







struct ParaT // Data transter parameters
{
  int rep_begin; // replica range
  int rep_end; // replica range
  int n_rep; // replica number

  int record_sz; // record size
};


#endif // DOCK_H

