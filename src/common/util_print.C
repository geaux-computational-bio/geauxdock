#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <size.h>
#include <toggle.h>
#include <geauxdock.h>
#include <util_print.h>



using namespace std;



// arg = 1      print title
// arg = 2      print content
// arg = 3      print all

void
PrintStepTrack (const ReplicaMC *rep, const int arg)
{

  char names[11][30] = {
    "step",
    "lig_conf",
    "prt_conf",
    "temp_idx",
    "total",
    "movevector"
  };

  int a = arg & 0x1;
  int b = (arg >> 1) & 0x1;

  if (a == 1) {
    for (int i = 0; i < 11; ++i)
      printf (",%s", names[i]);
    printf ("\n");
  }

  if (b == 1) {
    printf (",%d", rep->step);
    printf (",%d", rep->idx_prt);
    printf (",%d", rep->idx_lig);
    printf (",%d", rep->idx_tmp);
    printf (",%f", rep->energy[MAXWEI - 1]);
    for (int i = 0; i < 6; i++)
      printf (",%f", rep->movematrix[i]);
    printf ("\n");
  }
}


void
PrintTrack (Record * record, int r, int iter_begin, int iter_end, int arg)
{
  // print title
  PrintStepTrack (NULL, 1);

  for (int s = iter_begin; s <= iter_end; ++s) {
    const ReplicaMC *rep = &record[r].replica[s];
    PrintStepTrack (rep, arg);
  }
}






void
PrintCsv (const float* energy, const int idx_rep, const int step, const int arg)
{
  char names[MAXWEI][8] = {
    "vdw", // 0
    "ele", // 1
    "pmf", // 2
    "psp", // 3
    "hdb", // 4
    "hpc", // 5
    "kde", // 6
    "lhm", // 7
    "dst", // 8
    "total" // 9
  };

  int a = arg & 0x1;
  int b = (arg >> 1) & 0x1;

  if (a == 1) {
    printf ("rep step");
    for (int i = 0; i < MAXWEI; ++i)
      printf (" %s", names[i]);
    printf ("\n");
  }

  if (b == 1) {
    printf ("%3d %d", idx_rep,  step);
    for (int i = 0; i < MAXWEI; ++i)
      // printf (" %+14.10f"  energy[i]);
      printf (" %.4f",  energy[i]);
    printf ("\n");
  }
}





void
PrintMoveVector (const float m[6], const int step)
{
  printf ("\t  %3d\t\t\t", step);
  for (int i = 0; i < 6; ++i) {
    printf (" %+f\t", m[i]);
  }
  printf ("\n");
}



void
PrintMoveRecord (const Record * record, const int steps_per_dump, const int r,
		 const int iter_begin, const int iter_end, const int arg)
{
  for (int s = iter_begin; s <= iter_end; ++s) {
    const ReplicaMC *rep = &record[r].replica[s];
    PrintMoveVector (rep->movematrix, rep->step);
  }

}


// arg = 1      print title
// arg = 2      print content
// arg = 3      print all

void
PrintRecord (Record * record, int steps_per_dump, int r, int iter_begin, int iter_end, int arg)
{
  // print title
  //PrintCsv (NULL, 0, 0, 1);

  for (int s = iter_begin; s <= iter_end; ++s) {
    const ReplicaMC *rep = &record[r].replica[s];
    PrintCsv (rep->energy, r, rep->step, arg);
  }

}





void
PrintRepRecord (const Record * record, const int steps_per_dump, const int rep_begin,
		const int rep_end, const int iter_begin, const int iter_end, const int arg)
{
  printf ("\treplicas\n");

  printf ("step|\t");

  for (int r = rep_begin; r <= rep_end; ++r)
    printf ("%2d\t", r);
  putchar ('\n');

  printf ("----+");

  for (int r = rep_begin; r <= rep_end; ++r)
    printf ("--------");
  putchar ('\n');

  for (int s = iter_begin; s <= iter_end; ++s) {
    printf ("%3d |\t", s);

    for (int r = rep_begin; r <= rep_end; ++r) {
      const ReplicaMC *myrep = &record[r].replica[s];
      //printf ("%2d ", myrep->idx_prt);
      //printf ("%2d ", myrep->idx_tmp);
      //printf ("%2d ", myrep->idx_lig);

      printf ("%2d ", myrep->idx_rep);

      printf ("\t");
    }
    putchar ('\n');
  }

}





// print all temperature replicas of the same lig & prt
void
PrintRepRecord2 (Record * record, ComplexSize complexsize,
		 int steps_per_dump, int idx_prt, int idx_lig,
		 int iter_begin, int iter_end, int arg)
{
  printf ("temperature replicas with lig %d prt %d\n", idx_lig, idx_prt);

  printf ("MC step |\t");

  for (int t = 0; t < complexsize.n_tmp; ++t)
    printf ("%2d\t", t);
  putchar ('\n');

  printf ("--------+----");

  for (int t = 0; t < complexsize.n_tmp; ++t)
    printf ("--------");
  putchar ('\n');

  for (int s = iter_begin; s <= iter_end; ++s) {
    printf ("%5d   |\t", s);

    for (int t = 0; t < complexsize.n_tmp; ++t) {
      const int r =
	complexsize.n_tmp * complexsize.n_lig * idx_prt + complexsize.n_lig * t + idx_lig;
      const ReplicaMC *myrep = &record[r].replica[s];
      //printf ("%2d ", myrep->idx_prt);
      printf ("%2d ", myrep->idx_tmp);
      //printf ("%2d ", myrep->idx_lig);
      //printf ("%2d ", myrep->idx_rep);

      printf ("\t");
    }
    putchar ('\n');
  }
}







/*
void
PrintLigand (const Ligand * lig)
{
  printf ("center:\t\t%+10.6f\t%+10.6f\t%+10.6f\n", lig->center[0], lig->center[1], lig->center[2]);
  printf ("lig_natom:\t\t%d\n", lig->lig_natom);

  printf ("x \t\ty \t\tz \t\tc \t\t t \t n \tindex\n");
  printf ("-----------------------------------------------\n");
  const int lig_natom = lig->lig_natom;
  for (int i = 0; i < lig_natom; ++i) {
    printf ("%+10.6f\t", lig->x[i]);
    printf ("%+10.6f\t", lig->y[i]);
    printf ("%+10.6f\t", lig->z[i]);
    printf ("%+10.6f\t", lig->c[i]);
    printf ("%2d\t", lig->t[i]);
    printf ("%2d\t", lig->n[i]);
    printf ("%3d\n", i);
  }
}
*/

/*
void
PrintProtein (const Protein * prt)
{
  printf ("prt_npoint:\t\t%d\n", prt->prt_npoint);

  printf ("x \t\ty \t\tz \t\t t \t c \t d \tindex\n");
  printf ("-----------------------------------------------\n");
  const int prt_npoint = prt->prt_npoint;
  for (int i = 0; i < prt_npoint; ++i) {
    printf ("%+10.6f\t", prt->x[i]);
    printf ("%+10.6f\t", prt->y[i]);
    printf ("%+10.6f\t", prt->z[i]);
    printf ("%2d\t", prt->t[i]);
    printf ("%2d\t", prt->c[i]);
    printf ("%4d\n", i);
  }

}
*/


/*
void
PrintDataSize (const ParaH ph)
{
  float lig_sz = sizeof (Ligand) / 1024;
  float prt_sz = sizeof (Protein) / 1024;
  float psp_sz = sizeof (Psp) / 1024;
  float kde_sz = sizeof (Kde) / 1024;
  float mcs_sz = sizeof (Mcs) / 1024;
  float enepara_sz = sizeof (EnePara) / 1024;

  printf ("lig \t\tprt \t\tpsp \t\tkde \t\tmcs \t\tenepara\n");
  printf ("%f \t%f \t%f \t%f \t%f \t%f\t\t",
	  lig_sz, prt_sz, psp_sz, kde_sz, mcs_sz, enepara_sz);
  printf ("KB per struct\n");


  printf ("%f \t%f \t%f \t%f \t%f \t%f\t\t",
	  lig_sz * ph.complexsize->n_lig,
          prt_sz * ph.complexsize->n_prt,
          psp_sz,
          kde_sz,
          mcs_sz * ph.complexsize->mcs_nrow,
	  enepara_sz);
  printf ("KB total\n");
}
*/





void
PrintSummary (const Complex * c)
{
  putchar ('\n');

  printf ("===============================================================\n");
  printf ("Inputs and Outputs\n");
  printf ("===============================================================\n");
  printf ("ligand file\t\t\t%s\n", c->lig_file);
  printf ("protein file\t\t\t%s\n", c->prt_file);
  printf ("lhm file\t\t\t%s\n", c->lhm_file);
  printf ("enepara file\t\t\t%s\n", c->enepara_file);
  printf ("weight file\t\t\t%s\n", c->weight_file);
  
  printf ("output directory\t\t%s\n", c->mcpara.outputdir);
  printf ("out file (HDF5)\t\t\t%s/%s_XXXX.h5\n", c->mcpara.outputdir, c->mcpara.outputfile);

  printf ("steps_per_dump\t\t\t%d\n", c->mcpara.steps_per_dump);

  printf ("===============================================================\n");
  printf ("Sizes\n");
  const float sz_prt = (float) sizeof (Protein) * MAX_CONF_PRT / 1024 / 1024;
  const float sz1 = (float) sizeof (Complex) / 1024 / 1024;
  const float sz2 = (float) sizeof (Record) * MAX_REP / 1024 / 1024;
  const float sz3 = (float) sizeof (Record) * c->size.n_rep / 1024 / 1024;
  printf ("size of protein\t\t\t\t\t%.3f MB\n", sz_prt);
  printf ("size of each complex\t\t\t\t%.3f MB\n", sz1);
  printf ("record size (memory allocation):\t\t%.3f MB\n", sz2);
  printf ("record size (memory size of each dump file):\t%.3f MB\n", sz3);
  printf ("===============================================================\n");


  printf ("Monte Carlo parameters\n");
  printf ("===============================================================\n");
  printf ("steps_total\t\t\t%d\n", c->mcpara.steps_total);
  printf ("steps_per_dump\t\t\t%d\n", c->mcpara.steps_per_dump);


  printf ("translational scale\t\t");
  for (int i = 0; i < 3; ++i)
    printf ("%.8f ", c->mcpara.move_scale[i]);
  printf ("\n");

  printf ("rotational scale\t\t");
  for (int i = 3; i < 6; ++i)
    printf ("%.8f ", c->mcpara.move_scale[i]);
  printf ("\n");

  printf ("ligand conformations\t\t%d\n", c->size.n_lig);
  printf ("prt conformations\t\t%d\n", c->size.n_prt);
  printf ("temperatures\t\t\t%d\n", c->size.n_tmp);
  printf ("replica ensembles\t\t%d\n", c->size.n_rep);

  printf ("size_lig\t\t\t%d\n", c->size.lig_natom);
  printf ("size_prt\t\t\t%d\n", c->size.prt_npoint);
  printf ("size_pnk\t\t\t%d\n", c->size.kde_npoint);
  printf ("size_mcs\t\t\t%d\n", c->size.mcs_nrow);
  printf ("===============================================================\n");

}


