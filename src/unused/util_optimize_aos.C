#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "size.h"
#include "toggle.h"
#include "geauxdock.h"
#include "util_optimize.h"



void
OptimizeLigand (const Ligand0 * lig0, Ligand * lig, const ComplexSize complexsize)
{

  // data structure translation
  for (int i = 0; i < complexsize.n_lig; ++i) {
    const Ligand0 *src = &lig0[i];
    Ligand *dst = &lig[i];

    for (int l = 0; l < MAXLIG; ++l) {
      LigandPoint* lp = &dst->a[l];
      lp->x = src->coord_orig.x[l];
      lp->y = src->coord_orig.y[l];
      lp->z = src->coord_orig.z[l];
      lp->t = src->t[l];
      lp->c = src->c[l];
      lp->n = src->n[l];
    }

    dst->center[0] = src->pocket_center[0];
    dst->center[1] = src->pocket_center[1];
    dst->center[2] = src->pocket_center[2];

    dst->lig_natom = src->lig_natom;
  }


  // to do??
  // sort lig in increament order of t
}


void
CopyProteinResidue (const Protein0 * src, Protein * dst, const int residue_src,
		    const int residue_dst, const EnePara0 * enepara0)
{
  ProteinPoint* pp = &dst->a[residue_dst];

  pp->x = src->x[residue_src];
  pp->y = src->y[residue_src];
  pp->z = src->z[residue_src];

  pp->t = src->t[residue_src];
  pp->c = src->c[residue_src];

  int d = src->d[residue_src];
  int dt = pp->t == 0 ? d + 30 : pp->t;
  pp->ele = enepara0->ele[dt];

  pp->seq3r = src->seq3[src->r[residue_src]];

  pp->c0_and_d12_or_c2 =
    (src->c[residue_src] == 0 && src->d[residue_src] == 12) || (src->c[residue_src] == 2);

  pp->hpp = enepara0->hpp[d];
}

void
OptimizeProtein (const Protein0 * prt0, Protein * prt, const EnePara0 * enepara0,
		 const Ligand0 * lig0, const ComplexSize complexsize)
{
  // pocket center
  const float cx = lig0[0].pocket_center[0];
  const float cy = lig0[0].pocket_center[1];
  const float cz = lig0[0].pocket_center[2];

  for (int i = 0; i < complexsize.n_prt; ++i) {
    const Protein0 *src = &prt0[i];
    Protein *dst = &prt[i];
    const int prt_npoint = src->prt_npoint;

    // sort protein in increament order of t
    int *t = (int *) malloc (sizeof (int) * prt_npoint);
    int *order = (int *) malloc (sizeof (int) * prt_npoint);
    for (int j = 0; j < prt_npoint; ++j) {
      t[j] = src->t[j];
      order[j] = j;
    }

#if 0
    Sort (t, order, 0, prt_npoint - 1);
#endif

#if 0
    for (int j = 0; j < prt_npoint; ++j)
      printf ("%d ", t[j]);
    putchar ('\n');

    for (int j = 0; j < prt_npoint; ++j)
      printf ("%d ", src->t[order[j]]);
    putchar ('\n');
#endif

    dst->prt_npoint = prt_npoint;
    for (int j = 0; j < prt_npoint; ++j)
      CopyProteinResidue (src, dst, order[j], j, enepara0);

    free (t);
    free (order);

    // assign the pocket center from the ligand structure to the protein sturcture
    dst->pocket_center[0] = cx;
    dst->pocket_center[1] = cy;
    dst->pocket_center[2] = cz;

  }

}





void
OptimizePsp (const Psp0 * psp0, Psp * psp, const Ligand * lig, const Protein * prt)
{
  for (int i = 0; i < MAXPRO; ++i) {
    for (int j = 0; j < MAXLIG; ++j) {
      psp->psp[j][i] = psp0->psp[i][j];
    }
  }

/*
  int lig_lna = lig[0].lig_natom;
  int prt_pnp = prt[0].prt_npoint;
  int count = 0;
*/

/*
  for (int i = 0; i < lig_lna; ++i) {
    int lig_t = lig[0].t[i];
    for (int j = 0; j < prt_pnp; ++j) {
      int pspidx2 = prt[0].seq3r[j];

      if (prt[0].c[j] == 2) {
	printf ("potentialy accessed psp \t%2d\t%2d\n", lig_t, pspidx2);
	count++;

      if (psp->psp[i][j] != 0)
        printf ("lig %3d \tprt %3d \t\tpsp\t%f\n", i, j, psp->psp[j][i]);

      }
    }
  }
*/

/*
  printf ("percentage %f / %f = %f", count, lig_lna * prt_pnp,
	  (float) count / (lig_lna * prt_pnp));
*/

}

void
OptimizeKde (const Kde0 * kde0, Kde * kde)
{
  for (int i = 0; i < MAXKDE; ++i) {
    KdePoint* kp = &kde->a[i];
    kp->x = kde0->x[i];
    kp->y = kde0->y[i];
    kp->z = kde0->z[i];
    kp->t = kde0->t[i];
  }
  kde->kde_npoint = kde0->kde_npoint;
}


void
OptimizeMcs (const Mcs0 * mcs0, Mcs * mcs, const ComplexSize complexsize)
{
  // mcs_nrow
  for (int i = 0; i < complexsize.mcs_nrow; ++i) {
    mcs[i].tcc = mcs0[i].tcc;

    // n_mcs
    for (int j = 0; j < MAX_MCS_COL; ++j) {
      McsPoint* mp = &mcs[i].a[j];
      mp->x = mcs0[i].x[j];
      mp->y = mcs0[i].y[j];
      mp->z = mcs0[i].z[j];
    }
  }

}

void
OptimizeEnepara (const EnePara0 * enepara0, EnePara * enepara)
{
  const float sqrt_2_pi = sqrtf (2.0f * PI);

  for (int i = 0; i < MAXTP2; ++i) {	// lig
    for (int j = 0; j < MAXTP1; ++j) {	// prt
      const float tmp = enepara0->vdw[j][i][0] * enepara0->lj[2];
      enepara->p12[i][j][0] = 2.0f * enepara0->vdw[j][i][1] * powf (tmp, 9.0f);
      enepara->p12[i][j][1] = 3.0f * enepara0->vdw[j][i][1] * powf (tmp, 6.0f);
    }
  }
  enepara->lj0 = enepara0->lj[0];
  enepara->lj1 = enepara0->lj[1];

  enepara->el1 = enepara0->el[1];
  enepara->el0 = enepara0->el[0];
  enepara->a1 = 4.0f - 3.0f * enepara0->el[0];
  enepara->b1 = 2.0f * enepara0->el[0] - 3.0f;


  for (int i = 0; i < MAXTP2; ++i) {	// lig
    for (int j = 0; j < MAXTP1; ++j) {	// prt
      enepara->pmf[i][j][0] = enepara0->pmf[j][i][0];
      enepara->pmf[i][j][1] = enepara0->pmf[j][i][1];
      enepara->hdb[i][j][0] = enepara0->hdb[j][i][0];
      enepara->hdb[i][j][1] = 1.0f / enepara0->hdb[j][i][1];
    }
  }

  for (int i = 0; i < MAXTP4; ++i)
    enepara->hpp[i] = enepara0->hpp[i];
  for (int i = 0; i < MAXTP2; ++i) {
    enepara->hpl[i][0] = enepara0->hpl[i][0];
    enepara->hpl[i][1] = enepara0->hpl[i][1];
    enepara->hpl[i][2] = logf (1.0f / (enepara->hpl[i][1] * sqrt_2_pi));
  }

  enepara->kde2 = -0.5f / (enepara0->kde * enepara0->kde);
  enepara->kde3 = powf (enepara0->kde * sqrt_2_pi, 3.0f);

  for (int i = 0; i < MAXWEI; ++i) {
    enepara->w[i] = enepara0->w[i];
    // cout << enepara->w[i] << endl;
  }

  for (int i = 0; i < MAXWEI; ++i) {
    enepara->ab_para[i][0] = enepara0->a_para[i];
    enepara->ab_para[i][1] = enepara0->b_para[i];
    // cout << enepara->w[i] << endl;
  }
}


