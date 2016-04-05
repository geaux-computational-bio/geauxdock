#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <size.h>
#include <toggle.h>
#include <geauxdock.h>


struct DataSort
{
    int id;
    double val;
};

static bool
DataSortInc (DataSort a, DataSort b)
{
	return (a.val < b.val);
}


void
CopyLigandPoint (const Ligand0 * lig0, Ligand * lig, const int i0, const int i)
{
  lig->x1[i] = lig0->coord_orig.x[i0];
  lig->y1[i] = lig0->coord_orig.y[i0];
  lig->z1[i] = lig0->coord_orig.z[i0];
  lig->t1[i] = lig0->t[i0];
  lig->c1[i] = lig0->c[i0];
}


void
OptimizeLigand (const Ligand0 * lig0, const Kde * kde, Ligand * lig, const int n_lig)
{
  for (int i = 0; i < n_lig; ++i) {
    const Ligand0 *src = &lig0[i];
    Ligand *dst = &lig[i];
    const int lig_natom = src->lig_natom;

    dst->center[0] = src->pocket_center[0];
    dst->center[1] = src->pocket_center[1];
    dst->center[2] = src->pocket_center[2];
    dst->lig_natom = lig_natom;


    // keep ligand points in increament order of "n"


    // sort ligand in increament order of: t 
    // this violate the order of "n"
    std::vector <DataSort> ds;
    for (int j = 0; j < lig_natom; ++j) {
      DataSort d = {j, src->t[j]};
      ds.push_back (d);
    }

#if 1
    // PRT compute:
    // CPU/MIC: sort
    // GPU:     no sort
#if TARGET_DEVICE == TARGET_CPU || TARGET_DEVICE == TARGET_MIC
    //if (i == 0)
    //printf ("sort ligand\n");
    std::sort (ds.begin(), ds.end(), DataSortInc);
#endif
#endif
#if 0
    printf ("index:   ");
    for (int j = 0; j < lig_natom; ++j)
      printf ("%2d ", j);
    putchar ('\n');
    printf ("before:  ");
    for (int j = 0; j < lig_natom; ++j)
      printf ("%2d ", src->t[j]);
    putchar ('\n');
    printf ("after:   ");
    for (int j = 0; j < lig_natom; ++j)
      printf ("%2d ", src->t[ds[j].id]);
    printf ("\n\n");
#endif



    // xyz, might sort
    for (int j = 0; j < lig_natom; ++j)
      CopyLigandPoint (src, dst, ds[j].id, j);



    // xyz2, no sort
    for (int j = 0; j < lig_natom; ++j) {
      dst->x2[j] = src->coord_orig.x[j];
      dst->y2[j] = src->coord_orig.y[j];
      dst->z2[j] = src->coord_orig.z[j];
      dst->t2[j] = src->t[j];
      dst->c2[j] = src->c[j];
      //dst->c2[j] = 1.1f; // mode
    }

    /*
    if (i == 0)
        for (int j = 0; j < lig_natom; ++j)
            printf ("lig_c[%2d] = %26.18f\n", j, src->c[j]);
    */



    // sort ligand in increament order of t 
    std::vector <DataSort> ds2;
    for (int j = 0; j < lig_natom; ++j) {
      DataSort d = {j, src->t[j]};
      ds2.push_back (d);
    }
    std::sort (ds2.begin(), ds2.end(), DataSortInc);
    for (int j = 0; j < lig_natom; ++j) {
        const int j2 = ds2[j].id;
        dst->x3[j] = src->coord_orig.x[j2];
        dst->y3[j] = src->coord_orig.y[j2];
        dst->z3[j] = src->coord_orig.z[j2];
        dst->t3[j] = src->t[j2];
        dst->c3[j] = src->c[j2];
    }
  }











  // index for access KDE points


/*
  find kde_begin_idx[], and kde_end_idx[]
  kde is always sort by t

  state machine; 0-begin, 1-match, 2-quit

  0  0  0  0  0  1  1  1  1  2  2  2  2
  --------------============-----------
               begin        end
*/

  for (int i = 0; i < n_lig; ++i) {
      Ligand *mylig = &lig[i];
      const int lig_natom = mylig->lig_natom;

      for (int l = 0; l < lig_natom; ++l) {
          mylig->kde_begin_idx[l] = 0;
          mylig->kde_end_idx[l] = 0;
          int match_kde_status = 0;
          const int lig__t = mylig->t3[l];
          const int kde_npoint = kde->kde_npoint;

#if 0
          if (i == 0)
              printf ("lig_t %2d %2d  ", l, lig__t);
#endif
          for (int k = 0; k < kde_npoint; ++k) {
              const int kde__t = kde->t[k];
              if (lig__t == kde__t && match_kde_status == 0) {
                  match_kde_status = 1;
                  mylig->kde_begin_idx[l] = k;
              }
              if (lig__t != kde__t && match_kde_status == 1) {
                  match_kde_status = 2;
                  mylig->kde_end_idx[l] = k;
              }
          }
#if 0
          if (i == 0)
              printf ("kde %2d: %d %d\n", l, mylig->kde_begin_idx[l], mylig->kde_end_idx[l]);
#endif
      }
  }


}










void
CopyProteinPoint (const Protein0 * prt0, Protein * prt,
		  const int i0, const int i, const EnePara0 * enepara0)
{
  prt->x[i] = prt0->x[i0];
  prt->y[i] = prt0->y[i0];
  prt->z[i] = prt0->z[i0];
  prt->t[i] = prt0->t[i0];
  prt->c[i] = prt0->c[i0];

  int t = prt0->t[i0];
  int d = prt0->d[i0];
  int dt = t == 0 ? d + 30 : t;

  prt->ele[i] = enepara0->ele[dt];
  prt->seq3r[i] = prt0->seq3r[i0];
  prt->cdc[i] = prt0->cdc[i0];
  prt->hpp[i] = enepara0->hpp[d];

}

void
OptimizeProtein (Protein0 * prt0, Protein * prt, const EnePara0 * enepara0,
		 const int n_prt)
{
  for (int i = 0; i < n_prt; ++i) {
    Protein0 *src = &prt0[i];
    Protein *dst = &prt[i];
    const int prt_npoint = src->prt_npoint;
    dst->prt_npoint = prt_npoint;

    for (int j = 0; j < prt_npoint; ++j) {
      src->cdc[j] = ((src->c[j] == 0 && src->d[j] == 12) || (src->c[j] == 2));
      src->seq3r[j] = src->seq3[src->r[j]];
    }



    // sort protein in increament order of: c
    std::vector <DataSort> ds;
    for (int j = 0; j < prt_npoint; ++j) {
      DataSort d = {j, src->c[j]};
      ds.push_back (d);
    }


#if 1
    //if (i == 0)
    //printf ("sort protein\n");
    std::sort (ds.begin(), ds.end(), DataSortInc);
#endif
#if 0
    printf ("index:   ");
    for (int j = 0; j < prt_npoint; ++j)
      printf ("%2d ", j);
    putchar ('\n');
    printf ("before:  ");
    for (int j = 0; j < prt_npoint; ++j)
      printf ("%2d ", src->seq3r[j]);
    putchar ('\n');
    printf ("after:   ");
    for (int j = 0; j < prt_npoint; ++j)
      printf ("%2d ", src->seq3r[ds[j].id]);
    printf ("\n\n");
#endif



    for (int j = 0; j < prt_npoint; ++j)
      CopyProteinPoint (src, dst, ds[j].id, j, enepara0);
  }
}










void
OptimizePsp (const Psp0 * psp0, Psp * psp) // const Ligand * lig, const Protein * prt)
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
CopyKdePoint (const Kde0 * kde0, Kde * kde, const int i0, const int i)
{
  kde->x[i] = kde0->x[i0];
  kde->y[i] = kde0->y[i0];
  kde->z[i] = kde0->z[i0];
  kde->t[i] = kde0->t[i0];
}


void
OptimizeKde (const Kde0 * kde0, Kde * kde)
{
  const int kde_npoint = kde0->kde_npoint;
  kde->kde_npoint = kde_npoint;


  // sort kde in increament order of: t
  std::vector <DataSort> ds;
  for (int j = 0; j < kde_npoint; ++j) {
    DataSort d = {j, kde0->t[j]};
    ds.push_back (d);
  }


#if 1
  //printf ("sort kde\n");
  std::sort (ds.begin(), ds.end(), DataSortInc);
#endif
#if 0
    printf ("index:   ");
    for (int j = 0; j < kde_npoint; ++j)
      printf ("%2d ", j);
    putchar ('\n');
    printf ("before:  ");
    for (int j = 0; j < kde_npoint; ++j)
      printf ("%2d ", kde0->t[j]);
    putchar ('\n');
    printf ("after:   ");
    for (int j = 0; j < kde_npoint; ++j)
      printf ("%2d ", kde0->t[ds[j].id]);
    printf ("\n\n");
#endif


  for (int j = 0; j < kde_npoint; ++j)
    CopyKdePoint (kde0, kde, ds[j].id, j);
}





// sorting MCS is incorrect
void
OptimizeMcs (const Mcs0 * mcs0, Mcs * mcs, Mcs_R * mcs_r, Mcs_ELL * mcs_ell, Mcs_CSR *mcs_csr, const int nrow)
{
    //printf ("mcs_nrow = %d, MAX_MCS_COL= %d\n", nrow, MAX_MCS_COL); // nrow 11, ncol 128


#if 1
    // mcs
    for (int i = 0; i < nrow; ++i) {
        const Mcs0 *src = &mcs0[i];
        Mcs *dst = &mcs[i];
        for (int j = 0; j < MAX_MCS_COL; ++j) {
            dst->x[j] = src->x[j];
            dst->y[j] = src->y[j];
            dst->z[j] = src->z[j];
        }
        mcs[i].ncol = mcs0[i].ncol;
        mcs[i].tcc = mcs0[i].tcc;
    }
#endif


#if 0
    // mcs_r
    for (int i = 0; i < nrow; ++i) {
        const Mcs0 *src = &mcs0[i];
        for (int j = 0; j < MAX_MCS_COL; ++j) {
            mcs_r->x[j][i] = src->x[j];
            mcs_r->y[j][i] = src->y[j];
            mcs_r->z[j][i] = src->z[j];
        }
        mcs_r->tcc[i] = src->tcc;
    }
#endif


    // mcs_ell
#if 1
    for (int i = 0; i < nrow; ++i) {
        const Mcs0 *src = &mcs0[i];
        const int ncol = src->ncol;
        for (int j = 0; j < ncol; ++j) {
            mcs_ell->i[i][j] = src->idx_col[j];
            mcs_ell->x[i][j] = src->x2[j];
            mcs_ell->y[i][j] = src->y2[j];
            mcs_ell->z[i][j] = src->z2[j];
        }
        mcs_ell->ncol[i] = src->ncol;
        mcs_ell->tcc[i] = src->tcc;
    }
#endif





#if 1
    // mcs_csr, mcs_coo
    int row_ptr[MAX_MCS_ROW]; // row pointer for CSR sparse matrix
    row_ptr[0] = 0;
    for (int i = 0; i < nrow; ++i) {
        const Mcs0 *src = &mcs0[i];
        const int ncol = src->ncol;
        for (int j = 0; j < ncol; ++j) {
            const int linear_ptr = row_ptr[i] + j;
            mcs_csr->idx_col[linear_ptr] = src->idx_col[j];
            mcs_csr->idx_row[linear_ptr] = i;
            mcs_csr->x[linear_ptr] = src->x2[j];
            mcs_csr->y[linear_ptr] = src->y2[j];
            mcs_csr->z[linear_ptr] = src->z2[j];
        }
        mcs_csr->tcc[i] = src->tcc;
        row_ptr[i + 1] = row_ptr[i] + ncol;
    }
    mcs_csr->npoint = row_ptr[nrow];

/*
    mcs_csr->nrow = nrow;
    mcs_csr->row_ptr[0] = row_ptr[0];
    for (int i = 0; i < nrow; ++i) {
        mcs_csr->row_ptr[i + 1] = row_ptr[i + 1];
*/



#if 0
    int npoint = row_ptr[nrow];
    printf ("%d mcs rows\n", nrow);
    printf ("%d total mcs points\n", npoint);
    for (int i = 0; i < npoint; ++i) {
        printf ("%2d ", mcs_csr->i[i]);
    }
    printf ("\n");
#endif
#endif
}













#if 1
void
OptimizeEnepara (const EnePara0 * enepara0, EnePara * enepara)
{
  const float sqrt_2_pi = sqrtf (2.0f * PI);

  for (int i = 0; i < MAXTP2; ++i) {	// lig
    for (int j = 0; j < MAXTP1; ++j) {	// prt
      const float tmp = enepara0->vdw[j][i][0] * enepara0->lj[2];
      enepara->p1a[i][j] = 2.0f * enepara0->vdw[j][i][1] * powf (tmp, 9.0f);
      enepara->p2a[i][j] = 3.0f * enepara0->vdw[j][i][1] * powf (tmp, 6.0f);
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
      enepara->pmf0[i][j] = enepara0->pmf[j][i][0];
      enepara->pmf1[i][j] = enepara0->pmf[j][i][1];
      enepara->hdb0[i][j] = enepara0->hdb[j][i][0];
      enepara->hdb1[i][j] = 1.0f / enepara0->hdb[j][i][1];
    }
  }

  for (int i = 0; i < MAXTP4; ++i)
    enepara->hpp[i] = enepara0->hpp[i];
  for (int i = 0; i < MAXTP2; ++i) {
    enepara->hpl0[i] = enepara0->hpl[i][0];
    enepara->hpl1[i] = enepara0->hpl[i][1];
    enepara->hpl2[i] = logf (1.0f / (enepara->hpl1[i] * sqrt_2_pi));
  }

  enepara->kde2 = -0.5f / (enepara0->kde * enepara0->kde);
  enepara->kde3 = powf (enepara0->kde * sqrt_2_pi, 3.0f);

  for (int i = 0; i < MAXWEI; ++i) {
    enepara->w[i] = enepara0->w[i];
    // cout << enepara->w[i] << endl;
  }

  for (int i = 0; i < MAXWEI; ++i) {
    enepara->a_para[i] = enepara0->a_para[i];
    enepara->b_para[i] = enepara0->b_para[i];
    // cout << enepara->w[i] << endl;
  }
}
#endif


#if 0
void
OptimizeEnepara (const EnePara0 * enepara0, EnePara * enepara)
{
  const float sqrt_2_pi = sqrtf (2.0f * PI);

  for (int i = 0; i < MAXTP2; ++i) {	// lig
    for (int j = 0; j < MAXTP1; ++j) {	// prt
      const float tmp = enepara0->vdw[j][i][0] * enepara0->lj[2];
      enepara->p1a[i][j] = 2.0f * enepara0->vdw[i][j][1] * powf (tmp, 9.0f);
      enepara->p2a[i][j] = 3.0f * enepara0->vdw[i][j][1] * powf (tmp, 6.0f);
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
      enepara->pmf0[i][j] = enepara0->pmf[i][j][0];
      enepara->pmf1[i][j] = enepara0->pmf[i][j][1];
      enepara->hdb0[i][j] = enepara0->hdb[i][j][0];
      enepara->hdb1[i][j] = 1.0f / enepara0->hdb[i][j][1];
    }
  }

  for (int i = 0; i < MAXTP4; ++i)
    enepara->hpp[i] = enepara0->hpp[i];
  for (int i = 0; i < MAXTP2; ++i) {
    enepara->hpl0[i] = enepara0->hpl[i][0];
    enepara->hpl1[i] = enepara0->hpl[i][1];
    enepara->hpl2[i] = logf (1.0f / (enepara->hpl1[i] * sqrt_2_pi));
  }

  enepara->kde2 = -0.5f / (enepara0->kde * enepara0->kde);
  enepara->kde3 = powf (enepara0->kde * sqrt_2_pi, 3.0f);

  for (int i = 0; i < MAXWEI; ++i) {
    enepara->w[i] = enepara0->w[i];
    // cout << enepara->w[i] << endl;
  }

  for (int i = 0; i < MAXWEI; ++i) {
    enepara->a_para[i] = enepara0->a_para[i];
    enepara->b_para[i] = enepara0->b_para[i];
    // cout << enepara->w[i] << endl;
  }
}
#endif

