// Intel MKL VSL Random Number Generator
#include <mkl_vsl.h>

float MyRand_d ()
{
  float randdd = (float) rand () / RAND_MAX;
  //float randdd = 0.002f;
  return randdd;
}


#define MYRAND (MyRand_d ())





void
MonteCarlo_d (ParaD pd, const int s1, const int s2max)
{

#pragma omp parallel for
  for (int r = pd.rep_begin; r <= pd.rep_end; ++r) {
    // GPU shared vars
    ReplicaMC *rep = &pd.replica[r];
    const Ligand * const lig = &pd.lig[rep->idx_lig];
    const Protein * const prt = &pd.prt[rep->idx_prt];
    const Psp * const psp = pd.psp;
    const Kde * const kde = pd.kde;
    const Mcs * const mcs = pd.mcs;
    const EnePara * const enepara = pd.enepara;



    float movematrix[6]; // translation x y z, rotation x y z
    float rot[3][3]; // rotz roty rotx
    float ligxyz[MAXLIG][3];


    float enepara_p12[MAXTP2][MAXTP1][2];
    float enepara_pmf[MAXTP2][MAXTP1][2];
    float enepara_hdb[MAXTP2][MAXTP1][2];
    float enepara_hpl[MAXTP2][3];

    float enepara_lj0;
    float enepara_lj1;
    float enepara_el1;
    float enepara_el0;
    float enepara_a1;
    float enepara_b1;
    float enepara_kde2;
    float enepara_kde3;

    float sqrt_2_pi_m1; // -1.0f / sqrtf (2.0f * PI)

    int lig_natom, prt_npoint, kde_npoint, mcs_nrow;
    int is_accept_s;


    for (int l = 0; l < MAXTP2; ++l) {
      for (int p = 0; p < MAXTP1; ++p) {
	enepara_p12[l][p][0] = enepara->p12[l][p][0];
	enepara_p12[l][p][1] = enepara->p12[l][p][1];
	enepara_pmf[l][p][0] = enepara->pmf[l][p][0];
	enepara_pmf[l][p][1] = enepara->pmf[l][p][1];
	enepara_hdb[l][p][0] = enepara->hdb[l][p][0];
	enepara_hdb[l][p][1] = enepara->hdb[l][p][1];
      }
    }

    for (int l = 0; l < MAXTP2; ++l) {
      for (int k = 0; k < MAXTP2; ++k) {
	enepara_hpl[l][k] = enepara->hpl[l][k];
      }
    }

    enepara_lj0 = enepara->lj0;
    enepara_lj1 = enepara->lj1;
    enepara_el1 = enepara->el1;
    enepara_el0 = enepara->el0;
    enepara_a1 = enepara->a1;
    enepara_b1 = enepara->b1;
    enepara_kde2 = enepara->kde2;
    enepara_kde3 = enepara->kde3;

    sqrt_2_pi_m1 = -1.0f / sqrtf (2.0f * PI);
    lig_natom = lig->lig_natom;
    prt_npoint = prt->prt_npoint;
    kde_npoint = pd.kde_npoint;
    mcs_nrow = pd.mcs_nrow;
    is_accept_s = rep->is_accept;

    pd.record[r - pd.rep_begin].next_entry = 0; // reset the record's next_entry




    for (int s2 = 0; s2 < s2max; ++s2) {

      ////////////////////////////////////////////////////////////////////////
      // record old states
      if (is_accept_s == 1) {
	rep->step = s1 + s2;

	const int rr = r - pd.rep_begin;
	const int next_entry = pd.record[rr].next_entry;
	pd.record[rr].replica[next_entry] = *rep;
	pd.record[rr].next_entry = next_entry + 1;
      }



      ////////////////////////////////////////////////////////////////////////
      // move

      //#pragma unroll (6)
      for (int bidx = 0; bidx < 6; ++bidx) {
#if IS_AWAY == 0
	const float fixed_var = 0.0f;
#elif IS_AWAY == 1
	const float fixed_var = 44.5f;
#endif
	float moveamount = s2max != 1 ? MYRAND : fixed_var;
	movematrix[bidx] = pd.move_scale[bidx] * moveamount + rep->movematrix[bidx];
      }

      // http://en.wikipedia.org/wiki/Euler_angles
      // http://upload.wikimedia.org/math/e/9/c/e9cf817bce9c1780216921cd93233459.png
      // http://upload.wikimedia.org/math/f/4/e/f4e55dcos2c9581007648967d29b15121e.png
      const float sin1 = sinf (movematrix[3]);
      const float cos1 = cosf (movematrix[3]);
      const float sin2 = sinf (movematrix[4]);
      const float cos2 = cosf (movematrix[4]);
      const float sin3 = sinf (movematrix[5]);
      const float cos3 = cosf (movematrix[5]);
      rot[0][0] = cos1 * cos2;
      rot[0][1] = cos1 * sin2 * sin3 - cos3 * sin1;
      rot[0][2] = sin1 * sin3 + cos1 * cos3 * sin2;
      rot[1][0] = cos2 * sin1;
      rot[1][1] = cos1 * cos3 + sin1 * sin2 * sin3;
      rot[1][2] = cos3 * sin1 * sin2 - cos1 * sin3;
      rot[2][0] = -1 * sin2;
      rot[2][1] = cos2 * sin3;
      rot[2][2] = cos2 * cos3;

      const float cx = lig->center[0];
      const float cy = lig->center[1];
      const float cz = lig->center[2];

      // rotation and translation, and apply coordinate system transformation
      for (int l = 0; l < lig_natom; ++l) {
	const float x = lig->a[l].x;
	const float y = lig->a[l].y;
	const float z = lig->a[l].z;
	ligxyz[l][0] = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z + movematrix[0] + cx;
	ligxyz[l][1] = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z + movematrix[1] + cy;
	ligxyz[l][2] = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z + movematrix[2] + cz;
      }





      ////////////////////////////////////////////////////////////////////////
      // calcenergy

#if 1

      float evdw = 0.0f;
      float eele = 0.0f;
      float epmf = 0.0f;
      float epsp = 0.0f;
      float ehdb = 0.0f;
      float ehpc = 0.0f;
      float ekde = 0.0f;
      float elhm = 0.0f;


      // lig loop, ~30
      //#pragma unroll 8
      for (int l = 0; l < lig_natom; ++l) {
	float ehpc1 = 0.0f;
	const int lig_t = lig->a[l].t;

	// prt loop, ~300
	//#pragma unroll (4)
	for (int p = 0; p < prt_npoint; ++p) {
	  const int prt_t = prt->a[p].t;
	
	  const float dx = ligxyz[l][0] - prt->a[p].x;
	  const float dy = ligxyz[l][1] - prt->a[p].y;
	  const float dz = ligxyz[l][2] - prt->a[p].z;
	  const float dst_pow2 = dx * dx + dy * dy + dz * dz;
	  const float dst_pow4 = dst_pow2 * dst_pow2;
	  const float dst = sqrtf (dst_pow2);


#if 1
	  /* hydrophobic potential */
	  if (prt->a[p].c0_and_d12_or_c2 == 1 && dst_pow2 <= 81.0f) {
	    ehpc1 += prt->a[p].hpp *
	      (1.0f - (3.5f / 81.0f * dst_pow2 -
		       4.5f / 81.0f / 81.0f * dst_pow4 +
		       2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
		       0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
	  }
#endif
      

#if 1
	  // p1a[MAXTP2][MAXTP1]
	  // p2a[MAXTP2][MAXTP1]
      
	  /* L-J potential */
	  const float p1 = enepara_p12[lig_t][prt_t][0] / (dst_pow4 * dst_pow4 * dst);
	  const float p2 = enepara_p12[lig_t][prt_t][1] / (dst_pow4 * dst_pow2);
	  const float p4 = p1 * enepara_lj0 * (1.0f + enepara_lj1 * dst_pow2) + 1.0f;
	  evdw += (p1 - p2) / p4;
#endif




#if 1
	  /* electrostatic potential */
	  const float s1 = enepara_el1 * dst;
	  float g1;
	  if (s1 < 1)
	    g1 = enepara_el0 + enepara_a1 * s1 * s1 + enepara_b1 * s1 * s1 * s1;
	  else
	    g1 = 1.0f / s1;
	  eele += lig->a[l].c * prt->a[p].ele * g1;
#endif


#if 1
	  // pmf0[MAXTP2][MAXTP1]
	  // pmf1[MAXTP2][MAXTP1]
	  // psp[MAXTP2][MAXPRO]

	  /* contact potential */
	  const float dst_minus_pmf0 = dst - enepara_pmf[lig_t][prt_t][0];

	  epmf +=
	    enepara_pmf[lig_t][prt_t][1] /
	    (1.0f + expf ((-0.5f * dst + 6.0f) * dst_minus_pmf0));


	  /* pocket-specific potential */
	  // the senmatics do not match with the original program:
	  // if (found psp[][])
	  //   accumulate to epsp;
	  // else
	  //   do nothing
	  if (prt->a[p].c == 2 && dst_minus_pmf0 <= 0) {
	    const int i1 = prt->a[p].seq3r;
	    epsp += psp->psp[lig_t][i1]; // sparse matrix
	  }
#endif


#if 1
	  // hdb0[MAXTP2][MAXTP1]
	  // hdb1[MAXTP2][MAXTP1]

	  /* hydrogen bond potential */
	  const float hdb0 = enepara_hdb[lig_t][prt_t][0];
	  if (hdb0 > 0.1f) {
	    const float hdb1 = enepara_hdb[lig_t][prt_t][1];
	    const float hdb3 = (dst - hdb0) * hdb1;
	    ehdb += hdb1 * expf (-0.5f * hdb3 * hdb3);
	  }
#endif

	} // prt loop


#if 1
	// not performance critical
	// hpl0[MAXTP2]
	// hpl1[MAXTP2]
	// hpl2[MAXTP2]

	/* hydrophobic restraits*/
	// transpose may help improve the performance
	const float hpc2 = (ehpc1 - enepara_hpl[lig_t][0]) / enepara_hpl[lig_t][1];
	ehpc += 0.5f * hpc2 * hpc2 - enepara_hpl[lig_t][2];
#endif

      } // lig loop
#endif

      ehdb *= sqrt_2_pi_m1;



#if 1
      /* kde potential */

      // lig loop, ~30
      for (int l = 0; l < lig_natom; ++l) {
	float kde_val = 0.0f;
	int kde_sz = 0;
	const int lig_t = lig->a[l].t;

	// kde loop, ~400
	for (int k = 0; k < kde_npoint; ++k) {
	  if (lig_t == kde->a[k].t) {
	    const float dx = ligxyz[l][0] - kde->a[k].x;
	    const float dy = ligxyz[l][1] - kde->a[k].y;
	    const float dz = ligxyz[l][2] - kde->a[k].z;
	    const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
	    kde_val += expf (enepara_kde2 * kde_dst_pow2);
	    kde_sz++;
	  }
	} // kde loop

	if (kde_sz != 0)
	  ekde += (kde_val / (float) kde_sz);

      } // lig loop

      ekde /= (enepara_kde3 * lig_natom);
#endif


#if 1

      /* position restraints */




      // lhm loop, ~11
      //#pragma unroll 2
      for (int m = 0; m < mcs_nrow; ++m) {
	float lhm_val = 0.0f;
	int lhm_sz = 0;

	// lig loop, ~30
	for (int l = 0; l < lig_natom; ++l) {
	  const int lig_n = lig->a[l].n + 1;
	  if (mcs[m].a[lig_n].x != MCS_INVALID_COORD) {
	    const float dx = ligxyz[l][0] - mcs[m].a[lig_n].x;
	    const float dy = ligxyz[l][1] - mcs[m].a[lig_n].y;
	    const float dz = ligxyz[l][2] - mcs[m].a[lig_n].z;
	    lhm_val += dx * dx + dy * dy + dz * dz;
	    lhm_sz++;
	  }
	} // lig loop

	if (lhm_sz != 0)
	  elhm += mcs[m].tcc * sqrtf (lhm_val / (float) lhm_sz);

      } // lhm loop
      elhm = logf (elhm / mcs_nrow);

#endif





      // energy edst e[8]
#if 1
      const float dx = lig->center[0] + movematrix[0]- prt->pocket_center[0];
      const float dy = lig->center[1] + movematrix[1]- prt->pocket_center[1];
      const float dz = lig->center[2] + movematrix[2]- prt->pocket_center[2];
      const float edst = sqrtf (dx * dx + dy * dy + dz * dz);
#endif




      float e[MAXWEI];
      e[0] = evdw / lig_natom; // 0 - vdw 
      e[1] = eele / lig_natom; // 1 - ele
      e[2] = epmf / lig_natom; // 2 - pmf (CP)
      e[3] = epsp / lig_natom; // 3 - psp (PS CP)
      e[4] = ehdb / lig_natom; // 4 - hdb (HB)
      e[5] = ehpc / lig_natom; // 5 - hpc (HP)
      e[6] = ekde; // 6 - kde (PHR)
      e[7] = elhm; // 7 - lhm (MCS)
      e[8] = edst; // 8 - dst (DST)

      // normalization
      for (int i = 0; i < MAXWEI - 1; ++i)
	e[i] = enepara->ab_para[i][0] * e[i] + enepara->ab_para[i][1];

      // calculate the total energy from energy terms using linear model
      float etotal = 0.0f;
      for (int i = 0; i < MAXWEI - 1; ++i)
	etotal += enepara->w[i] * e[i];
      e[MAXWEI - 1] = etotal;




      ////////////////////////////////////////////////////////////////////////
      // accept
      const float delta_energy = etotal - rep->energy[MAXWEI -1];
      const float beta = pd.temp[rep->idx_tmp].minus_beta;
      float rand = s2max != 1 ? MYRAND : 0.0f; // force to accept if s2max == 1
      const int is_accept = (rand < expf (delta_energy * beta));  // mybeta < 0

      if (is_accept == 1) {
	//#pragma unroll MAXWEI
	for (int i = 0; i < MAXWEI; ++i)
	  rep->energy[i] = e[i];
	//#pragma unroll (6)
	for (int i = 0; i < 6; ++i)
	  rep->movematrix[i] = movematrix[i];
      }
      is_accept_s = is_accept;

    } // s2 loop

    rep->is_accept = is_accept_s;


  } // replica loop

}

