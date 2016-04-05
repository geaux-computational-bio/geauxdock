ty = threadIdx.x / bdx_prt;
tx = threadIdx.x % bdx_prt;



// lig loop, ~30
//#pragma unroll 8
for (int j = 0; j < lig_natom; j += bdy_prt) { // y loop
    float ehpc1 = 0.0f;

    const int l = j + ty;
    if (l < lig_natom) {
        const int lig__t = lig_t[l];

//#pragma unroll 4
        // prt loop, ~300
        for (int p = tx; p < prt_npoint; p += bdx_prt) { // x loop
            const int prt__t = CUDA_LDG_D (prt->t[p]);

            const float dx = lig_x2[l] - CUDA_LDG_D (prt->x[p]);
            const float dy = lig_y2[l] - CUDA_LDG_D (prt->y[p]);
            const float dz = lig_z2[l] - CUDA_LDG_D (prt->z[p]);
            const float dst_pow2 = dx * dx + dy * dy + dz * dz;
            const float dst_pow4 = dst_pow2 * dst_pow2;
            const float dst = sqrtf (dst_pow2);
            //const float dst = 1.111f;


            /* hydrophobic potential */
            // slower on GK110
#if 0
            if (CUDA_LDG_D (prt->cdc[p]) == 1 && dst_pow2 <= 81.0f OROR1) {
                ehpc1 += CUDA_LDG_D (prt->hpp[p]) *
                    (1.0f - (3.5f / 81.0f * dst_pow2 -
                             4.5f / 81.0f / 81.0f * dst_pow4 +
                             2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
                             0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
            }
#endif
            // faster on GK110
#if 1
            const int mask_hpc = (CUDA_LDG_D (prt->cdc[p]) == 1 && dst_pow2 <= 81.0f);
            // constant folding optimization is working here
            const float temp_hpc = CUDA_LDG_D (prt->hpp[p]) *
                (1.0f - (3.5f / 81.0f * dst_pow2 -
                         4.5f / 81.0f / 81.0f * dst_pow4 +
                         2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
                         0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
            ehpc1 += temp_hpc * mask_hpc;
#endif





#if 1
            /* L-J potential */
            // p1a[MAXTP2][MAXTP1]
            // p2a[MAXTP2][MAXTP1]
            const float p1 = enepara_p1a[lig__t][prt__t] / (dst_pow4 * dst_pow4 * dst);
            const float p2 = enepara_p2a[lig__t][prt__t] / (dst_pow4 * dst_pow2);
            const float p4 = p1 * enepara_lj0 * (1.0f + enepara_lj1 * dst_pow2) + 1.0f;
            evdw += (p1 - p2) / p4;
#endif




            /* electrostatic potential */
            // slower on GK110
#if 0
            const float s1 = enepara_el1 * dst;
            float g1;
                g1 = 1.0f / s1;
            if (s1 < 1.0f)
                g1 = enepara_el0 + enepara_a1 * s1 * s1 + enepara_b1 * s1 * s1 * s1;
            eele += lig_c[l] * CUDA_LDG_D (prt->ele[p]) * g1;
#endif
            // faster on GK110
#if 1
            const float s1 = enepara_el1 * dst;
            const float g1_0 = enepara_el0 + enepara_a1 * s1 * s1 + enepara_b1 * s1 * s1 * s1;
            const float g1_1 = 1.0f / s1;
            const float g1 = (s1 < 1.0f) ? g1_0 : g1_1;
            eele += lig_c[l] * CUDA_LDG_D (prt->ele[p]) * g1;
#endif



            /* contact potential */
            // pmf0[MAXTP2][MAXTP1]
            // pmf1[MAXTP2][MAXTP1]
            // psp[MAXTP2][MAXPRO]

            const float dst_minus_pmf0 = dst - enepara_pmf0[lig__t][prt__t];

            epmf += enepara_pmf1[lig__t][prt__t] / (1.0f + expf ((-0.5f * dst + 6.0f) * dst_minus_pmf0));



            /* pocket-specific potential */
            // the senmatics do not match  with the original program:
            // if (found psp[][])
            //   accumulate to epsp
            // else
            //   do nothing
            // significantly better on GK110
#if 1
            if (CUDA_LDG_D (prt->c[p]) == 2 && dst_minus_pmf0 <= 0 OROR1) {
                const int i1 = CUDA_LDG_D (prt->seq3r[p]);
                epsp += CUDA_LDG_D (psp->psp[lig__t][i1]); // sparse matrix, indirect dereference
                // performance profiling:
                //epsp += float (lig__t + i1); // no difference

            }
#endif
            // signigicantly worse on GK110, because of "psp"
#if 0
            const int mask_psp = (CUDA_LDG_D (prt->c[p]) == 2 && dst_minus_pmf0 <= 0 OROR1);
            const int i1 = CUDA_LDG_D (prt->seq3r[p]);
            const float temp_psp = CUDA_LDG_D (psp->psp[lig__t][i1]); // sparse matrix, indirect dereference
            epsp += temp_psp * mask_psp;
#endif




            /* hydrogen bond potential */
            // hdb0[MAXTP2][MAXTP1]
            // hdb1[MAXTP2][MAXTP1]
            // a little bit faster on GK110
#if 1
            const float hdb0 = enepara_hdb0[lig__t][prt__t];
            if (hdb0 > 0.1f OROR1) {
                const float hdb1 = enepara_hdb1[lig__t][prt__t];
                const float hdb3 = (dst - hdb0) * hdb1;
                ehdb += hdb1 * expf (-0.5f * hdb3 * hdb3);
            }
#endif
            // a little bit slower on GK110
#if 0
            const float hdb0 = enepara_hdb0[lig__t][prt__t];
            const int mask_hdb = (hdb0 > 0.1f);
            const float hdb1 = enepara_hdb1[lig__t][prt__t];
            const float hdb3 = (dst - hdb0) * hdb1;
            const float temp_hdb = hdb1 * expf (-0.5f * hdb3 * hdb3);
            ehdb += temp_hdb * mask_hdb;
#endif


        } // prt loop
    } // if (l < lig_natom)



    /* hydrophobic restraits*/
    // hpl0[MAXTP2]
    // hpl1[MAXTP2]
    // hpl2[MAXTP2]

    //BlockReduceSum_2D_1_d <float> (bdy_prt, bdx_prt, ehpc1);
    BlockReduceSum_2D_1_d_2 (bdy_prt, bdx_prt, ehpc1);
    if (threadIdx.x < bdy_prt) {
        const int l = j + threadIdx.x;
        if (l < lig_natom) {
            const int lig__t = lig_t[l];
            const float hpc2 = (ehpc1 - enepara_hpl0[lig__t]) / enepara_hpl1[lig__t]; // div hpl is faster than mul hpl_inv
            ehpc += 0.5f * hpc2 * hpc2 - enepara_hpl2[lig__t];
        }
    }



} // lig loop


//BlockReduceSum_5_d <float, float, float, float, float> (evdw, eele, epmf, epsp, ehdb);
BlockReduceSum_5_d_2 (evdw, eele, epmf, epsp, ehdb);
WarpReduceSum_1_d_2 (ehpc);


