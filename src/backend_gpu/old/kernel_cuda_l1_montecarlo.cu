
#define MYRAND (MyRand_d(curandstate_d))

//#define MYRAND 0.54f


// templating s2max cause a slow down from 337 to 331



// python: sqrt_2_pi_inv = -1.0 / math.sqrt (3.1415926535897932384626433 * 2.0)
#define SQRT_2_PI_INV  -0.398942f



__global__ void
//__launch_bounds__(512, 1)
MonteCarlo_d (Complex * __restrict__ complex,
    Record * __restrict__ record,
    const int s1, const int s2max, curandState *curandstate_d)
{
    const int bidx = blockDim.x * threadIdx.y + threadIdx.x;

#if 1
    for (int r = complex->rep_begin + blockIdx.x; r <= complex->rep_end; r += gridDim.x) {

        // constant
        // pointers
        __shared__ ReplicaMC * __restrict__ rep;
        __shared__ Ligand * __restrict__ lig;
        __shared__ Protein * __restrict__ prt;
        __shared__ Psp * __restrict__ psp;
        __shared__ Kde * __restrict__ kde;
        __shared__ Mcs * __restrict__ mcs;
        __shared__ EnePara * __restrict__ enepara;


        // temporary write-read
        __shared__ float move_scale[6]; // translation x y z, rotation x y z
        __shared__ float movematrix[6]; // translation x y z, rotation x y z
        __shared__ float rot[3][3]; // rotz roty rotx
        __shared__ float lig_x1[MAXLIG];
        __shared__ float lig_y1[MAXLIG];
        __shared__ float lig_z1[MAXLIG];
        __shared__ float lig_x2[MAXLIG];
        __shared__ float lig_y2[MAXLIG];
        __shared__ float lig_z2[MAXLIG];


        // constant
        // ligand
        __shared__ int lig_t[MAXLIG];
        __shared__ float lig_c[MAXLIG];
        __shared__ float lig_center[3];


        // constant
        // protein
        __shared__ float prt_pocket_center[3];


        // constant
        // mcs
        __shared__ float mcs_tcc[MAX_MCS_ROW];


        // constant
        // enepara
        __shared__ float enepara_p1a[MAXTP2][MAXTP1];
        __shared__ float enepara_p2a[MAXTP2][MAXTP1];
        __shared__ float enepara_pmf0[MAXTP2][MAXTP1];
        __shared__ float enepara_pmf1[MAXTP2][MAXTP1];
        __shared__ float enepara_hdb0[MAXTP2][MAXTP1];
        __shared__ float enepara_hdb1[MAXTP2][MAXTP1];
        __shared__ float enepara_hpl0[MAXTP2];
        __shared__ float enepara_hpl1[MAXTP2];
        __shared__ float enepara_hpl2[MAXTP2];
        __shared__ float enepara_a_para[MAXWEI];
        __shared__ float enepara_b_para[MAXWEI];
        __shared__ float enepara_w[MAXWEI];
        __shared__ float enepara_lj0;
        __shared__ float enepara_lj1;
        __shared__ float enepara_el1;
        __shared__ float enepara_el0;
        __shared__ float enepara_a1;
        __shared__ float enepara_b1;
        __shared__ float enepara_kde2;
        __shared__ float enepara_kde3_inv;


        // constant
        // scalars
        __shared__ float sqrt_2_pi_inv;
        __shared__ int lig_natom, prt_npoint, kde_npoint, mcs_nrow;


        // temporary write-read vars
        __shared__ int is_accept_s;



        if (bidx == 0) {
            rep = &complex->replica[r];
            lig = &complex->lig[rep->idx_lig];
            prt = &complex->prt[rep->idx_prt];
            psp = &complex->psp;
            kde = &complex->kde;
            mcs = &complex->mcs[0];
            enepara = &complex->enepara;
        }

        __syncthreads ();


        if (bidx < 6)
            move_scale[bidx] = complex->mcpara.move_scale[bidx];

        for (int l = threadIdx.y; l < MAXTP2; l += blockDim.y) {
            for (int p = threadIdx.x; p < MAXTP1; p += blockDim.x) {
                enepara_p1a[l][p] = enepara->p1a[l][p];
                enepara_p2a[l][p] = enepara->p2a[l][p];
                enepara_pmf0[l][p] = enepara->pmf0[l][p];
                enepara_pmf1[l][p] = enepara->pmf1[l][p];
                enepara_hdb0[l][p] = enepara->hdb0[l][p];
                enepara_hdb1[l][p] = enepara->hdb1[l][p];
            }
        }

        for (int l = bidx; l < MAXTP2; l += blockDim.x * blockDim.y) {
            enepara_hpl0[l] = enepara->hpl0[l];
            enepara_hpl1[l] = enepara->hpl1[l];
            enepara_hpl2[l] = enepara->hpl2[l];
        }

        if (bidx < MAXWEI - 1) {
            enepara_a_para[bidx] = enepara->a_para[bidx];
            enepara_b_para[bidx] = enepara->b_para[bidx];
            enepara_w[bidx] = enepara->w[bidx];
        }

        if (bidx == 0) {
            enepara_lj0 = enepara->lj0;
            enepara_lj1 = enepara->lj1;
            enepara_el1 = enepara->el1;
            enepara_el0 = enepara->el0;
            enepara_a1 = enepara->a1;
            enepara_b1 = enepara->b1;
            enepara_kde2 = enepara->kde2;
            enepara_kde3_inv = 1.0f / enepara->kde3;

            sqrt_2_pi_inv = -1.0f / sqrtf (2.0f * PI);
            lig_natom = lig->lig_natom;
            prt_npoint = prt->prt_npoint;
            kde_npoint = complex->size.kde_npoint;
            mcs_nrow = complex->size.mcs_nrow;
            is_accept_s = rep->is_accept;

            record[r - complex->rep_begin].next_entry = 0; // reset the record's entry point
        }

        __syncthreads ();

        for (int l = bidx; l < lig_natom; l += blockDim.x * blockDim.y) {
            lig_t[l] = lig->t[l];
            lig_c[l] = lig->c[l];
        }

        if (bidx < 3) {
            lig_center[bidx] = lig->center[bidx];
            prt_pocket_center[bidx] = prt->pocket_center[bidx];
        }

        for (int m = bidx; m < mcs_nrow; m += blockDim.x * blockDim.y) {
            mcs_tcc[m] = mcs[m].tcc;
        }

        __syncthreads ();







        for (int s2 = 0; s2 < s2max; ++s2) {


            /////////////////////////////////////////////////////////////////////////////
            // record old states
            // 1.0% time
            if (bidx == 0 && is_accept_s == 1) {
                rep->step = s1 + s2;

                const int rr = r - complex->rep_begin;
                const int next_entry = record[rr].next_entry;
                record[rr].replica[next_entry] = *rep;
                record[rr].next_entry = next_entry + 1;
            }






            /////////////////////////////////////////////////////////////////////////////
            // move

            if (bidx < 6) {

#if IS_AWAY == 0
                const float fixed_var = 0.0f;
#elif IS_AWAY == 1
                const float fixed_var = 44.5f;
#endif

#if 1
                float moveamount;
                if (s2max != 1)
                    moveamount = MYRAND;
                else
                    moveamount = fixed_var;
#endif

#if 0
                float moveamount = (s2max != 1) ? MYRAND : fixed_var;
#endif




                movematrix[bidx] = move_scale[bidx] * moveamount + rep->movematrix[bidx];
            }

            __syncthreads ();


            if (bidx == 0) {
                // http://en.wikipedia.org/wiki/Euler_angles
                // http://upload.wikimedia.org/math/e/9/c/e9cf817bce9c1780216921cd93233459.png
                // http://upload.wikimedia.org/math/f/4/e/f4e55dc2c9581007648967d29b15121e.png
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
            }

            __syncthreads ();

            // rotation, translation, coordinate system transformation
            for (int l = bidx; l < lig_natom; l += blockDim.x * blockDim.y) {
                const float x1 = lig->x[l];
                const float y1 = lig->y[l];
                const float z1 = lig->z[l];
                lig_x1[l] = rot[0][0] * x1 + rot[0][1] * y1 + rot[0][2] * z1 + movematrix[0] + lig_center[0];
                lig_y1[l] = rot[1][0] * x1 + rot[1][1] * y1 + rot[1][2] * z1 + movematrix[1] + lig_center[1];
                lig_z1[l] = rot[2][0] * x1 + rot[2][1] * y1 + rot[2][2] * z1 + movematrix[2] + lig_center[2];

                const float x2 = lig->x2[l];
                const float y2 = lig->y2[l];
                const float z2 = lig->z2[l];
                lig_x2[l] = rot[0][0] * x2 + rot[0][1] * y2 + rot[0][2] * z2 + movematrix[0] + lig_center[0];
                lig_y2[l] = rot[1][0] * x2 + rot[1][1] * y2 + rot[1][2] * z2 + movematrix[1] + lig_center[1];
                lig_z2[l] = rot[2][0] * x2 + rot[2][1] * y2 + rot[2][2] * z2 + movematrix[2] + lig_center[2];
            }

            __syncthreads ();





            /////////////////////////////////////////////////////////////////////////////
            // calcenergy

            float evdw = 0.0f;
            float eele = 0.0f;
            float epmf = 0.0f;
            float epsp = 0.0f;
            float ehdb = 0.0f;
            float ehpc = 0.0f;


#if CALC_PRT == 1

            // lig loop, ~30
            //#pragma unroll 8
            for (int j = 0; j < lig_natom; j += blockDim.y) {
                float ehpc1 = 0.0f;

                {
                    const int l = j + threadIdx.y;
                    if (l < lig_natom) {
                        const int lig__t = lig_t[l];

                        //#pragma unroll 4
                        // prt loop, ~300
                        for (int p = threadIdx.x; p < prt_npoint; p += blockDim.x) {
                            const int prt__t = CUDA_LDG_D (prt->t[p]);

                            const float dx = lig_x1[l] - CUDA_LDG_D (prt->x[p]);
                            const float dy = lig_y1[l] - CUDA_LDG_D (prt->y[p]);
                            const float dz = lig_z1[l] - CUDA_LDG_D (prt->z[p]);
                            const float dst_pow2 = dx * dx + dy * dy + dz * dz;
                            const float dst_pow4 = dst_pow2 * dst_pow2;
                            const float dst = sqrtf (dst_pow2);


                            /* hydrophobic potential */
                            {
                                // worse
#if 0
                                if (CUDA_LDG_D (prt->cdc[p]) == 1 && dst_pow2 <= 81.0f OROR1) {
                                    ehpc1 += CUDA_LDG_D (prt->hpp[p]) *
                                        (1.0f - (3.5f / 81.0f * dst_pow2 -
                                                 4.5f / 81.0f / 81.0f * dst_pow4 +
                                                 2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
                                                 0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
                                }
#endif
                                // better
#if 1
                                const int mask_hpc = (CUDA_LDG_D (prt->cdc[p]) == 1 && dst_pow2 <= 81.0f OROR1);
                                const float temp_hpc = CUDA_LDG_D (prt->hpp[p]) *
                                    (1.0f - (3.5f / 81.0f * dst_pow2 -
                                             4.5f / 81.0f / 81.0f * dst_pow4 +
                                             2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
                                             0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
                                ehpc1 += temp_hpc * mask_hpc;
#endif
                            }





#if 1
                            /* L-J potential */
                            // p1a[MAXTP2][MAXTP1]
                            // p2a[MAXTP2][MAXTP1]
                            {
                                const float p1 = enepara_p1a[lig__t][prt__t] / (dst_pow4 * dst_pow4 * dst);
                                const float p2 = enepara_p2a[lig__t][prt__t] / (dst_pow4 * dst_pow2);
                                const float p4 = p1 * enepara_lj0 * (1.0f + enepara_lj1 * dst_pow2) + 1.0f;
                                evdw += (p1 - p2) / p4;
                            }
#endif




                            /* electrostatic potential */
                            {
                                // better
#if 1
                                const float s1 = enepara_el1 * dst;
                                float g1;
                                if (s1 < 1.0f OROR1)
                                    g1 = enepara_el0 + enepara_a1 * s1 * s1 + enepara_b1 * s1 * s1 * s1;
                                else
                                    g1 = 1.0f / s1;
                                eele += lig_c[l] * CUDA_LDG_D (prt->ele[p]) * g1;
#endif
                                // worse
#if 0
                                const float s1 = enepara_el1 * dst;
                                const float g1_0 = enepara_el0 + enepara_a1 * s1 * s1 + enepara_b1 * s1 * s1 * s1;
                                const float g1_1 = 1.0f / s1;
                                const float g1 = (s1 < 1) ? g1_0 : g1_1;
                                eele += lig_c[l] * CUDA_LDG_D (prt->ele[p]) * g1;
#endif
                            }



                            /* contact potential */
                            // pmf0[MAXTP2][MAXTP1]
                            // pmf1[MAXTP2][MAXTP1]
                            // psp[MAXTP2][MAXPRO]

                            const float dst_minus_pmf0 = dst - enepara_pmf0[lig__t][prt__t];

                            epmf += enepara_pmf1[lig__t][prt__t] / (1.0f + expf ((-0.5f * dst + 6.0f) * dst_minus_pmf0));



                            /* pocket-specific potential */
                            // the senmatics do not match with the original program:
                            // if (found psp[][])
                            //   accumulate to epsp
                            // else
                            //   do nothing
                            {
                                // better
#if 1
                                if (CUDA_LDG_D (prt->c[p]) == 2 && dst_minus_pmf0 <= 0 OROR1) {
                                    const int i1 = CUDA_LDG_D (prt->seq3r[p]);
                                    epsp += CUDA_LDG_D (psp->psp[lig__t][i1]); // sparse matrix, indirect dereference
                                    // performance profiling:
                                    //epsp += float (lig__t + i1); // improve from 336 to 352, not worth doing

                                }
#endif
                                // worse
#if 0
                                const int mask_psp = (CUDA_LDG_D (prt->c[p]) == 2 && dst_minus_pmf0 <= 0 OROR1);
                                const int i1 = CUDA_LDG_D (prt->seq3r[p]);
                                const float temp_psp = CUDA_LDG_D (psp->psp[lig__t][i1]); // sparse matrix, indirect dereference
                                epsp += temp_psp * mask_psp;
#endif
                            }




                            /* hydrogen bond potential */
                            // hdb0[MAXTP2][MAXTP1]
                            // hdb1[MAXTP2][MAXTP1]
                            {
                                // better
#if 1
                                const float hdb0 = enepara_hdb0[lig__t][prt__t];
                                if (hdb0 > 0.1f OROR1) {
                                    const float hdb1 = enepara_hdb1[lig__t][prt__t];
                                    const float hdb3 = (dst - hdb0) * hdb1;
                                    ehdb += hdb1 * expf (-0.5f * hdb3 * hdb3);
                                }
#endif
                                // worse
#if 0
                                const float hdb0 = enepara_hdb0[lig__t][prt__t];
                                const int mask_hdb = (hdb0 > 0.1f OROR1);
                                const float hdb1 = enepara_hdb1[lig__t][prt__t];
                                const float hdb3 = (dst - hdb0) * hdb1;
                                const float temp_hdb = hdb1 * expf (-0.5f * hdb3 * hdb3);
                                ehdb += temp_hdb * mask_hdb;
#endif
                            }


                        } // prt loop
                    } // if (l < lig_natom)
                } // end of scope "const int l = j + threadIdx.y"



                /* hydrophobic restraits*/
                // hpl0[MAXTP2]
                // hpl1[MAXTP2]
                // hpl2[MAXTP2]

                BlockReduceSum_2D_1_d <float> (bidx, ehpc1);
                if (threadIdx.y == 0 && threadIdx.x < blockDim.y) {
                    const int l = j + threadIdx.x;
                    if (l < lig_natom) {
                        const int lig__t = lig_t[l];
                        const float hpc2 = (ehpc1 - enepara_hpl0[lig__t]) / enepara_hpl1[lig__t]; // div hpl is faster than mul hpl_inv
                        ehpc += 0.5f * hpc2 * hpc2 - enepara_hpl2[lig__t];
                    }
                }



            } // lig loop


            BlockReduceSum_5_d <float, float, float, float, float> (bidx, evdw, eele, epmf, epsp, ehdb);
            WarpReduceSum_1_d <float> (ehpc);

#endif

            __shared__ float e[MAXWEI];
            if (bidx == 0) {
                e[0] = evdw; // 0 - vdw 
                e[1] = eele; // 1 - ele
                e[2] = epmf; // 2 - pmf (CP)
                e[3] = epsp; // 3 - psp (PS CP)
                e[4] = ehdb * SQRT_2_PI_INV; // 4 - hdb (HB)
                //e[4] = ehdb * sqrt_2_pi_inv; // 4 - hdb (HB)
                e[5] = ehpc; // 5 - hpc (HP)
            }








#if CALC_KDE == 1
            /* kde potential */
            // fully optimized

            float ekde = 0.0f;

            // lig loop, ~30
            for (int j = 0; j < lig_natom; j += blockDim.y) {
                float ekde1 = 0.0f;
                int ekde2 = 0;

                {
                    const int l = j + threadIdx.y;
                    if (l < lig_natom) {

                        //#pragma unroll 2
                        // kde loop, ~400
                        for (int k = threadIdx.x; k < kde_npoint; k += blockDim.x) {
#if 1
                            if (lig_t[l] == kde->t[k] OROR1) {
                                const float dx = lig_x1[l] - kde->x[k];
                                const float dy = lig_y1[l] - kde->y[k];
                                const float dz = lig_z1[l] - kde->z[k];
                                const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
                                ekde1 += expf (enepara_kde2 * kde_dst_pow2);
                                ekde2++;
#endif

#if 0
                                const int mask_kde = (lig_t[l] == kde->t[k] OROR1);
                                const float dx = lig_x1[l] - kde->x[k];
                                const float dy = lig_y1[l] - kde->y[k];
                                const float dz = lig_z1[l] - kde->z[k];
                                const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
                                const float temp_kde = expf (enepara_kde2 * kde_dst_pow2);
                                ekde1 += temp_kde * mask_kde;
                                ekde2 += mask_kde;
#endif


                            }
                        } // kde loop
                    } // if (l < lig_natom)
                }

                BlockReduceSum_2D_2_d <float, int> (bidx, ekde1, ekde2);
                if (threadIdx.y == 0 && threadIdx.x < blockDim.y) {
                    const int l = j + threadIdx.x;
                    if (l < lig_natom && ekde2 != 0)
                        ekde += ekde1 / (float) ekde2;
                }

            } // lig loop

            WarpReduceSum_1_d <float> (ekde);
            if (bidx == 0)
                e[6] = ekde * enepara_kde3_inv;

#endif


#if CALC_MCS == 1
            /* position restraints */
            // fully optimized

            float elhm = 0.0f;

            // lhm loop, ~11
            // unrolling outer loop does not help
            //#pragma unroll 4
            for (int j = 0; j < mcs_nrow; j += blockDim.y) {
                float elhm1 = 0.0f;
                int elhm2 = 0;

                {
                    const int m = j + threadIdx.y;
                    if (m < mcs_nrow) {
                        // lig loop, ~30

                        for (int l = threadIdx.x; l < lig_natom; l += blockDim.x) {
                            if (CUDA_LDG_D (mcs[m].x[l + 1]) != MCS_INVALID_COORD OROR1) {
                                const float dx = lig_x2[l] - CUDA_LDG_D (mcs[m].x[l + 1]);
                                const float dy = lig_y2[l] - CUDA_LDG_D (mcs[m].y[l + 1]);
                                const float dz = lig_z2[l] - CUDA_LDG_D (mcs[m].z[l + 1]);
                                elhm1 += dx * dx + dy * dy + dz * dz;
                                elhm2++;
                            }
                        } // lig loop
                    } // if (m < mcs_nrow)
                }

                BlockReduceSum_2D_2_d <float, int> (bidx, elhm1, elhm2);
                if (threadIdx.y == 0 && threadIdx.x < blockDim.y) {
                    const int m = j + threadIdx.x;
                    if (m < mcs_nrow && elhm2 != 0)
                        elhm += mcs_tcc[m] * sqrtf (elhm1 / (float) elhm2);
                }

            } // lhm loop

            WarpReduceSum_1_d <float> (elhm);
            if (bidx == 0)
                e[7] = logf (elhm / mcs_nrow);

#endif


#if CALC_DST == 1
            // fully optimized

            {
                float dst;
                if (bidx < 3) {
                    dst = lig_center[bidx] + movematrix[bidx] - prt_pocket_center[bidx];
                    dst = dst * dst;
                }
                if (bidx == 0) {
                    dst += __shfl (dst, 1) + __shfl (dst, 2);
                    e[8] = sqrtf (dst);
                }
            }
#endif


            __syncthreads ();

            // normalization
            if (bidx < 7)
                e[bidx] = e[bidx] / lig_natom;
            if (bidx < MAXWEI - 1)
                e[bidx] = enepara_a_para[bidx] * e[bidx] + enepara_b_para[bidx];
            float etotal = 0.0f;
            if (bidx < MAXWEI - 1)
                etotal = enepara->w[bidx] * e[bidx]; // enepara->w is faster than enepara_w
            WarpReduceSum_1_d <float> (etotal);
            if (bidx == 0)
                e[MAXWEI - 1] = etotal;



            __syncthreads ();


            ////////////////////////////////////////////////////////////////////////
            // accept

            if (bidx == 0) {
                const float delta_energy = e[MAXWEI - 1] - rep->energy[MAXWEI -1];
                const float beta = complex->temp[rep->idx_tmp].minus_beta;
#if 1
                float rand;
                if (s2max != 1)
                    rand = MYRAND;
                else
                    rand = 0.0f; // force to accept if s2max == 1
#endif
#if 0
                float rand = (s2max != 1) ? MYRAND : 0.0f;
#endif

                is_accept_s = (rand < expf (delta_energy * beta));  // mybeta < 0
            }
            __syncthreads ();
            if (is_accept_s == 1) {
                if (bidx < MAXWEI)
                    rep->energy[bidx] = e[bidx];
                if (bidx < 6)
                    rep->movematrix[bidx] = movematrix[bidx];
            }


        } // s2 loop



        if (bidx == 0)
            rep->is_accept = is_accept_s;

    } // replica loop
#endif
}

