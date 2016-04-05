//#include <cuda_profiler_api.h>


#define MYRAND (MyRand_d(curandstate_d))

//#define MYRAND 0.54f


// templating s2max cause a slow down from 337 to 331



// python3 -c "import math; print (-1.0 / math.sqrt (3.1415926535897932384626433 * 2.0))"
#define SQRT_2_PI_INV  -0.3989422804f

/*
#if __CUDA_ARCH__ >= 500
#define MY_KERNEL_MIN_BLOCKS 2
#elif __CUDA_ARCH__ >= 350
#define MY_KERNEL_MIN_BLOCKS 1
#endif
*/

//#define MY_KERNEL_MIN_BLOCKS 2




__global__ void
__launch_bounds__(BD, MC_BperMP)
MonteCarlo_d (Complex * __restrict__ complex,
    Record * __restrict__ record,
    const int s1, const int s2max, curandState *curandstate_d)
{

#if 1
    // moving "r" to shared memory reduce performance by 5%
    for (int r = complex->rep_begin + blockIdx.x; r <= complex->rep_end; r += gridDim.x) {

        // constant
        // pointers
        __shared__ ReplicaMC * __restrict__ rep;
        __shared__ Ligand * __restrict__ lig;
        __shared__ Protein * __restrict__ prt;
        __shared__ Psp * __restrict__ psp;
        __shared__ Kde * __restrict__ kde;
        __shared__ Mcs * __restrict__ mcs; // verified faster than accessing complex->mcs...
        //__shared__ Mcs_R * __restrict__ mcs_r;
        //__shared__ Mcs_ELL * __restrict__ mcs_ell;
        //__shared__ Mcs_CSR * __restrict__ mcs_csr;
        __shared__ EnePara * __restrict__ enepara;


        // temporary write-read
        __shared__ float move_scale[6]; // translation x y z, rotation x y z
        __shared__ float movematrix[6]; // translation x y z, rotation x y z
        __shared__ float rot[3][3]; // rotz roty rotx

        // might sorted
        //__shared__ float lig_x1[MAXLIG];
        //__shared__ float lig_y1[MAXLIG];
        //__shared__ float lig_z1[MAXLIG];

        // never sorted, for MCS compute
        __shared__ float lig_x2[MAXLIG];
        __shared__ float lig_y2[MAXLIG];
        __shared__ float lig_z2[MAXLIG];

        // sort by t, for KDE compute
        __shared__ float lig_x3[MAXLIG];
        __shared__ float lig_y3[MAXLIG];
        __shared__ float lig_z3[MAXLIG];



        // kde is always sorted by t
        //__shared__ int kde_begin_idx[MAXLIG];
        //__shared__ int kde_end_idx[MAXLIG];


        // constant
        // ligand
        __shared__ int   lig_t[MAXLIG];
        __shared__ float lig_c[MAXLIG];
        //__shared__ int   lig_t2[MAXLIG];
        __shared__ int   lig_t3[MAXLIG];
        //__shared__ float lig_c2[MAXLIG];
        __shared__ float lig_center[3];


        // constant
        // protein
        __shared__ float prt_pocket_center[3];


        // constant
        // mcs
        //__shared__ int mcs_ncol[MAX_MCS_ROW]; // row length in the sparse matrix representation
        __shared__ float mcs_tcc[MAX_MCS_ROW];
        //__shared__ float elhm1s[MAX_MCS_ROW];


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
        //__shared__ float enepara_a_para[MAXWEI]; // global memory is faster here
        //__shared__ float enepara_b_para[MAXWEI];
        //__shared__ float enepara_w[MAXWEI];
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
        //__shared__ float sqrt_2_pi_inv;
        __shared__ int lig_natom;
        __shared__ int prt_npoint;
        __shared__ int kde_npoint;
        __shared__ int mcs_nrow;
        //__shared__ int mcs_npoint;


        // temporary write-read vars
        __shared__ int is_accept_s;



        if (threadIdx.x == 0) {
            rep = &complex->replica[r];
            lig = &complex->lig[rep->idx_lig];
            prt = &complex->prt[rep->idx_prt];
            psp = &complex->psp;
            kde = &complex->kde;
            mcs = &complex->mcs[0];
            //mcs_r = &complex->mcs_r;
            //mcs_ell = &complex->mcs_ell;
            //mcs_csr = &complex->mcs_csr;
            enepara = &complex->enepara;
        }

        __syncthreads ();


        if (threadIdx.x < 6)
            move_scale[threadIdx.x] = complex->mcpara.move_scale[threadIdx.x];

        for (int l = 0; l < MAXTP2; ++l) {
            for (int p = threadIdx.x; p < MAXTP1; p += blockDim.x) {
                enepara_p1a[l][p] = enepara->p1a[l][p];
                enepara_p2a[l][p] = enepara->p2a[l][p];
                enepara_pmf0[l][p] = enepara->pmf0[l][p];
                enepara_pmf1[l][p] = enepara->pmf1[l][p];
                enepara_hdb0[l][p] = enepara->hdb0[l][p];
                enepara_hdb1[l][p] = enepara->hdb1[l][p];
            }
        }

        for (int l = threadIdx.x; l < MAXTP2; l += blockDim.x) {
            enepara_hpl0[l] = enepara->hpl0[l];
            enepara_hpl1[l] = enepara->hpl1[l];
            enepara_hpl2[l] = enepara->hpl2[l];
        }

        /*
           if (threadIdx.x < MAXWEI - 1) {
           enepara_a_para[threadIdx.x] = enepara->a_para[threadIdx.x];
           enepara_b_para[threadIdx.x] = enepara->b_para[threadIdx.x];
           enepara_w[threadIdx.x] = enepara->w[threadIdx.x];
           }
         */

        if (threadIdx.x == 0) {
            enepara_lj0 = enepara->lj0;
            enepara_lj1 = enepara->lj1;
            enepara_el1 = enepara->el1;
            enepara_el0 = enepara->el0;
            enepara_a1 = enepara->a1;
            enepara_b1 = enepara->b1;
            enepara_kde2 = enepara->kde2;
            enepara_kde3_inv = 1.0f / enepara->kde3;

            //sqrt_2pi_inv = -1.0f / sqrtf (2.0f * PI);
            lig_natom = lig->lig_natom;
            prt_npoint = prt->prt_npoint;
            kde_npoint = complex->size.kde_npoint;
            mcs_nrow = complex->size.mcs_nrow;
            //mcs_npoint = mcs_csr->npoint;

            is_accept_s = rep->is_accept;

            record[r - complex->rep_begin].next_entry = 0; // reset the record's entry point
        }

        __syncthreads ();

        for (int l = threadIdx.x; l < lig_natom; l += blockDim.x) {
            lig_t[l] = lig->t2[l]; // used in PRT
            lig_c[l] = lig->c2[l]; // used in PRT
            //lig_t2[l] = lig->t2[l]; // used in PRT
            //lig_c2[l] = lig->c2[l]; // used in PRT
            lig_t3[l] = lig->t3[l]; // used in KDE
            //kde_begin_idx[l] = lig->kde_begin_idx[l];
            //kde_end_idx[l] = lig->kde_end_idx[l];
        }

        if (threadIdx.x < 3) {
            lig_center[threadIdx.x] = lig->center[threadIdx.x];
            prt_pocket_center[threadIdx.x] = prt->pocket_center[threadIdx.x];
        }

        for (int m = threadIdx.x; m < mcs_nrow; m += blockDim.x) {
            //mcs_ncol[m] = mcs[m].ncol;
            mcs_tcc[m] = mcs[m].tcc;
        }



        // shapes of the 2D thread block
        // ty = threadIdx.x / bdx;
        // tx = threadIdx.x % bdx;

        // here __share__ is 10% slower
        int bdy_prt, bdx_prt, bdy_kde, bdx_kde, bdy_mcs, bdx_mcs;
        int ty, tx;

#if 0
__syncthreads ();
if (r == 0 && s2max == 1 && threadIdx.x == 0) {
    for (int i = 0; i < lig_natom; ++i) {
        printf ("lig_c[%2d] = %26.18f\n", i, lig_c[i]);
    }
}
#endif


#if 1
        bdy_prt = 8;
        bdx_prt = blockDim.x / bdy_prt;
        bdx_kde = 32;
        bdy_kde = blockDim.x / bdx_kde;
        //bdx_mcs = lig_natom >= 48 ? 64 : 32;
        bdx_mcs = 32;
        bdy_mcs = blockDim.x / bdx_mcs;
#endif


#if 0
        CalcThreadBlockShape (prt_npoint, lig_natom, bdx_prt, bdy_prt);
        CalcThreadBlockShape (kde_npoint, lig_natom, bdx_kde, bdy_kde);
        CalcThreadBlockShape (lig_natom, mcs_nrow, bdx_mcs, bdy_mcs);
#endif

#if 0
        if (r == 0 && threadIdx.x == 3) {
            printf ("data: %3d x %3d     TB: %3d x %3d\n", prt_npoint, lig_natom, bdx_prt, bdy_prt);
            printf ("data: %3d x %3d     TB: %3d x %3d\n", kde_npoint, lig_natom, bdx_kde, bdy_kde);
            printf ("data: %3d x %3d     TB: %3d x %3d\n", lig_natom, mcs_nrow, bdx_mcs, bdy_mcs);
        }
#endif










        for (int s2 = 0; s2 < s2max; ++s2) {


            /////////////////////////////////////////////////////////////////////////////
            // record old states
            // 1.0% time
            if (threadIdx.x == 0 && is_accept_s == 1) {
                rep->step = s1 + s2;

                const int rr = r - complex->rep_begin;
                const int next_entry = record[rr].next_entry;
                record[rr].replica[next_entry] = *rep;
                record[rr].next_entry = next_entry + 1;
            }






            /////////////////////////////////////////////////////////////////////////////
            // move

            if (threadIdx.x < 6) {

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

                movematrix[threadIdx.x] = move_scale[threadIdx.x] * moveamount + rep->movematrix[threadIdx.x];
            }



            __syncthreads ();


            {
                if (threadIdx.x == 0) {
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
            }

            __syncthreads ();

            // rotation, translation, coordinate system transformation
            for (int l = threadIdx.x; l < lig_natom; l += blockDim.x) {
                //const float x1 = lig->x1[l];
                //const float y1 = lig->y1[l];
                //const float z1 = lig->z1[l];
                const float x2 = lig->x2[l];
                const float y2 = lig->y2[l];
                const float z2 = lig->z2[l];
                const float x3 = lig->x3[l];
                const float y3 = lig->y3[l];
                const float z3 = lig->z3[l];
                //lig_x1[l] = rot[0][0] * x1 + rot[0][1] * y1 + rot[0][2] * z1 + movematrix[0] + lig_center[0];
                //lig_y1[l] = rot[1][0] * x1 + rot[1][1] * y1 + rot[1][2] * z1 + movematrix[1] + lig_center[1];
                //lig_z1[l] = rot[2][0] * x1 + rot[2][1] * y1 + rot[2][2] * z1 + movematrix[2] + lig_center[2];
                lig_x2[l] = rot[0][0] * x2 + rot[0][1] * y2 + rot[0][2] * z2 + movematrix[0] + lig_center[0];
                lig_y2[l] = rot[1][0] * x2 + rot[1][1] * y2 + rot[1][2] * z2 + movematrix[1] + lig_center[1];
                lig_z2[l] = rot[2][0] * x2 + rot[2][1] * y2 + rot[2][2] * z2 + movematrix[2] + lig_center[2];
                lig_x3[l] = rot[0][0] * x3 + rot[0][1] * y3 + rot[0][2] * z3 + movematrix[0] + lig_center[0];
                lig_y3[l] = rot[1][0] * x3 + rot[1][1] * y3 + rot[1][2] * z3 + movematrix[1] + lig_center[1];
                lig_z3[l] = rot[2][0] * x3 + rot[2][1] * y3 + rot[2][2] * z3 + movematrix[2] + lig_center[2];
            }

            __syncthreads ();



            //cudaProfilerInitialize ();
            //cudaProfilerStart ();

            /////////////////////////////////////////////////////////////////////////////
            // calcenergy

            float evdw = 0.0f;
            float eele = 0.0f;
            float epmf = 0.0f;
            float epsp = 0.0f;
            float ehdb = 0.0f;
            float ehpc = 0.0f;


#if CALC_PRT == 1
#include <energy_prt_v1.cu> // correct,
#endif

            __shared__ float e_s[MAXWEI];
            if (threadIdx.x == 0) {
                e_s[0] = evdw; // 0 - vdw
                e_s[1] = eele; // 1 - ele
                e_s[2] = epmf; // 2 - pmf (CP)
                e_s[3] = epsp; // 3 - psp (PS CP)
                e_s[4] = ehdb * SQRT_2_PI_INV; // 4 - hdb (HB)

                // constant folding optimization not work
                //e_s[4] = ehdb * (-1.0 / sqrtf (3.1415926535897932384626433 * 2.0));
                e_s[5] = ehpc; // 5 - hpc (HP)
            }



#if CALC_KDE == 1
#include <energy_kde_v1.cu> // correct, the paper, the best non-sparse implementation, CUDA thread 128*8
//#include <energy_kde_v2.cu> // correct, CUDA_LDG_D, performance varies on different code revisions
//#include <energy_kde_v3.cu> // correct, not computing ekde2, no faster
//#include <energy_kde_v4.cu> // correct, sparse format, the fastest, CUDA thread 32*X
#endif


#if CALC_MCS == 1
#include <energy_mcs_v1.cu> // correct, the paper, row: lig, colum: mcs.  the fastest
//#include <energy_mcs_v2.cu> // correct, not computing elhm2, no faster, bacause need extra data "mcs_ncol[]"
//#include <energy_mcs_v3_r.cu> // wrong, transpose x y, reduction is inefficient, also consume too much SMEM
//#include <energy_mcs_v4_ell.cu> // correct, sparse matrix, ELLPACK format, no faster than version 1
//#include <energy_mcs_v5_coo.cu> // correct, sparse matrix, COO format, slow due to atomicAdd
#endif


            //cudaProfilerStop ();


#if CALC_DST == 1
            // fully optimized
            {
                float dst;
                if (threadIdx.x < 3) {
                    dst = lig_center[threadIdx.x] +
                        movematrix[threadIdx.x] -
                        prt_pocket_center[threadIdx.x];
                    dst = dst * dst;
                }
                if (threadIdx.x == 0) {
                    dst += __shfl (dst, 1) + __shfl (dst, 2);
                    e_s[8] = sqrtf (dst);
                }
            }
#endif


            __syncthreads ();


            // normalization
            float e = 0.0f;
            if (threadIdx.x < MAXWEI - 1)
                e = e_s[threadIdx.x];
            if (threadIdx.x < 7)
                e /= lig_natom;
            if (threadIdx.x < MAXWEI - 1) {
                e = enepara->a_para[threadIdx.x] * e + enepara->b_para[threadIdx.x];
                e_s[threadIdx.x] = e;
                e *= enepara->w[threadIdx.x];
            }
            e += __shfl_xor (e, 8);
            e += __shfl_xor (e, 4);
            e += __shfl_xor (e, 2);
            e += __shfl_xor (e, 1);
            if (threadIdx.x == 0)
                e_s[MAXWEI - 1] = e;





            ////////////////////////////////////////////////////////////////////////
            // accept

            if (threadIdx.x == 0) {
                const float delta_energy = e - rep->energy[MAXWEI -1];
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
                if (threadIdx.x < MAXWEI)
                    rep->energy[threadIdx.x] = e_s[threadIdx.x];
                if (threadIdx.x < 6)
                    rep->movematrix[threadIdx.x] = movematrix[threadIdx.x];
            }


#if 0
            if (r == 0 && s2max == 1 && threadIdx.x == 0) {
                printf ("k1  ");
                for (int i = 0; i < MAXWEI; ++i)
                    printf ("%.8f ", e_s[i]);
                printf ("\n");


                printf ("ko    ");
                for (int i = 0; i < MAXWEI; ++i)
                    printf ("%.4f ", e_s[i]);
                printf ("\n");

            }
#endif


        } // s2 loop



        if (threadIdx.x == 0)
            rep->is_accept = is_accept_s;

    } // replica loop
#endif
}

