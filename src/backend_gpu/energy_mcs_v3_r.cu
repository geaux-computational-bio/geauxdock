
// slow and WRONG!
// transpose te MCS matrix
// the x dimension is large enough
// the x dimension is tiled
// no simple method to reduce on Y dimension


float elhm = 0.0f;
ty = threadIdx.x / bdx_mcs;
tx = threadIdx.x % bdx_mcs;
__shared__ float elhm1_s[BD];
__shared__ int elhm2_s[BD];

for (int i = 0; i < mcs_nrow; i += bdx_mcs) { // x
    float elhm1 = 0.0f;
    int elhm2 = 0;

    const int m = i + tx;
    if (m < mcs_nrow) {
        for (int l = ty; l < lig_natom; l += bdy_mcs) { // y


            if (CUDA_LDG_D (mcs_r->x[l][m]) != MCS_INVALID_COORD OROR1) {
                const float dx = lig_x2[l] - CUDA_LDG_D (mcs_r->x[l][m]);
                const float dy = lig_y2[l] - CUDA_LDG_D (mcs_r->y[l][m]);
                const float dz = lig_z2[l] - CUDA_LDG_D (mcs_r->z[l][m]);
                elhm1 += dx * dx + dy * dy + dz * dz;
                elhm2++;
            }

        }
        elhm1_s[threadIdx.x] = elhm1;
        elhm2_s[threadIdx.x] = elhm2;
    }

    __syncthreads ();
    if (ty == 0) {
        for (int iter = 0; iter < bdy_mcs; ++iter) {
            elhm1 += elhm1_s[tx + bdx_mcs * iter]; // bank conflict
            elhm2 += elhm2_s[tx + bdx_mcs * iter];
        }
        elhm += mcs_tcc[m] * sqrtf (elhm1 / (float) elhm2);
    }

}

BlockReduceSum_2D_1_d_2 (bdy_mcs, bdx_mcs, elhm);
if (threadIdx.x == 0)
    e_s[7] = logf (elhm / mcs_nrow);





