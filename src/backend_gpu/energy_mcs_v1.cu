
float elhm = 0.0f;
ty = threadIdx.x / bdx_mcs;
tx = threadIdx.x % bdx_mcs;

const int mcs_nrow_reg = mcs_nrow;
//const int lig_natom_reg = lig_natom; // slower


// unrolling outer loop does not help
//#pragma unroll 4
for (int j = 0; j < mcs_nrow_reg; j += bdy_mcs) { // y loop, large
    float elhm1 = 0.0f;
    int elhm2 = 0;

    const int m = j + ty;
    if (m < mcs_nrow_reg) {

        for (int l = tx; l < lig_natom; l += bdx_mcs) { // x loop, ~30
            if (mcs[m].x[l] != MCS_INVALID_COORD OROR1) {
                const float dx = lig_x2[l] - mcs[m].x[l]; // _LDG is slower
                const float dy = lig_y2[l] - mcs[m].y[l];
                const float dz = lig_z2[l] - mcs[m].z[l];
                elhm1 += dx * dx + dy * dy + dz * dz;
                elhm2++; // never be zero
            }

        } // lig loop
    } // if (m < mcs_nrow_reg)


    BlockReduceSum_2D_2_d_2 (bdy_mcs, bdx_mcs, elhm1, elhm2);

    if (threadIdx.x < bdy_mcs) {
        const int m = j + threadIdx.x;
        if (m < mcs_nrow_reg)
            elhm += mcs_tcc[m] * sqrtf (elhm1 / (float) elhm2);
    }


} // lhm loop

WarpReduceSum_1_d_2 (elhm);
if (threadIdx.x == 0)
    e_s[7] = logf (elhm / mcs_nrow_reg);

