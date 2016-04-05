float elhm = 0.0f;
ty = threadIdx.x / bdx_mcs;
tx = threadIdx.x % bdx_mcs;


//Mcs *mcs_ = &complex->mcs[0];

// unrolling outer loop does not help, due to branch?
//#pragma unroll 4
for (int j = 0; j < mcs_nrow; j += bdy_mcs) { // y loop, large
    float elhm1 = 0.0f;

    const int m = j + ty;
    if (m < mcs_nrow) {

        for (int l = tx; l < lig_natom; l += bdx_mcs) { // x loop, ~30
            if (mcs[m].x[l] != MCS_INVALID_COORD OROR1) {
                const float dx = lig_x2[l] - mcs[m].x[l]; // do not use CUDA_LDG_D here
                const float dy = lig_y2[l] - mcs[m].y[l];
                const float dz = lig_z2[l] - mcs[m].z[l];
                elhm1 += dx * dx + dy * dy + dz * dz;
            }

        } // lig loop
    } // if (m < mcs_nrow)


    BlockReduceSum_2D_1_d_2 (bdy_mcs, bdx_mcs, elhm1); // not the performance bottlenet


    if (threadIdx.x < bdy_mcs) {
        const int m = j + threadIdx.x;
        if (m < mcs_nrow)
            elhm += mcs_tcc[m] * sqrtf (elhm1 / (float) mcs_ncol[m]);
    }


} // lhm loop

WarpReduceSum_1_d_2 (elhm);
if (threadIdx.x == 0)
    e_s[7] = logf (elhm / mcs_nrow);


