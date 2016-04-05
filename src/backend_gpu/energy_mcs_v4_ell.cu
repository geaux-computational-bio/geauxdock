
// sparse matrix in ELLPACK format
// no performance improvement



float elhm = 0.0f;
ty = threadIdx.x / bdx_mcs;
tx = threadIdx.x % bdx_mcs;


//#pragma unroll 4
for (int j = 0; j < mcs_nrow; j += bdy_mcs) { // y loop
    float elhm1 = 0.0f;
    int elhm2 = 0;

    const int m = j + ty;
    if (m < mcs_nrow) {
        for (int i = tx; i < mcs_ncol[m]; i += bdx_mcs) { // x loop
            const int l = mcs_ell->i[m][i];
            const float dx = lig_x2[l] - mcs_ell->x[m][i]; // do not use __LDG
            const float dy = lig_y2[l] - mcs_ell->y[m][i];
            const float dz = lig_z2[l] - mcs_ell->z[m][i];
            elhm1 += dx * dx + dy * dy + dz * dz;
            elhm2++; // never zero
        }
    }


    BlockReduceSum_2D_2_d_2 (bdy_mcs, bdx_mcs, elhm1, elhm2);

    if (threadIdx.x < bdy_mcs) {
        const int m = j + threadIdx.x;
        if (m < mcs_nrow)
            elhm += mcs_tcc[m] * sqrtf (elhm1 / (float) elhm2);
    }


} // lhm loop

WarpReduceSum_1_d_2 (elhm);
if (threadIdx.x == 0)
    e_s[7] = logf (elhm / mcs_nrow);


