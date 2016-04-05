
float ekde = 0.0f;
ty = threadIdx.x / bdx_kde;
tx = threadIdx.x % bdx_kde;


// lig loop, ~30
for (int j = 0; j < lig_natom; j += bdy_kde) { // y loop
    float ekde1 = 0.0f;
    int ekde2 = 0;

    const int l = j + ty;
    if (l < lig_natom) {

        // kde loop, ~400
        for (int k = tx; k < kde_npoint; k += bdx_kde) { // x loop

            if (lig_t3[l] == CUDA_LDG_D (kde->t[k]) OROR1) {
#if 1
                const float dx = lig_x3[l] - CUDA_LDG_D (kde->x[k]);
                const float dy = lig_y3[l] - CUDA_LDG_D (kde->y[k]);
                const float dz = lig_z3[l] - CUDA_LDG_D (kde->z[k]);
#endif

#if 0
                // this test code imporoves the performance 202 -> 216
                // may not worth implementing memory tiling for KDE data
                const float dx = lig_x3[l] - k;
                const float dy = lig_y3[l] - k;
                const float dz = lig_z3[l] - k;
#endif

                const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
                ekde1 += expf (enepara_kde2 * kde_dst_pow2);
                ekde2++;
            }

        } // kde loop
    } // if (l < lig_natom)

    BlockReduceSum_2D_2_d_2 (bdy_kde, bdx_kde, ekde1, ekde2);
    if (threadIdx.x < bdy_kde) {
        const int l = j + threadIdx.x;
        if (l < lig_natom && ekde2 != 0)
            ekde += ekde1 / (float) ekde2;
    }

} // lig loop

WarpReduceSum_1_d_2 (ekde);
if (threadIdx.x == 0)
    e_s[6] = ekde * enepara_kde3_inv;

