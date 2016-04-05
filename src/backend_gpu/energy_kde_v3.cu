

// computing kde_sz in a different way



float ekde = 0.0f;
ty = threadIdx.x / bdx_kde;
tx = threadIdx.x % bdx_kde;


// lig loop, ~30
for (int j = 0; j < lig_natom; j += bdy_kde) { // y loop
    float ekde1 = 0.0f;
    const int l = j + ty;

    if (l < lig_natom) {

        // kde loop, ~400
        for (int k = tx; k < kde_npoint; k += bdx_kde) { // x loop
            if (lig_t3[l] == CUDA_LDG_D (kde->t[k]) OROR1) {
                const float dx = lig_x3[l] - CUDA_LDG_D (kde->x[k]);
                const float dy = lig_y3[l] - CUDA_LDG_D (kde->y[k]);
                const float dz = lig_z3[l] - CUDA_LDG_D (kde->z[k]);
                const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
                ekde1 += expf (enepara_kde2 * kde_dst_pow2);
            }

        } // kde loop
    } // if (l < lig_natom)

    BlockReduceSum_2D_1_d_2 (bdy_kde, bdx_kde, ekde1);
    if (threadIdx.x < bdy_kde) {
        const int l = j + threadIdx.x;
        const float kde_sz = (float) (kde_end_idx[l] - kde_begin_idx[l]);
        if (l < lig_natom)
            ekde += ekde1 / kde_sz;
    }

} // lig loop

WarpReduceSum_1_d_2 (ekde);
if (threadIdx.x == 0)
    e_s[6] = ekde * enepara_kde3_inv;

