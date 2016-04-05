
float ekde = 0.0f;
ty = threadIdx.x / bdx_kde;
tx = threadIdx.x % bdx_kde;
const float enepara_kde2_reg = enepara_kde2;
const int lig_natom_reg = lig_natom;


for (int j = 0; j < lig_natom_reg; j += bdy_kde) { // y loop, ~30
    float ekde1 = 0.0f;

    const int l = j + ty;
    //const int begin = kde_begin_idx[l];
    //const int end = kde_end_idx[l];

    if (l < lig_natom_reg) {

//#pragma unroll 2
        for (int k = kde_begin_idx[l] + tx; k < kde_end_idx[l]; k += bdx_kde) { // x loop, ~400
#if 0
            const float dx = lig_x3[l] - CUDA_LDG_D (kde->x[k]);
            const float dy = lig_y3[l] - CUDA_LDG_D (kde->y[k]);
            const float dz = lig_z3[l] - CUDA_LDG_D (kde->z[k]);
#endif
#if 1
            const float dx = lig_x3[l] - kde->x[k];
            const float dy = lig_y3[l] - kde->y[k];
            const float dz = lig_z3[l] - kde->z[k];
#endif
            const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
            ekde1 += expf (enepara_kde2_reg * kde_dst_pow2);
        }
    }

    BlockReduceSum_2D_1_d_2 (bdy_kde, bdx_kde, ekde1);
    if (threadIdx.x < bdy_kde) {
        const int l = j + threadIdx.x;
        const float kde_sz = (float) (kde_end_idx[l] - kde_begin_idx[l]);
        if (l < lig_natom_reg)
            ekde += ekde1 / kde_sz;
    }

} // lig loop

WarpReduceSum_1_d_2 (ekde);
if (threadIdx.x == 0)
    e_s[6] = ekde * enepara_kde3_inv;

