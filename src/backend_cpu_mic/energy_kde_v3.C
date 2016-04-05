
float ekde = 0.0f;


//#pragma loop count (32)
for (int l = 0; l < lig_natom; ++l) {	// lig loop, ~30
    float ekde1 = 0.0f;

    const int begin = kde_begin_idx[l];
    const int end = kde_end_idx[l];


    //#pragma unroll (2)
    //#pragma loop count (1024)


#ifndef NOVEC
    //#pragma simd reduction(+:ekde1)
    //#pragma simd private(dx, dy, dz, kde_dst_pow2)
#else
#pragma novector
#endif


    for (int k = begin; k < end; ++k) {	// kde loop, ~400
        const float dx = lig_x3[l] - kde->x[k];
        const float dy = lig_y3[l] - kde->y[k];
        const float dz = lig_z3[l] - kde->z[k];
        const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
        ekde1 += expf (enepara_kde2 * kde_dst_pow2);
    }

    float kde_sz = (float) (end - begin);
    ekde += ekde1 / kde_sz;

} // lig loop

e_s[6] = ekde * enepara_kde3_inv;
