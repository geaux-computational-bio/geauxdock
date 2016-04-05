
float ekde = 0.0f;
ty = threadIdx.x / bdx_kde;
tx = threadIdx.x % bdx_kde;


// lig loop, ~30
for (int j = 0; j < lig_natom; j += bdy_kde) { // y loop
    float ekde1 = 0.0f;
    int ekde2 = 0;

    const int l = j + ty;
    if (l < lig_natom) {

        /*
        // this slightly reduce performance
        // which implies compute bound
        float lx = lig_x3[l];
        float ly = lig_y3[l];
        float lz = lig_z3[l];
        float lt = lig_[l];
         */

        //#pragma unroll 2
        // kde loop, ~400
        for (int k = tx; k < kde_npoint; k += bdx_kde) { // x loop
#if 1
            if (lig_t3[l] == kde->t[k] OROR1) {

                const float dx = lig_x3[l] - kde->x[k];
                const float dy = lig_y3[l] - kde->y[k];
                const float dz = lig_z3[l] - kde->z[k];

                /*
                // this test code imporoves the performance by 3%
                // may not worth implementing memory tiling
                const float dx = lig_x3[l] - 1.1f; //kde->x[k];
                const float dy = lig_y3[l] - 1.2f; //kde->y[k];
                const float dz = lig_z3[l] - 1.3f; //kde->z[k];
                 */

                const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;

                // a little bit worse on GK110
#if 0
                // underflow protection
                // for IEEE Std 754-1985 double,
                // 709.8 < x implies exp(x) has overflowed
                // x < -708.4 implies exp(x) has underflowed
                const float t = enepara_kde2 * kde_dst_pow2;
                if (t > -30) {
                    //printf ("%f ", t);
                    ekde1 += expf (t);
                }
                ekde2++;
#endif
                // a little bit better on GK110
#if 1
                ekde1 += expf (enepara_kde2 * kde_dst_pow2);
                ekde2++;
#endif
            }
#endif

            // slower on GK110
#if 0

            const int mask_kde = (lig_t3[l] == kde->t[k]);
            const float dx = lig_x3[l] - kde->x[k];
            const float dy = lig_y3[l] - kde->y[k];
            const float dz = lig_z3[l] - kde->z[k];
            const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
            const float temp_kde = expf (enepara_kde2 * kde_dst_pow2);
            ekde1 += temp_kde * mask_kde;
            ekde2 += mask_kde;
#endif


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

