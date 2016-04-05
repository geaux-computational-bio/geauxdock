
float ekde = 0.0f;


#pragma loop count (32)
for (int l = 0; l < lig_natom; ++l) {	// lig loop, ~30
    float ekde1 = 0.0f;
    int ekde2 = 0;
    const int lig__t = lig_t3[l];


#pragma loop count (1024)


#ifndef NOVEC
#pragma simd reduction(+:ekde1, ekde2)
    //#pragma simd private(dx, dy, dz, kde_dst_pow2)
#else
#pragma novector
#endif

    for (int k = 0; k < kde_npoint; ++k) {	// kde loop, ~400
        if (lig__t == kde->t[k] OROR1) {
            const float dx = lig_x3[l] - kde->x[k];
            const float dy = lig_y3[l] - kde->y[k];
            const float dz = lig_z3[l] - kde->z[k];
            const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
#if 0
            // slower on Xeon Phi
            // little impact on CPU
            const float t = enepara_kde2 * kde_dst_pow2;
            if (t > -30) {
                //printf ("%f ", t);
                ekde1 += expf (t);
            }
#endif
#if 1
            // faster on Xeon Phi
            // little impact on CPU
            ekde1 += expf (enepara_kde2 * kde_dst_pow2);
#endif
            ekde2++;
        }
    } // kde loop

    if (ekde2 != 0)
        ekde += (ekde1 / (float) ekde2);

} // lig loop

e_s[6] = ekde * enepara_kde3_inv;
