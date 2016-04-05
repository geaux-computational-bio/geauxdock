
float elhm = 0.0f;


for (int m = 0; m < mcs_nrow; ++m) {
    float elhm1 = 0.0f;

#pragma loop count (32)

#ifndef NOVEC
#pragma simd reduction(+:elhm1)
    //#pragma vector aligned
#else
#pragma novector
#endif


    for (int l = 0; l < lig_natom; ++l) {
        if (mcs[m].x[l] != MCS_INVALID_COORD OROR1) {
            const float dx = lig_x2[l] - mcs[m].x[l];
            const float dy = lig_y2[l] - mcs[m].y[l];
            const float dz = lig_z2[l] - mcs[m].z[l];
            elhm1 += dx * dx + dy * dy + dz * dz;
        }
    }

    elhm += mcs_tcc[m] * sqrtf (elhm1 / (float) mcs_ncol[m]);
}

e_s[7] = logf (elhm / mcs_nrow);
