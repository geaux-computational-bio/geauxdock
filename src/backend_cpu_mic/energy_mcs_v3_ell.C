
float elhm = 0.0f;


// icc outout: remark: unroll_and_jam pragma will be ignored due to (null)

//#pragma unroll_and_jam (4)



for (int m = 0; m < mcs_nrow; ++m) {
    float elhm1 = 0.0f;
    const int ncol = mcs_ncol[m];

    //#pragma loop count (32)

#ifndef NOVEC

    // default behavior of icc auto optimizer:
    // without SIMD pragma
    // icc outout: loop was not vectorized: vectorization possible but seems inefficient

#pragma simd reduction(+:elhm1)
#pragma vector aligned

#else

#pragma novector

#endif
    for (int i = 0; i < ncol; ++i) {
        const int l = mcs_ell->i[m][i];
        const float dx = lig_x2[l] - mcs_ell->x[m][i];
        const float dy = lig_y2[l] - mcs_ell->y[m][i];
        const float dz = lig_z2[l] - mcs_ell->z[m][i];
        elhm1 += dx * dx + dy * dy + dz * dz;
    }

    elhm += mcs_tcc[m] * sqrtf (elhm1 / (float) ncol);
}

e_s[7] = logf (elhm / mcs_nrow);
