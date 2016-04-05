
for (int l = 0; l < lig_natom; ++l) {
    float ehpc1 = 0.0f;
    const int lig__t = lig_t[l];

    for (int p = 0; p < prt_npoint; ++p) {
        const int prt__t = prt->t[p];

        const float dx = lig_x1[l] - prt->x[p];
        const float dy = lig_y1[l] - prt->y[p];
        const float dz = lig_z1[l] - prt->z[p];
        const float dst_pow2 = dx * dx + dy * dy + dz * dz;
        const float dst_pow4 = dst_pow2 * dst_pow2;
        const float dst = sqrtf (dst_pow2);

        /* hydrophobic potential */
        if (prt->cdc[p] == 1 && dst_pow2 <= 81.0f OROR1) {
            ehpc1 += prt->hpp[p] *
                (1.0f - (3.5f / 81.0f * dst_pow2 -
                         4.5f / 81.0f / 81.0f * dst_pow4 +
                         2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
                         0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
        }

        /* L-J potential */
        const float p1 =
            enepara_p1a[lig__t][prt__t] / (dst_pow4 * dst_pow4 * dst);
        const float p2 = enepara_p2a[lig__t][prt__t] / (dst_pow4 * dst_pow2);
        const float p4 =
            p1 * enepara_lj0 * (1.0f + enepara_lj1 * dst_pow2) + 1.0f;
        evdw += (p1 - p2) / p4;

        /* electrostatic potential */
        const float s1 = enepara_el1 * dst;
        float g1;
        if (s1 < 1 OROR1)
            g1 = enepara_el0 + enepara_a1 * s1 * s1 + enepara_b1 * s1 * s1 * s1;
        else
            g1 = 1.0f / s1;
        eele += lig_c[l] * prt->ele[p] * g1;

        /* contact potential */
        const float dst_minus_pmf0 = dst - enepara_pmf0[lig__t][prt__t];
        epmf +=
            enepara_pmf1[lig__t][prt__t] /
            (1.0f + expf ((-0.5f * dst + 6.0f) * dst_minus_pmf0));

        /* pocket-specific potential */
        if (prt->c[p] == 2 && dst_minus_pmf0 <= 0 OROR1) {
            const int i1 = prt->seq3r[p];
            epsp += psp->psp[lig__t][i1];
        }

        /* hydrogen bond potential */
        const float hdb0 = enepara_hdb0[lig__t][prt__t];
        if (hdb0 > 0.1f OROR1) {
            const float hdb1 = enepara_hdb1[lig__t][prt__t];
            const float hdb3 = (dst - hdb0) * hdb1;
            ehdb += hdb1 * expf (-0.5f * hdb3 * hdb3);
        }

    }

    /* hydrophobic restraits */
    const float hpc2 = (ehpc1 - enepara_hpl0[lig__t]) / enepara_hpl1[lig__t];
    ehpc += 0.5f * hpc2 * hpc2 - enepara_hpl2[lig__t];
}



