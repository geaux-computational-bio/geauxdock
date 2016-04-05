float ekde = 0.0f;

float ekde1;
int lig__tk;

int begin;
int end;

float dxk;
float dyk;
float dzk;
float kde_dst_pow2;



#pragma scop

for (int l = 0; l < lig_natom; ++l) {	// lig loop, ~30
    ekde1 = 0.0f;
    lig__tk = lig_t[l];

    begin = kde_begin_idx[l];
    end = kde_end_idx[l];

  for (int k = begin; k < end; ++k) {	// kde loop, ~400
      dxk = lig_x1[l] - kde->x[k];
      dyk = lig_y1[l] - kde->y[k];
      dzk = lig_z1[l] - kde->z[k];
      kde_dst_pow2 = dxk * dxk + dyk * dyk + dzk * dzk;
      ekde1 += expf (enepara_kde2 * kde_dst_pow2);
  }

  ekde += ekde1 / (float) (end - begin);
}

e_s[6] = ekde * enepara_kde3_inv;


#pragma endscop


