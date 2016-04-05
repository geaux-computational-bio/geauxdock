
// sparse matrix in Coordinate list (COO) format


for (int i = threadIdx.x; i < mcs_nrow; i += blockDim.x)
elhm1s[i] = 0.0f;

for (int i = threadIdx.x; i < mcs_npoint; i += blockDim.x) {
    const int idx_col = mcs_csr->idx_col[i];
    const int idx_row = mcs_csr->idx_row[i];
    const float dx = lig_x2[idx_col] - mcs_csr->x[i];
    const float dy = lig_y2[idx_col] - mcs_csr->y[i];
    const float dz = lig_z2[idx_col] - mcs_csr->z[i];
    const float tmp = dx * dx + dy * dy + dz * dz;
    atomicAdd (&elhm1s[idx_row], tmp);
}

__syncthreads ();
float elhm;
float elhm_sum = 0.0f;



for (int i = 0; i < mcs_nrow; i += blockDim.x) {
    const int m = i + threadIdx.x;
    elhm = 0.0f;
    if (m < mcs_nrow)
        elhm = mcs_tcc[m] * sqrtf (elhm1s[m] / (float) mcs_ncol[m]);

    BlockReduceSum_1_d_2 (elhm);

    if (threadIdx.x == 0)
        elhm_sum += elhm;
}

if (threadIdx.x == 0)
    e_s[7] = logf (elhm_sum / mcs_nrow);


