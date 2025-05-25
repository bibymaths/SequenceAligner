import numpy as np

def downsample(mat, target_rows, target_cols):
    rstep = max(1, mat.shape[0] // target_rows)
    cstep = max(1, mat.shape[1] // target_cols)
    trimmed = mat[:(mat.shape[0]//rstep)*rstep, :(mat.shape[1]//cstep)*cstep]
    ds = trimmed.reshape(-1, rstep, trimmed.shape[1]//cstep, cstep).mean(axis=(1,3))
    return ds.astype(np.intc)