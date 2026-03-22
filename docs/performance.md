### Parallelization

- MPI for distributed computation
- OpenMP for shared-memory parallelism

### SIMD

- Uses AVX instructions via `immintrin.h`

### Optimization Strategies

- Seed filtering reduces DP workload
- Binary matrix output avoids I/O overhead
- Chunked processing of sequences

### Progress Monitoring

- Rank 0 prints progress bar
- ETA estimation included