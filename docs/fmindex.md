The FM-index enables fast substring search in the target sequence.

### Components

- Suffix Array
- Burrows-Wheeler Transform (BWT)
- C table (cumulative counts)
- Occ table (rank queries)

### Backward Search

Given a pattern:

- Iterate from right to left
- Narrow interval [l, r) in suffix array
- Final interval gives match positions

### Complexity

- Query time: O(m)
- Preprocessing: O(n log n)

### Implementation Notes

- Stored in memory-efficient structures
- Supports serialization
- Used for seed generation