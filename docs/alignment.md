The tool implements affine gap dynamic programming.

### Gap Model

- GAP_OPEN = -5
- GAP_EXTEND = -1

### States

- S: match/mismatch
- E: gap in X
- F: gap in Y

### Recurrence

Each cell considers:

- Match/mismatch
- Gap open
- Gap extension

### Modes

- Global (Needleman–Wunsch)
- Local (Smith–Waterman)

Local alignment resets negative scores to zero.

### Traceback

- Encoded via pointer matrix
- Supports reconstruction of alignment path