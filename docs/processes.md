# Mathematical Foundations and Implementation Model

This page describes the actual mathematics implemented in `SequenceAligner`, based on the codebase rather than a generic
alignment summary. The implementation combines classical dynamic programming with indexed seed discovery and chaining.
The main algorithmic parts are:

1. suffix array construction
2. FM-index construction and backward search
3. k-mer seed generation
4. seed chaining with affine-style gap costs
5. affine-gap global alignment
6. affine-gap local alignment in candidate windows
7. longest common subsequence (LCS)

The code supports both DNA and protein scoring modes through EDNAFULL and BLOSUM62-style substitution tables,
respectively.

---

## 1. Sequences, notation, and scoring

Let

- \(X = x_1x_2\cdots x_m\) be the query sequence
- \(Y = y_1y_2\cdots y_n\) be the target sequence

with lengths \(m = |X|\) and \(n = |Y|\).

The implementation defines a substitution score function

\[
\sigma(a,b)
\]

which depends on the selected mode:

- DNA mode: EDNAFULL-derived scoring
- protein mode: BLOSUM62-derived scoring

In the code, this is dispatched through:

\[
\sigma(a,b)=
\begin{cases}
\sigma_{\mathrm{DNA}}(a,b), & \text{DNA mode}\\[4pt]
\sigma_{\mathrm{protein}}(a,b), & \text{protein mode}
\end{cases}
\]

Gap penalties are affine:

- gap open: \(g_o = \texttt{GAP\_OPEN} = -5\)
- gap extend: \(g_e = \texttt{GAP\_EXTEND} = -1\)

so that a gap of length \(L \ge 1\) has score

\[
g(L) = g_o + (L-1)g_e.
\]

These constants are explicitly defined in the code.

---

## 2. Suffix array construction

For a string \(T = t_0t_1\cdots t_{N-1}\), the suffix array is the permutation

\[
\mathrm{SA}[0],\mathrm{SA}[1],\ldots,\mathrm{SA}[N-1]
\]

such that

\[
T[\mathrm{SA}[0]:] < T[\mathrm{SA}[1]:] < \cdots < T[\mathrm{SA}[N-1]:]
\]

in lexicographic order.

The implementation uses a rank-doubling construction. Initially, each suffix starting at position \(i\) is ranked by the
character code of \(T_i\). At iteration \(k=1,2,4,\dots\), suffixes are sorted by the pair

\[
\bigl(r_k(i),\, r_k(i+k)\bigr),
\]

where \(r_k(i)\) is the current rank of the suffix beginning at \(i\). If \(i+k \ge N\), the second component is treated
as \(-1\).

After sorting, new ranks are assigned:

\[
r_{2k}( \mathrm{SA}[0]) = 0,
\]

and for \(j \ge 1\),

\[
r_{2k}(\mathrm{SA}[j]) =
r_{2k}(\mathrm{SA}[j-1]) +
\mathbf{1}\!\left[
\bigl(r_k(\mathrm{SA}[j]), r_k(\mathrm{SA}[j]+k)\bigr)
\neq
\bigl(r_k(\mathrm{SA}[j-1]), r_k(\mathrm{SA}[j-1]+k)\bigr)
\right].
\]

The process stops once all ranks are unique. This is the standard \(O(N\log N)\) rank-doubling method implemented in the
code.

---

## 3. FM-index construction

To index the target sequence, the implementation builds an FM-index on

\[
T = Y\$,
\]

where \(\$\) is a sentinel smaller than all biological alphabet symbols.

### 3.1 Burrows-Wheeler transform

Given the suffix array \(\mathrm{SA}\), the Burrows-Wheeler transform \(B\) is

\[
B[i] =
\begin{cases}
T[N-1], & \mathrm{SA}[i]=0,\\[4pt]
T[\mathrm{SA}[i]-1], & \mathrm{SA}[i]>0.
\end{cases}
\]

### 3.2 C-table

For each character \(c\), the FM-index stores

\[
C[c] = \#\{\,a \in T : a < c\,\},
\]

the total number of symbols in \(T\) that are lexicographically smaller than \(c\).

### 3.3 Occ table

For each character \(c\), the occurrence table stores

\[
\mathrm{Occ}(c,i)=\#\{\,k < i : B[k]=c\,\},
\]

that is, the number of occurrences of \(c\) in the prefix \(B[0:i-1]\).

### 3.4 Backward search

For a pattern \(P = p_1p_2\cdots p_\ell\), matching intervals are computed right-to-left. Start with

\[
[l,r) = [0,N).
\]

Then for \(i=\ell,\ell-1,\dots,1\),

\[
l \leftarrow C[p_i] + \mathrm{Occ}(p_i,l),
\qquad
r \leftarrow C[p_i] + \mathrm{Occ}(p_i,r).
\]

If at any step \(l \ge r\), the pattern does not occur. Otherwise, after the last character, the suffix-array
interval \([l,r)\) identifies all matching locations of the pattern in the target sequence. The code then maps that
interval back to text coordinates using `sa[i]`.

---

## 4. Seed generation from query k-mers

The code generates exact-match seeds by sliding a k-mer window across the query. For a chosen seed length \(k\), define
query k-mers

\[
K_i = X[i:i+k-1], \qquad i=0,1,\dots,m-k.
\]

For each \(K_i\), FM-index backward search returns all matching target positions

\[
\mathcal{P}_i = \{t_1,t_2,\dots\}.
\]

Each match yields a seed

\[
s = (q,t,k),
\]

where

- \(q\) = query start position
- \(t\) = target start position
- \(k\) = seed length

So the raw seed set is

\[
\mathcal{S} = \{(q,t,k): X[q:q+k-1] = Y[t:t+k-1]\}.
\]

This is exactly what `generate_raw_seeds` constructs. 

---

## 5. Seed chaining

A seed \(s_i\) has coordinates

\[
s_i = (q_i,t_i,\ell_i),
\]

with query end and target end

\[
q_i^{\mathrm{end}} = q_i + \ell_i - 1,
\qquad
t_i^{\mathrm{end}} = t_i + \ell_i - 1.
\]

The seeds are sorted by query position, then target position. The chaining step computes a best monotone chain under
overlap, gap, and diagonal-consistency constraints.

### 5.1 Feasibility constraints

A seed \(s_j\) may precede \(s_i\) only if:

1. **no overlap in query**
   \[
   q_j^{\mathrm{end}} + \delta_{\min} < q_i
   \]

2. **no overlap in target**
   \[
   t_j^{\mathrm{end}} + \delta_{\min} < t_i
   \]

3. **bounded inter-seed gaps**
   \[
   0 \le d_q \le \delta_{\max},
   \qquad
   0 \le d_t \le \delta_{\max}
   \]
   where
   \[
   d_q = q_i - q_j^{\mathrm{end}} - 1,
   \qquad
   d_t = t_i - t_j^{\mathrm{end}} - 1
   \]

4. **diagonal consistency**
   \[
   \left| (q_i-t_i) - (q_j-t_j) \right| \le \Delta_{\max}.
   \]

These correspond to `min_diag_gap_val`, `max_diag_gap_val`, and `max_offset_dev_val` in the code.

### 5.2 Chain score

Each seed contributes its length:

\[
w_i = \ell_i.
\]

The gap cost between two consecutive seeds is computed independently in query and target using the same affine penalty
parameters:

\[
\gamma_q(d_q) =
\begin{cases}
0, & d_q=0\\
g_o + (d_q-1)g_e, & d_q>0
\end{cases}
\]

\[
\gamma_t(d_t) =
\begin{cases}
0, & d_t=0\\
g_o + (d_t-1)g_e, & d_t>0
\end{cases}
\]

and the total inter-seed cost is

\[
\Gamma(j,i)=\gamma_q(d_q)+\gamma_t(d_t).
\]

The DP recurrence implemented in `find_best_seed_chain` is

\[
D[i] = w_i + \max\left(
0,\;
\max_{j<i,\; s_j \prec s_i}
\bigl(D[j] - \Gamma(j,i)\bigr)
\right).
\]

Equivalently, in the code’s update form,

\[
D[i] = \max\left(D[i],\; D[j] + w_i - \Gamma(j,i)\right).
\]

The best chain score is

\[
\max_i D[i].
\]

Backtracking through the predecessor array reconstructs the ordered chain of anchors.
{index=9}

---

## 6. Affine-gap global alignment

Global alignment is implemented in `align_segment_globally` using three DP states.

For prefixes \(X_{1:i}\) and \(Y_{1:j}\), define:

- \(S_{i,j}\): best score ending in a match/mismatch state
- \(E_{i,j}\): best score ending with a gap in sequence 1 relative to sequence 2 in the code’s traceback convention
- \(F_{i,j}\): best score ending with a gap in sequence 2 relative to sequence 1

The code uses row-compressed score arrays and a full traceback matrix.

### 6.1 Initialization

The affine boundary conditions are:

\[
S_{0,0}=0,\qquad E_{0,0}=-\infty,\qquad F_{0,0}=-\infty.
\]

Along the first row:

\[
E_{0,j}=
\begin{cases}
g_o, & j=1,\\
E_{0,j-1}+g_e, & j>1,
\end{cases}
\qquad
S_{0,j}=E_{0,j},
\qquad
F_{0,j}=-\infty.
\]

Along the first column:

\[
F_{i,0}=
\begin{cases}
g_o, & i=1,\\
F_{i-1,0}+g_e, & i>1,
\end{cases}
\qquad
S_{i,0}=F_{i,0},
\qquad
E_{i,0}=-\infty.
\]

### 6.2 Recurrence

For \(i \ge 1, j \ge 1\),

\[
M_{i,j} =
\max\{S_{i-1,j-1}, E_{i-1,j-1}, F_{i-1,j-1}\}

+ \sigma(x_i,y_j).
  \]

Gap state updates are

\[
E_{i,j} = \max\{S_{i,j-1}+g_o,\; E_{i,j-1}+g_e\},
\]

\[
F_{i,j} = \max\{S_{i-1,j}+g_o,\; F_{i-1,j}+g_e\}.
\]

Then

\[
S_{i,j} = \max\{M_{i,j}, E_{i,j}, F_{i,j}\}.
\]

The final global score is

\[
\mathrm{Score}_{\mathrm{global}} = S_{m,n}.
\]

### 6.3 Traceback

The traceback matrix stores one of:

- `M` for match/mismatch transition
- `E` or `e` for gap-open or gap-extend in one direction
- `F` or `f` for gap-open or gap-extend in the other direction

Backtracking from \((m,n)\) reconstructs aligned strings

\[
\hat X,\hat Y
\]

with insertions represented as `-`. This is exactly the mechanism implemented by `align_segment_globally`.

---

## 7. Affine-gap local alignment in windows

Local alignment is implemented in `perform_sw_in_window`. It uses a Smith–Waterman-style local objective with affine
gaps and three score matrices \(S,E,F\). The local alignment is not necessarily run on the full sequences; instead, the
code uses candidate windows derived from seeds/anchors. 

### 7.1 State definition

For a query substring \(X'\) and target substring \(Y'\),

- \(S_{i,j}\): best local score ending at \((i,j)\)
- \(E_{i,j}\): best local score ending with a horizontal affine gap state
- \(F_{i,j}\): best local score ending with a vertical affine gap state

### 7.2 Initialization

All local matrices start at zero:

\[
S_{0,j}=0,\quad S_{i,0}=0,
\qquad
E_{0,j}=0,\quad E_{i,0}=0,
\qquad
F_{0,j}=0,\quad F_{i,0}=0.
\]

### 7.3 Recurrence

For \(i,j \ge 1\),

\[
M_{i,j} =
\max\{S_{i-1,j-1}, E_{i-1,j-1}, F_{i-1,j-1}\}

+ \sigma(x_i,y_j).
  \]

Gap states are floored at zero:

\[
E_{i,j} = \max\{0,\; S_{i,j-1}+g_o,\; E_{i,j-1}+g_e\},
\]

\[
F_{i,j} = \max\{0,\; S_{i-1,j}+g_o,\; F_{i-1,j}+g_e\}.
\]

The local score state is also floored at zero:

\[
S_{i,j} = \max\{0,\; M_{i,j},\; E_{i,j},\; F_{i,j}\}.
\]

The best local alignment score is

\[
\mathrm{Score}_{\mathrm{local}} = \max_{i,j} S_{i,j}.
\]

### 7.4 Local traceback

Traceback starts from the cell \((i^\*,j^\*)\) where \(S_{i,j}\) is maximal and stops once the running state reaches 0.
This produces the optimal local alignment within that window. The implementation then maps window-relative coordinates
back to original sequence coordinates using stored offsets.

---

## 8. Longest common subsequence (LCS)

The code also implements exact longest common subsequence using standard dynamic programming in
`compute_lcs_for_segment`. This is distinct from score-based alignment: it maximizes subsequence length, not
substitution score. 

### 8.1 DP table

Let \(L_{i,j}\) be the LCS length of prefixes \(X_{1:i}\) and \(Y_{1:j}\). Then

\[
L_{0,j}=0,\qquad L_{i,0}=0.
\]

For \(i,j \ge 1\),

\[
L_{i,j}=
\begin{cases}
L_{i-1,j-1}+1, & x_i=y_j,\\[4pt]
\max(L_{i-1,j},L_{i,j-1}), & x_i\neq y_j.
\end{cases}
\]

### 8.2 Backpointer matrix

The implementation stores directions:

- `D` for diagonal match
- `U` for up
- `L` for left

Backtracking reconstructs:

1. the LCS string itself
2. a gapped representation of both sequences consistent with the chosen traceback

Thus the final outputs are:

- `lcs_string`
- `gapped_seq1`
- `gapped_seq2`

with

\[
| \text{lcs\_string} | = L_{m,n}.
\]

### 8.3 Anchored LCS

For long inputs, the code optionally uses FM-index-derived anchors first. The full problem is decomposed into
inter-anchor segments:

\[
(X_0,Y_0), (X_1,Y_1), \dots, (X_r,Y_r)
\]

LCS is computed on each segment independently, while the anchor substrings themselves are inserted directly into the
final result. If anchor \(a_k\) has length \(\ell_k\), then the final total LCS length is

\[
L_{\mathrm{total}}=
\sum_{s=0}^{r} L^{(s)}

+

\sum_{k=1}^{r} \ell_k,
\]

where \(L^{(s)}\) is the LCS length of segment \(s\). This is exactly how reconstruction is performed in the anchored
LCS routine. 

---

## 9. Anchored decomposition of alignment problems

For both global and LCS workflows, the code may first find a best seed chain and then split the full sequences into
inter-anchor segments.

Suppose the chosen anchor chain is

\[
a_1,a_2,\dots,a_r
\]

with anchor \(a_k=(q_k,t_k,\ell_k)\).

Define segment boundaries recursively:

- segment 0: from sequence start to just before \(a_1\)
- segment \(k\): between \(a_k\) and \(a_{k+1}\)
- final segment: after \(a_r\)

These segment pairs are aligned independently, then concatenated with the exact anchor substrings themselves. If \((\hat
X^{(s)},\hat Y^{(s)})\) is the alignment of segment \(s\), the final reconstructed global-style alignment is

\[
\hat X =
\hat X^{(0)} \,\Vert\, A_1 \,\Vert\, \hat X^{(1)} \,\Vert\, A_2 \,\Vert\, \cdots \,\Vert\, A_r \,\Vert\, \hat X^{(r)},
\]

\[
\hat Y =
\hat Y^{(0)} \,\Vert\, A_1 \,\Vert\, \hat Y^{(1)} \,\Vert\, A_2 \,\Vert\, \cdots \,\Vert\, A_r \,\Vert\, \hat Y^{(r)},
\]

where \(A_k\) is the exact anchor substring, copied identically into both aligned outputs.
{index=16}

---

## 10. Matrix output and serialization

The code supports exporting DP and traceback-style matrices in both text and binary form.

### 10.1 Integer DP matrix

If \(D \in \mathbb{Z}^{r \times c}\) is a DP matrix, the binary writer stores:

1. number of rows
2. number of columns
3. row-major matrix entries

so the serialized payload is conceptually

\[
(r,c,D_{0,0},D_{0,1},\dots,D_{r-1,c-1}).
\]

### 10.2 Character matrix

For traceback matrices \(T \in \Sigma^{r \times c}\), the character writer stores the same leading dimensions and then
the row-major flattened character data, padding rows if needed. This matches the matrix export utilities in the code. 