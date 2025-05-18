
## Theoretical Background

### 1. Longest Common Subsequence (LCS)

The LCS problem is solved using a dynamic programming table to compute the length and path of the longest subsequence common to both strings. The solution provides insights into evolutionary or structural similarity at a character level.

* Time Complexity: $O(m \cdot n)$
* Space Complexity: $O(m \cdot n)$
* Alignment Type: Character-based (not contiguous)

### 2. Global Alignment (Needleman-Wunsch)

This method aligns the entire sequences using a scoring matrix initialized with gap penalties. It is ideal for full-length, homologous sequences.

* Time Complexity: $O(m \cdot n)$
* Uses a traceback matrix to reconstruct the alignment path.
* Supports scoring for matches, mismatches, and gaps.

### 3. Local Alignment (Smith-Waterman)

Focuses on identifying the most similar subsequences. Suitable for aligning protein motifs or conserved regions.

* Time Complexity: $O(m \cdot n)$
* Incorporates zero in scoring to allow for local reset.
* Traceback stops when score becomes zero.

--- 

## Alignment Matrix and Mathematical Explanation

### 1. Longest Common Subsequence (LCS)

**Goal:** Find the longest subsequence common to two sequences $X = x_1x_2\ldots x_m$ and $Y = y_1y_2\ldots y_n$.

---

### Matrix Definition

Let $C[i][j]$ denote the length of the LCS of the prefixes $x_1\ldots x_i$ and $y_1\ldots y_j$.

**Recurrence Relation:**

$$
C[i][j] =
\begin{cases}
0 & \text{if } i = 0 \text{ or } j = 0 \\\\
C[i-1][j-1] + 1 & \text{if } x_i = y_j \\\\
\max(C[i-1][j], C[i][j-1]) & \text{otherwise}
\end{cases}
$$

**Backtracking Table:**

* `D`: Diagonal (match)
* `U`: Up (skip from X)
* `L`: Left (skip from Y)

**In Code:**

```cpp
if (x[i - 1] == y[j - 1]) {
    c[i][j] = c[i - 1][j - 1] + 1;
    b[i][j] = 'D';
} else if (c[i - 1][j] >= c[i][j - 1]) {
    c[i][j] = c[i - 1][j];
    b[i][j] = 'U';
} else {
    c[i][j] = c[i][j - 1];
    b[i][j] = 'L';
}
```

---

### 2. Global Alignment (Needleman–Wunsch)

**Goal:** Align full sequences $X$ and $Y$ from start to end, minimizing penalties or maximizing score.

---

### Matrix Definition

Let $S[i][j]$ denote the optimal alignment score between $x_1\ldots x_i$ and $y_1\ldots y_j$.

**Initialization:**

$$
S[i][0] = i \cdot \text{GAP},\quad S[0][j] = j \cdot \text{GAP}
$$

**Recurrence Relation:**

$$
S[i][j] = \max \begin{cases}
S[i-1][j-1] + \text{score}(x_i, y_j) \\
S[i-1][j] + \text{GAP} \\
S[i][j-1] + \text{GAP}
\end{cases}
$$

where:

$$
\text{score}(x_i, y_j) = 
\begin{cases}
\text{MATCH} & \text{if } x_i = y_j \\\\
\text{MISMATCH} & \text{otherwise}
\end{cases}
$$

**Traceback:** Built from a `traceback[][]` matrix storing `'D'`, `'U'`, `'L'` for path recovery.

---

### 3. Local Alignment (Smith–Waterman)

**Goal:** Find the highest scoring local region between two sequences, ideal for detecting motifs.

---

### Matrix Definition

Let $H[i][j]$ be the highest scoring alignment ending at $x_i, y_j$.

**Initialization:**

$$
H[i][0] = 0,\quad H[0][j] = 0
$$

**Recurrence Relation:**

$$
H[i][j] = \max \begin{cases}
0 \\
H[i-1][j-1] + \text{score}(x_i, y_j) \\
H[i-1][j] + \text{GAP} \\
H[i][j-1] + \text{GAP}
\end{cases}
$$

**Reset to 0** allows alignment to start anywhere and stop at the best score.

**Traceback:** Starts from the max score cell and moves until score is zero.

---

## Time and Space Complexities

| Algorithm   | Time Complexity | Space Complexity | Alignment Type |
| ----------- | --------------- | ---------------- | -------------- |
| LCS         | O(m·n)          | O(m·n)           | Character-only |
| Global (NW) | O(m·n)          | O(m·n)           | Full sequence  |
| Local (SW)  | O(m·n)          | O(m·n)           | Subsequence    |

Where $m = |X|$, $n = |Y|$

---


## Visualizing with `gnuplot`

### Why use it?

* DP matrix and traceback files are large, often matching the product of the sequence lengths, making manual inspection impractical.
* Visual representation helps interpret alignment decisions and traceback paths more easily.
* `gnuplot` is a robust and flexible tool for generating high-quality heatmaps and matrix plots from plain text data.
 
- Install `gnuplot` if not already installed. For example, on Ubuntu: 

```bash
sudo apt-get install gnuplot
```  
 
or on Fedora:  

```bash 
sudo dnf install gnuplot
``` 

or on macOS: 

```bash
brew install gnuplot
``` 
 
- Run the script to generate plot from traceback/DP files  
 
```bash 
chmod +x plotDP.sh 
./plotDP.sh <lcs_traceback_file> <global_dp_file> <local_dp_file> 
```


