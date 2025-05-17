## Output File Inspection 
 
### 1. Check the output file for the expected number of rows and columns 

```bash
printf "file\trows\tcols\n"                                                                  
for f in dna/*_dp_matrix.txt dna/lcs_traceback.txt; do
  r=$(grep -v '^[[:space:]]*$' "$f" | wc -l)
  c=$(grep -v '^[[:space:]]*$' "$f" | head -n 1 | awk '{print NF}')
  printf "%s\t%d\t%d\n" "$f" "$r" "$c"
done | column -t
```

### 2. **Check for non-numeric or invalid values in DP matrices**

Ensure all values are integers (no NaNs, garbage, or floating-point values):

```bash
for f in dna/*_dp_matrix.txt; do
  if grep -Pvq '^(\s*-?\d+\s*)+$' "$f"; then
    echo "Invalid values found in $f"
  else
    echo "Valid integers in $f"
  fi
done
```

---

### 3. **Check for uniform column count in each row**

Ensure each row in the matrix has the same number of columns:

```bash
for f in dna/*_dp_matrix.txt dna/lcs_traceback.txt; do
  awk '{print NF}' "$f" | uniq -c | wc -l | grep -q '^1$' \
    && echo "Uniform columns in $f" \
    || echo "Column inconsistency in $f"
done
```

---

### 4. **Check for empty or very small files**

Useful to detect broken outputs:

```bash
for f in dna/*_dp_matrix.txt dna/lcs_traceback.txt; do
  [[ ! -s "$f" ]] && echo "Empty file: $f"
done
```

---

## Alignment Logic Checks

### 5. **Check for alignment score consistency**

Ensure the score in `dp[m][n]` is equal to the final alignment score printed/stored elsewhere. You could compare it with your JSON stats if available:

```bash
tail -n 1 dna/global_dp_matrix.txt | awk '{print $NF}'  # score at bottom-right corner
```

Compare with `global_stats.json`. 

```bash 
jq -r '.global_score' dna/global_stats.json
``` 
If they match, you can be more confident in the correctness of your DP matrix.

```bash 
if [[ $(tail -n 1 dna/global_dp_matrix.txt | awk '{print $NF}') == $(jq -r '.global_score' dna/global_stats.json) ]]; then
  echo "Alignment score matches"
else
  echo "Alignment score mismatch"
fi
```

---

## Traceback-Specific Checks

### 6. **Check only valid traceback characters (U, L, D, space)**

```bash
for f in dna/lcs_traceback.txt; do
  if grep -Pvq '^[ULD ]*$' "$f"; then
    echo "Invalid characters in traceback: $f"
  else
    echo "Valid characters in traceback: $f"
  fi
done
```

---

### 7. **Count direction types in traceback matrix**

```bash
awk '{for(i=1;i<=NF;++i) a[$i]++} END{for(k in a) print k, a[k]}' dna/lcs_traceback.txt
```

Helpful to understand distribution of U/L/D (Up/Left/Diagonal) directions.

---

### 8. **Check matrix is square if expected**

Only applies if your DP matrices must be square:

```bash
for f in dna/*_dp_matrix.txt; do
  rows=$(wc -l < "$f")
  cols=$(awk 'NR==1{print NF;exit}' "$f")
  if [ "$rows" -eq "$cols" ]; then
    echo "✅ Square matrix in $f"
  else
    echo "ℹ️  Non-square matrix in $f: ${rows}x${cols}"
  fi
done
```

