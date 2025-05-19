#!/usr/bin/env bash
set -euo pipefail

# ─── CONFIG ──────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MATRIX="$SCRIPT_DIR/../results/local_dp_matrix.txt"
TILE_DIR="tiles"                 # where to write tile .dat files
PNG_DIR="tiles/pngs"             # where to write tile .png files
FINAL="$SCRIPT_DIR/../results/local_dp_heatmap.png"

TR=500         # tile rows
TC=500         # tile cols
GP_CORES=$(lscpu -p=CORE \
  | grep -v '^#' \
  | sort -u \
  | wc -l)

# ─── PREREQS CHECK ───────────────────────────────────────────────
command -v parallel >/dev/null || { echo "Install GNU parallel"; exit 1; }
command -v gnuplot  >/dev/null || { echo "Install gnuplot";  exit 1; }
command -v montage  >/dev/null || { echo "Install ImageMagick"; exit 1; }

# ─── PREPARE DIRS ────────────────────────────────────────────────
rm -rf "$TILE_DIR" "$PNG_DIR"    # clean old intermediates
mkdir -p "$TILE_DIR" "$PNG_DIR"

# ─── READ DIMENSIONS ─────────────────────────────────────────────
# N = number of columns, M = number of rows
N=$(head -n1 "$MATRIX" | awk '{print NF}')
M=$(wc -l < "$MATRIX")

RTILES=$(( (M + TR - 1) / TR ))
CTILES=$(( (N + TC - 1) / TC ))

echo "Matrix: ${M}×${N}, splitting into ${RTILES}×${CTILES} tiles..."

# ─── BUILD TILE SPECS ────────────────────────────────────────────
tilespecs=()
for i in $(seq 0 $((RTILES-1))); do
  r1=$(( i*TR + 1 ))
  r2=$(( (i+1)*TR < M ? (i+1)*TR : M ))
  for j in $(seq 0 $((CTILES-1))); do
    c1=$(( j*TC + 1 ))
    c2=$(( (j+1)*TC < N ? (j+1)*TC : N ))
    tilespecs+=("${r1}:${r2}:${c1}:${c2}:${i}:${j}")
  done
done

export MATRIX TILE_DIR

# ─── SPLIT INTO TILES IN PARALLEL ────────────────────────────────
printf "%s\n" "${tilespecs[@]}" | parallel --bar -j "$GP_CORES" --colsep ':' '
  r1={1}; r2={2}; c1={3}; c2={4}; i={5}; j={6}
  out="$TILE_DIR/tile_${i}_${j}.dat"
  # extract rows r1–r2, then columns c1–c2
  sed -n "${r1},${r2}p" "$MATRIX" | cut -d" " -f"${c1}-${c2}" > "$out"
'

# ─── GNUPLOT TEMPLATE ────────────────────────────────────────────
GP_FILE=$(mktemp)
cat > "$GP_FILE" <<'EOF'
set terminal pngcairo size 600,600
set palette viridis
unset key; unset xtics; unset ytics; unset border
set view map
plot ARG1 matrix with image
EOF

# ─── PLOT TILES IN PARALLEL ──────────────────────────────────────
export GP_FILE PNG_DIR
find "$TILE_DIR" -name 'tile_*.dat' | \
  parallel --bar -j "$GP_CORES" '
    tile={}
    base=$(basename "$tile" .dat)
    gnuplot -e "ARG1=\"$tile\"; set output=\"$PNG_DIR/${base}.png\"" "$GP_FILE"
  '

rm "$GP_FILE"

# ─── STITCH PNGs ─────────────────────────────────────────────────
montage "$PNG_DIR"/tile_*.png \
  -tile ${CTILES}x${RTILES} -geometry +0+0 "$FINAL"

# ─── CLEANUP INTERMEDIATES ───────────────────────────────────────
rm -rf "$TILE_DIR" "$PNG_DIR"

echo "Done: $FINAL"