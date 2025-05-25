#!/bin/bash
set -euo pipefail

# Argument Check
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <lcs_traceback_file> <global_dp_matrix.txt> <local_dp_matrix.txt>"
  exit 1
fi

# Paths and Setup
LCS_TRACEBACK_FILE="$1"
GLOBAL_DP_MATRIX="$2"
LOCAL_DP_MATRIX="$3"
OUTDIR="$4"
OUTPREFIX="${OUTDIR}/plot"
STATS_DIR=$(dirname "$GLOBAL_DP_MATRIX")

mkdir -p "$OUTDIR"

# Flip + Optional Downsampling
flip_and_downsample() {
  local infile="$1"
  local outfile="$2"
  local label="$3"
  local max_dim=50000

  local flip_tmp="${outfile}_flip.tmp"
  awk '{ buf[NR] = $0 } END { for (i = NR; i > 0; i--) print buf[i] }' "$infile" > "$flip_tmp"

  local rows cols
  rows=$(wc -l < "$flip_tmp")
  cols=$(awk 'NR==1 { print NF }' "$flip_tmp")

  local rskip=$(( (rows + max_dim - 1) / max_dim ))
  local cskip=$(( (cols + max_dim - 1) / max_dim ))

  if [ "$rskip" -gt 1 ] || [ "$cskip" -gt 1 ]; then
    echo "Downsampling $label matrix (${rows}x${cols}) → (~$((rows/rskip))x$((cols/cskip)))"
    awk -v rskip="$rskip" -v cskip="$cskip" '
      NR % rskip == 0 {
        for (i = 1; i <= NF; i++) {
          if (i % cskip == 0) printf "%s ", $i;
        }
        print "";
      }
    ' "$flip_tmp" > "$outfile"
    rm "$flip_tmp"
  else
    mv "$flip_tmp" "$outfile"
  fi
}

# Prepare Matrices
awk 'NF {
  for (i=1; i<=NF; i++) {
    if ($i=="U") printf "1 ";
    else if ($i=="L") printf "2 ";
    else if ($i=="D") printf "3 ";
    else printf "0 ";
  }
  printf "\n"
}' "$LCS_TRACEBACK_FILE" > "${OUTPREFIX}_lcs_numeric.txt"

flip_and_downsample "${OUTPREFIX}_lcs_numeric.txt" "${OUTPREFIX}_lcs_final.txt" "LCS"
flip_and_downsample "$GLOBAL_DP_MATRIX" "${OUTPREFIX}_global_final.txt" "Global DP"
flip_and_downsample "$LOCAL_DP_MATRIX" "${OUTPREFIX}_local_final.txt" "Local DP"

flip_path_y() {
  local pathfile="$1"
  local flipped_path="$2"
  local max_y

  max_y=$(awk '{if ($2 > max) max=$2} END{print max}' "$pathfile")
  awk -v max="$max_y" '{ print $1, max - $2 + 1 }' "$pathfile" > "$flipped_path"
}

plot_matrix() {
  local infile="$1"
  local outfile="$2"
  local palette="$3"
  local cbrange="$4"
  local colorlabel="$5"
  local pathfile="${6:-}"

  local overlay=""
  if [ -n "$pathfile" ] && [ -f "$pathfile" ]; then
    overlay=", '${pathfile}' using 1:2 with linespoints lw 2 pt 7 ps 1.0 lc rgb \"black\""
  fi

  gnuplot <<EOF
set terminal pngcairo size 1000,1000 enhanced font 'Helvetica,10'
set output '${outfile}.png'
unset key; unset border; unset xtics; unset ytics; unset title
set colorbox
set cblabel "${colorlabel}" font ",12" offset 2,0
set size ratio -1
set margins 0,0,0,0
set view map
$palette
$cbrange
plot '${infile}' matrix with image${overlay}
EOF

  if [ ! -s "${outfile}.png" ]; then
    echo "Failed to create ${outfile}.png"
    exit 1
  fi
}

# Plot DP Matrices
flip_path_y "${STATS_DIR}/lcs_path.txt" "${OUTPREFIX}_lcs_path_flipped.txt"
plot_matrix "${OUTPREFIX}_lcs_final.txt" "${OUTPREFIX}_lcs" \
  "set palette defined (0 'white', 1 '#ADD8E6', 2 '#90EE90', 3 '#FFB6C1')" \
  "set cbrange [0:3]; set cbtics (' ' 0, 'U' 1, 'L' 2, 'D' 3)" \
  "Direction (U=1, L=2, D=3)" \
  "${OUTPREFIX}_lcs_path_flipped.txt"

flip_path_y "${STATS_DIR}/global_path.txt" "${OUTPREFIX}_global_path_flipped.txt"
plot_matrix "${OUTPREFIX}_global_final.txt" "${OUTPREFIX}_global" \
  "set palette rgb 33,13,10" "" "Score" \
  "${OUTPREFIX}_global_path_flipped.txt"

flip_path_y "${STATS_DIR}/local_path.txt" "${OUTPREFIX}_local_path_flipped.txt"
plot_matrix "${OUTPREFIX}_local_final.txt" "${OUTPREFIX}_local" \
  "set palette rgb 33,13,10" "" "Score" \
  "${OUTPREFIX}_local_path_flipped.txt"

# Format Alignment Stats
format_stats() {
  local file="$1"
  local title="$2"
  jq -r --arg title "$title" '
    "\($title)\n\n" +
    "Query:      \(.queryid)\n" +
    "Target:     \(.targetid)\n" +
    "Score:      \(.score)\n" +
    "Matches:    \(.matches)\n" +
    "Gaps:       \(.gaps)\n" +
    "Total:      \(.total)\n" +
    "Identity:   \((.identity * 100) | round)%\n" +
    "Coverage:   \((.coverage * 100) | round)%\n" +
    "Time (ms):  \(.time_ms)"
  ' "$file"
}

format_stats "${STATS_DIR}/global_stats.json" "Global Alignment Stats" > "${OUTPREFIX}_stats.txt"
echo -e "\n" >> "${OUTPREFIX}_stats.txt"
format_stats "${STATS_DIR}/local_stats.json" "Local Alignment Stats" >> "${OUTPREFIX}_stats.txt"

# Annotate Images
for type in global local; do
  label=$(echo "$type" | tr '[:lower:]' '[:upper:]')
  magick "${OUTPREFIX}_${type}.png" \
    -gravity north -pointsize 28 -annotate +0+40 "${label} DP Matrix" \
    "${OUTPREFIX}_${type}_labeled.png"
done

magick "${OUTPREFIX}_lcs.png" \
  -gravity north -background white -splice 0x80 \
  -gravity north -pointsize 28 -annotate +0+20 "LCS Traceback" \
  "${OUTPREFIX}_lcs_labeled.png"

magick -size 1000x1000 xc:white \
  -gravity northwest -pointsize 22 \
  -annotate +50+100 "@${OUTPREFIX}_stats.txt" \
  "${OUTPREFIX}_stats.png"

# Montage and Final Check
for img in "${OUTPREFIX}"_{global_labeled,local_labeled,lcs_labeled,stats}.png; do
  if [ ! -s "$img" ]; then
    echo "Missing image: $img"
    exit 1
  fi
done

montage \
  "${OUTPREFIX}_global_labeled.png" "${OUTPREFIX}_local_labeled.png" \
  "${OUTPREFIX}_lcs_labeled.png" "${OUTPREFIX}_stats.png" \
  -tile 2x2 -geometry +10+10 -background white \
  "${OUTDIR}/summary.png"

# Cleanup
rm -f "${OUTPREFIX}"_lcs_numeric.txt
rm -f "${OUTPREFIX}"_*_final.txt
rm -f "${OUTPREFIX}"_global.png "${OUTPREFIX}"_local.png "${OUTPREFIX}"_lcs.png
rm -f "${OUTPREFIX}"_global_labeled.png "${OUTPREFIX}"_local_labeled.png "${OUTPREFIX}"_lcs_labeled.png "${OUTPREFIX}_stats.png"
rm -f "${OUTPREFIX}"_*_path_flipped.txt