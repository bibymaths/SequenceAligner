#!/bin/bash

# Usage check
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <lcs_traceback_file> <global_dp_matrix.txt> <local_dp_matrix.txt>"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LCS_TRACEBACK_FILE="$1"
GLOBAL_DP_MATRIX="$2"
LOCAL_DP_MATRIX="$3"
OUTDIR="$SCRIPT_DIR/../results"
OUTPREFIX="${OUTDIR}/plot"
STATS_DIR=$(dirname "$GLOBAL_DP_MATRIX")

mkdir -p "$OUTDIR"

awk 'NF {
  for (i=1; i<=NF; i++) {
    if ($i=="U") printf "1 ";
    else if ($i=="L") printf "2 ";
    else if ($i=="D") printf "3 ";
    else printf "0 ";
  }
  printf "\n"
}' "$LCS_TRACEBACK_FILE" > "${OUTPREFIX}_lcs_numeric.txt"

awk '{ lines[NR] = $0 } END { for (i = NR; i > 0; i--) print lines[i] }' "${OUTPREFIX}_lcs_numeric.txt" > "${OUTPREFIX}_lcs_numeric_flip.txt"
awk '{ lines[NR] = $0 } END { for (i = NR; i > 0; i--) print lines[i] }' "$GLOBAL_DP_MATRIX" > "${OUTPREFIX}_global_flip.txt"
awk '{ lines[NR] = $0 } END { for (i = NR; i > 0; i--) print lines[i] }' "$LOCAL_DP_MATRIX" > "${OUTPREFIX}_local_flip.txt"

plot_matrix() {
  local infile="$1"
  local outfile="$2"
  local palette="$3"
  local cbrange="$4"
  local colorlabel="$5"

  gnuplot -persist <<EOF

# set terminal pngcairo size 1000,1000 enhanced font 'Arial,10'
set terminal pngcairo size 1000,1000 enhanced font 'DejaVu-Sans,10'
set output '${outfile}.png'
set datafile separator whitespace
unset key; unset border; unset xtics; unset ytics; unset title
set colorbox
set cblabel "${colorlabel}" font ",12" offset 2,0
set size ratio -1
set margins 0,0,0,0
set view map
$palette
$cbrange
plot '${infile}' matrix with image
EOF
}

plot_matrix "${OUTPREFIX}_lcs_numeric_flip.txt" "${OUTPREFIX}_lcs" \
  "set palette defined (0 'white', 1 'blue', 2 'green', 3 'red')" \
  "set cbrange [0:3]; set cbtics (' ' 0, 'U' 1, 'L' 2, 'D' 3)" \
  "Direction (U=1, L=2, D=3)"

plot_matrix "${OUTPREFIX}_global_flip.txt" "${OUTPREFIX}_global" \
  "set palette rgb 33,13,10" "" "Score"

plot_matrix "${OUTPREFIX}_local_flip.txt" "${OUTPREFIX}_local" \
  "set palette rgb 33,13,10" "" "Score"

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

#montage \
#  "${OUTPREFIX}_global_labeled.png" "${OUTPREFIX}_local_labeled.png" \
#  "${OUTPREFIX}_lcs_labeled.png"    "${OUTPREFIX}_stats.png" \
#  -tile 2x2 -geometry +10+10 -background white \
#  "${OUTDIR}/summary.png"
#
## --- Cleanup ---
#rm -f "${OUTPREFIX}"_*.png "${OUTPREFIX}"_*.txt