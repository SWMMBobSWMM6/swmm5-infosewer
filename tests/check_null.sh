#!/bin/bash
# check_null.sh - the null test: with stock routing options, swmm5-infosewer must give
# the same .out bytes and the same .rpt (apart from clock lines) as unmodified EPA SWMM
# 5.2.4 (bin/swmm5_epa, built by build.sh from epa524/).
# Usage: tests/check_null.sh [model.inp ...]   (default: tests/null/*.inp)
# Written by Claude (AI) for SWMM Enablement, 2026-11-09.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
W="$ROOT/tests/work_null"; mkdir -p "$W"
files=("$@"); [ ${#files[@]} -eq 0 ] && files=("$ROOT"/tests/null/*.inp)
n=0; bad=0
for f in "${files[@]}"; do
  b=$(basename "$f" .inp)
  "$ROOT/bin/swmm5_is"  "$f" "$W/$b.is.rpt"  "$W/$b.is.out"  >/dev/null 2>&1
  "$ROOT/bin/swmm5_epa" "$f" "$W/$b.epa.rpt" "$W/$b.epa.out" >/dev/null 2>&1
  n=$((n+1))
  if cmp -s "$W/$b.is.out" "$W/$b.epa.out" && \
     diff -q <(grep -v "Analysis begun\|Analysis ended\|Total elapsed" "$W/$b.is.rpt") \
             <(grep -v "Analysis begun\|Analysis ended\|Total elapsed" "$W/$b.epa.rpt") >/dev/null; then
    echo "$b: identical to EPA SWMM 5.2.4"
  else
    echo "$b: DIFFERS from EPA SWMM 5.2.4"; bad=$((bad+1))
  fi
done
echo "$((n-bad)) of $n models identical"
[ $bad -eq 0 ]
