#!/bin/bash
# Builds swmm5-infosewer into bin/swmm5_is, and an unmodified EPA SWMM 5.2.4 into
# bin/swmm5_epa for the control runs. The build id is the first 12 hex digits of the
# SHA-256 of every engine source file in name order; it is printed in the report.
# -ffp-contract=off keeps the stock routing paths byte-identical to EPA 5.2.4.
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/src"
B=$(find solver run -type f \( -name '*.c' -o -name '*.h' -o -name '*.dat' \) | LC_ALL=C sort | xargs cat | sha256sum | cut -c1-12)
mkdir -p "$ROOT/bin"
FLAGS="-O2 -g -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0 -ffp-contract=off -w"
gcc $FLAGS -DISEWER_BUILD="\"$B\"" -Isolver -Isolver/include -o "$ROOT/bin/swmm5_is" run/main.c solver/*.c -lm
echo "build $B: $ROOT/bin/swmm5_is"
if [ -d "$ROOT/epa524/solver" ] && [ ! -x "$ROOT/bin/swmm5_epa" ]; then
  gcc $FLAGS -I"$ROOT/epa524/solver" -I"$ROOT/epa524/solver/include" -o "$ROOT/bin/swmm5_epa" run/main.c "$ROOT"/epa524/solver/*.c -lm
  echo "control: $ROOT/bin/swmm5_epa (EPA SWMM 5.2.4, unmodified)"
fi
