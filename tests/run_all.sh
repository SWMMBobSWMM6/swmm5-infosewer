#!/bin/bash
# run_all.sh - build, then run every verification: steady state and EPS against the
# InfoSewer engine's own .out files, and the null test against EPA SWMM 5.2.4.
# Written by Claude (AI) for SWMM Enablement, 2026-11-09.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
bash "$ROOT/build.sh"
echo; echo "== steady state vs InfoSewer engine"; python3 -I "$ROOT/tests/check_steady.py" "$ROOT/tests/steady" "$ROOT/tests/pressure" || st=1
echo; echo "== EPS vs InfoSewer engine";          python3 -I "$ROOT/tests/check_eps.py" "$ROOT/tests/eps" "$ROOT/tests/eps_pressure" || ep=1
echo; echo "== null test vs EPA SWMM 5.2.4";      bash "$ROOT/tests/check_null.sh" || nu=1
echo
if [ -z "$st$ep$nu" ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED (steady=${st:-0} eps=${ep:-0} null=${nu:-0})"; exit 1; fi
