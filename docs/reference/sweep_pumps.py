#!/usr/bin/env python3
"""sweep_pumps.py - random pump-station cases for checking swmm5-infosewer against
the InfoSewer engine: random pump speeds (mixed within stations), parallel units
and, in EPS, wet well initial levels, applied to one InfoSewer input file.

  python3 -I sweep_pumps.py <engine> <base.txt> <outdir> <steady|eps> <seed> <count>

<engine> is the InfoSewer engine built from its source (ise_run <in> <rpt> <out> mode).
Each case is written as <outdir>/cNNN.txt with the engine's .out/.rpt beside it.
The script prints whether the engine stopped; run tests/check_steady.py or
tests/check_eps.py on the cases that ran, and compare the stops with swmm5-infosewer's.

Written by Claude (AI) for SWMM Enablement, 2026-11-10.
"""
import random, subprocess, os, sys

engine, base, out, mode, seed, count = sys.argv[1:7]
random.seed(int(seed))
src = open(base, encoding="latin-1").read().split("\n")
os.makedirs(out, exist_ok=True)
pumps = []
sec = None
for l in src:
    if l.startswith("["):
        sec = l.strip()
    elif sec == "[Pumps]" and l.strip():
        pumps.append(l.split("\t")[0])
for n in range(int(count)):
    sp = {p: (round(random.uniform(0.55, 1.25), 3) if random.random() < 0.5 else 1.0) for p in pumps}
    ex = {p: random.choice([0, 0, 0, 1, 2]) for p in pumps}
    lines, sec = [], None
    for l in src:
        if l.startswith("["):
            sec = l.strip()
        t = l.split("\t")
        if sec == "[Pumps]" and t[0] in sp:
            t[11] = "%.6f" % sp[t[0]]
            t[4] = str(ex[t[0]])
        if mode == "eps" and sec == "[Wet Wells]" and len(t) > 6 and random.random() < 0.6:
            t[5] = "%.6f" % random.uniform(float(t[3]), float(t[4]))
        lines.append("\t".join(t))
    b = os.path.join(out, "c%03d" % n)
    open(b + ".txt", "w", encoding="latin-1").write("\n".join(lines))
    r = subprocess.run([engine, b + ".txt", b + ".rpt", b + ".out", mode], capture_output=True, text=True)
    err = [l for l in (r.stdout + r.stderr).splitlines() if "Error" in l][:1]
    print("%s engine %s %s" % (os.path.basename(b), "stopped" if r.returncode else "ran", err[0] if err else ""))
