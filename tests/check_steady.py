#!/usr/bin/env python3
"""check_steady.py - run swmm5-infosewer in steady state on InfoSewer engine inputs
and compare every steady-state output array with the InfoSewer engine's own .out.

  python3 -I check_steady.py <dir-or-txt>... [--bin bin/swmm5_is] [-v]

For each X.txt with an engine result X.out beside it (or ref_X.out), the script
converts X.txt with tools/ise2inp.py, runs bin/swmm5_is with
INFOSEWER_RESULTS_FILE, and compares the 28 arrays (20 pipe, 6 manhole,
2 outlet) value by value. Prints the largest relative difference per model.

Written by Claude (AI) for SWMM Enablement, 2026-11-09.
"""
import sys, os, struct, subprocess, glob, csv

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PIPE_ARRAYS = ["flow", "unpeak", "peak", "cover", "infil", "storm", "flowtype", "full",
               "rescap", "vel", "dD", "qQ", "depth", "crit", "froude", "covcnt",
               "backwater", "adjdepth", "adjvel", "partial"]
NODE_ARRAYS = ["base", "static", "nstorm", "head", "surch", "jump"]


def read_out(path, P, N, O):
    # static .out layout: header, 20 pipe arrays, 6 loading-manhole arrays, 2 outlet arrays
    b = open(path, "rb").read()
    H = len(b) - 8 * (20 * P + 6 * N + 2 * O)
    arr = lambda off, cnt: list(struct.unpack("<%dd" % cnt, b[off:off + 8 * cnt]))
    R = {k: arr(H + i * P * 8, P) for i, k in enumerate(PIPE_ARRAYS)}
    R.update({k: arr(H + 20 * P * 8 + i * N * 8, N) for i, k in enumerate(NODE_ARRAYS)})
    o = arr(H + 20 * P * 8 + 6 * N * 8, 2 * O)
    R["oflow"], R["ohead"] = o[:O], o[O:]
    return R


def read_csv(path):
    pipes, nodes, outs = [], [], []
    hdr = None
    for row in csv.reader(open(path)):
        if row[0] == "kind":
            hdr = row
            continue
        d = dict(zip(hdr[2:], (float(x) for x in row[2:])))
        d["id"] = row[1]
        {"pipe": pipes, "node": nodes, "outlet": outs}[row[0]].append(d)
    R = {k: [p[k] for p in pipes] for k in PIPE_ARRAYS}
    R.update({k: [n[k] for n in nodes] for k in NODE_ARRAYS})
    R["oflow"] = [o["oflow"] for o in outs]
    R["ohead"] = [o["ohead"] for o in outs]
    ids = ([p["id"] for p in pipes], [n["id"] for n in nodes], [o["id"] for o in outs])
    return R, ids


def check(txt, ref, exe, work, verbose):
    base = os.path.splitext(os.path.basename(txt))[0]
    inp = os.path.join(work, base + ".inp")
    rpt = os.path.join(work, base + ".rpt")
    out = os.path.join(work, base + ".out")
    res = os.path.join(work, base + ".csv")
    if os.path.exists(res):
        os.unlink(res)
    subprocess.run([sys.executable, "-I", os.path.join(ROOT, "tools", "ise2inp.py"), txt, inp,
                    "--mode", "steady", "--results", res], check=True, capture_output=True)
    r = subprocess.run([exe, inp, rpt, out], capture_output=True, text=True)
    if not os.path.exists(res):
        print("%-14s RUN FAILED: %s" % (base, (r.stdout + r.stderr).strip()[-300:]))
        os.system("grep -i error '%s' | head -5" % rpt)
        return None
    R, (pid, nid, oid) = read_csv(res)
    E = read_out(ref, len(pid), len(nid), len(oid))
    worst, worst_key = 0.0, ""
    nbad = 0
    for key in PIPE_ARRAYS + NODE_ARRAYS + ["oflow", "ohead"]:
        ids = pid if key in PIPE_ARRAYS else nid if key in NODE_ARRAYS else oid
        for i, (x, y) in enumerate(zip(R[key], E[key])):
            d = abs(x - y)
            rel = d / abs(y) if abs(y) > 1e-12 else d
            if rel > worst:
                worst, worst_key = rel, "%s %s (%.12g vs %.12g)" % (key, ids[i], x, y)
            if rel > 0:
                nbad += 1
                if verbose:
                    print("   %-8s %-24s swmm5_is %.17g engine %.17g" % (key, ids[i], x, y))
    nval = sum(len(v) for v in R.values())
    si = any(l.split(":", 1)[1].strip().upper() in ("LPS", "MLD")
             for l in open(txt, encoding="latin-1") if l.lower().startswith("flow unit"))
    print("%-14s %6d values  %6d differ  max rel diff %.3e  %s%s" % (
        base, nval, nbad, worst, worst_key, "  (SI units: diameters in mm/1000 m, tolerance 1e-10)" if si else ""))
    return 0.0 if (si and worst < 1e-10) else worst


def main(argv):
    exe = os.path.join(ROOT, "bin", "swmm5_is")
    verbose = "-v" in argv
    args = [a for a in argv if a != "-v"]
    if "--bin" in args:
        i = args.index("--bin")
        exe = args[i + 1]
        del args[i:i + 2]
    txts = []
    for a in args:
        txts += sorted(glob.glob(os.path.join(a, "*.txt"))) if os.path.isdir(a) else [a]
    work = os.path.join(HERE, "work_steady")
    os.makedirs(work, exist_ok=True)
    worst = 0.0
    n = 0
    for t in txts:
        d, b = os.path.split(t)
        b = os.path.splitext(b)[0]
        ref = next((p for p in (os.path.join(d, b + ".out"), os.path.join(d, "ref_" + b + ".out"),
                                 os.path.join(d, "ref_" + b.split("_")[-1] + ".out"))
                    if os.path.exists(p)), None)
        if not ref:
            continue
        w = check(t, ref, exe, work, verbose)
        if w is None:
            worst = float("inf")
        else:
            worst = max(worst, w)
        n += 1
    print("%d models, overall max relative difference %.3e" % (n, worst))
    return 0 if worst == 0.0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
