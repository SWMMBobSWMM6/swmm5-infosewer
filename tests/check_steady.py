#!/usr/bin/env python3
"""check_steady.py - run swmm5-infosewer in steady state on InfoSewer engine inputs
and compare every steady-state output array with the InfoSewer engine's own .out.

  python3 -I check_steady.py <dir-or-txt>... [--bin bin/swmm5_is] [-v]

For each X.txt with an engine result X.out beside it (or ref_X.out), the script
converts X.txt with tools/ise2inp.py, runs bin/swmm5_is with
INFOSEWER_RESULTS_FILE, and compares every steady array value by value: 20 per
gravity pipe, 6 per loading manhole, 1 per junction chamber, 2 per outlet, 1 per
wet well, 8 per force main and 5 per pump. Prints the largest relative
difference per model.

Written by Claude (AI) for SWMM Enablement, 2026-11-09.
"""
import sys, os, struct, subprocess, glob, csv

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PIPE_ARRAYS = ["flow", "unpeak", "peak", "cover", "infil", "storm", "flowtype", "full",
               "rescap", "vel", "dD", "qQ", "depth", "crit", "froude", "covcnt",
               "backwater", "adjdepth", "adjvel", "partial"]
NODE_ARRAYS = ["base", "static", "nstorm", "head", "surch", "jump"]
FM_ARRAYS = ["fflow", "funpeak", "fpeak", "fcover", "finfil", "fstorm", "fvel", "floss"]
PUMP_ARRAYS = ["pflow", "phead", "ppower", "pusage", "pspeed"]
# (csv kind, arrays) in the engine's output order after the pipes
GROUPS = [("node", NODE_ARRAYS), ("chamber", ["chead"]), ("outlet", ["oflow", "ohead"]),
          ("wetwell", ["whead"]), ("fm", FM_ARRAYS), ("pump", PUMP_ARRAYS)]


def read_out(path, counts):
    # static .out layout (Fv_StaticOut he_staticrun.cpp:987): header, 20 gravity pipe
    # arrays, then each group of GROUPS; an empty group writes nothing
    b = open(path, "rb").read()
    groups = [("pipe", PIPE_ARRAYS)] + GROUPS
    total = sum(len(arrs) * counts[kind] for kind, arrs in groups)
    off = len(b) - 8 * total
    R = {}
    for kind, arrs in groups:
        n = counts[kind]
        for k in arrs:
            R[k] = list(struct.unpack("<%dd" % n, b[off:off + 8 * n]))
            off += 8 * n
    return R


def read_csv(path):
    rows = {kind: [] for kind in ["pipe"] + [g[0] for g in GROUPS]}
    hdr = None
    for row in csv.reader(open(path)):
        if row[0] == "kind":
            hdr = row
            continue
        d = dict(zip(hdr[2:], (float(x) for x in row[2:])))
        d["id"] = row[1]
        rows[row[0]].append(d)
    R, ids, kind_of = {}, {}, {}
    for kind, arrs in [("pipe", PIPE_ARRAYS)] + GROUPS:
        for k in arrs:
            R[k] = [r[k] for r in rows[kind]]
            kind_of[k] = kind
        ids[kind] = [r["id"] for r in rows[kind]]
    counts = {kind: len(v) for kind, v in rows.items()}
    return R, ids, counts, kind_of


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
    R, idmap, counts, kind_of = read_csv(res)
    E = read_out(ref, counts)
    worst, worst_key = 0.0, ""
    nbad = 0
    for key in R:
        ids = idmap[kind_of[key]]
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


# InfoSewer engine input/run errors and the swmm5-infosewer error that must match
ENGINE_TO_SWMM = {"1441": "152", "3900": "144", "1222": "150"}


def check_fail(txt, expect, exe, work):
    """A model the InfoSewer engine refuses: swmm5-infosewer must stop too."""
    base = os.path.splitext(os.path.basename(txt))[0]
    inp, rpt, out, res = (os.path.join(work, base + e) for e in (".inp", ".rpt", ".out", ".csv"))
    if os.path.exists(res):
        os.unlink(res)
    subprocess.run([sys.executable, "-I", os.path.join(ROOT, "tools", "ise2inp.py"), txt, inp,
                    "--mode", "steady", "--results", res], check=True, capture_output=True)
    subprocess.run([exe, inp, rpt, out], capture_output=True, text=True)
    want = ENGINE_TO_SWMM.get(expect, "?")
    errs = [l.strip() for l in open(rpt, encoding="latin-1") if "ERROR " in l]
    ok = not os.path.exists(res) and any(("ERROR %s:" % want) in e for e in errs)
    print("%-14s engine Error %s -> %s" % (base, expect, (errs[0] if errs else "no error!")[:110]))
    return ok


def unxz(path, work):
    """reference files may be stored xz-compressed (X.txt.xz, X.out.xz): give a
    plain copy in the work folder"""
    if not path.endswith(".xz"):
        return path
    import lzma
    os.makedirs(os.path.join(work, "ref"), exist_ok=True)
    plain = os.path.join(work, "ref", os.path.basename(path)[:-3])
    if not os.path.exists(plain) or os.path.getmtime(plain) < os.path.getmtime(path):
        with lzma.open(path) as f, open(plain, "wb") as g:
            g.write(f.read())
    return plain


def find_ref(d, b, ext=".out"):
    for name in (b + ext, "ref_" + b + ext, "ref_" + b.split("_")[-1] + ext):
        for p in (os.path.join(d, name), os.path.join(d, name + ".xz")):
            if os.path.exists(p):
                return p
    return None


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
        if os.path.isdir(a):
            txts += sorted(glob.glob(os.path.join(a, "*.txt")) + glob.glob(os.path.join(a, "*.txt.xz")))
            txts += sorted(glob.glob(os.path.join(a, "fail", "*.txt")) +
                           glob.glob(os.path.join(a, "fail", "*.txt.xz")))
        else:
            txts.append(a)
    work = os.path.join(HERE, "work_steady")
    os.makedirs(work, exist_ok=True)
    worst = 0.0
    n = 0
    for t in txts:
        d, b = os.path.split(t)
        b = b[:-7] if b.endswith(".txt.xz") else os.path.splitext(b)[0]
        t = unxz(t, work)
        expect = os.path.join(d, b + ".engine_error")
        if os.path.exists(expect):
            if not check_fail(t, open(expect).read().strip(), exe, work):
                worst = float("inf")
            n += 1
            continue
        ref = find_ref(d, b)
        if not ref:
            continue
        w = check(t, unxz(ref, work), exe, work, verbose)
        if w is None:
            worst = float("inf")
        else:
            worst = max(worst, w)
        n += 1
    print("%d models, overall max relative difference %.3e" % (n, worst))
    return 0 if worst == 0.0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
