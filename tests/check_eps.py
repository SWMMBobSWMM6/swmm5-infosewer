#!/usr/bin/env python3
"""check_eps.py - run swmm5-infosewer in EPS mode on InfoSewer engine inputs and
compare every EPS output array, at every report time, with the InfoSewer engine's
own binary .out, plus the engine's mass balance block in its .rpt.

  python3 -I check_eps.py <txt or dir>... [--bin bin/swmm5_is] [-v]

For each X.txt with an engine result beside it (X.out, ref_X.out or
ref_<last part>.out) the script converts X.txt with tools/ise2inp.py --mode eps,
runs bin/swmm5_is with INFOSEWER_RESULTS_FILE, and compares value by value.

Written by Claude (AI) for SWMM Enablement, 2026-11-09.
"""
import sys, os, struct, subprocess, glob, csv, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PIPE_ARRAYS = ["flow", "vel", "depth", "froude", "dD", "overflow", "backwater",
               "adjdepth", "adjvel", "qQ", "flowtype", "partial"]
NODE_ARRAYS = ["load", "head", "status", "jump"]


def parse_out(path, P, N, O, start_h, nrep):
    """engine EPS .out (layout from spec/eps_ref.py:parse_out)"""
    b = open(path, "rb").read()
    for off in range(0, len(b) - 12, 4):
        n = struct.unpack_from("<i", b, off)[0]
        if n != nrep or off + 4 + 20 * n > len(b):
            continue
        t = struct.unpack_from("<%dd" % n, b, off + 4)
        if t[0] != start_h or not all(t[i + 1] > t[i] for i in range(n - 1)):
            continue
        if b[off + 4 + 8 * n:off + 4 + 8 * n + 2].isdigit():
            break
    else:
        raise SystemExit("time block not found in " + path)
    p = off + 4 + 8 * n + 12 * n
    p += 8 * P                                   # full flow
    steps = []
    for k in range(n):
        d = {}
        for name in ["load", "storm", "totstorm", "head", "status", "jump"]:
            d[name] = struct.unpack_from("<%dd" % N, b, p); p += 8 * N
        d["outflow"] = struct.unpack_from("<%dd" % O, b, p); p += 8 * O
        d["outhead"] = struct.unpack_from("<%dd" % O, b, p); p += 8 * O
        for name in PIPE_ARRAYS:
            d[name] = struct.unpack_from("<%dd" % P, b, p); p += 8 * P
        steps.append(d)
    return steps


def parse_massbal(rpt):
    txt = open(rpt, encoding="latin-1").read()
    def g(label):
        m = re.search(re.escape(label) + r"[^\d\-]*([\-\d\.]+)", txt)
        return float(m.group(1)) if m else None
    return dict(inflow=g("Total external inflow:"), outflow=g("Total extern outflow:"),
                storage=g("Final Gravity Main Storage:"), check=g("Mass Balance Check:"))


def read_csv(path):
    rows = {"pipe": [], "node": [], "outlet": [], "max": [], "massbal": []}
    hdr = None
    for row in csv.reader(open(path)):
        if row[0] == "kind":
            hdr = row
            continue
        rows[row[0]].append(dict(zip(hdr, row)))
    return rows


def start_hours(txt):
    for ln in open(txt, encoding="latin-1"):
        if ln.lower().startswith("simulation start time"):
            v = ln.split(":", 1)[1].split()
            if ":" in v[0]:
                hh, mm, *ss = (float(x) for x in v[0].split(":"))
                return hh + mm / 60.0 + (ss[0] if ss else 0.0) / 3600.0
            x = float(v[0]); u = v[1].lower() if len(v) > 1 else "hour"
            return x / 3600.0 if u.startswith("sec") else x / 60.0 if u.startswith("min") else \
                x * 24.0 if u.startswith("day") else x
    return 0.0


def check(txt, ref, exe, work, verbose):
    base = os.path.splitext(os.path.basename(txt))[0]
    inp, rpt, out, res = (os.path.join(work, base + e) for e in (".inp", ".rpt", ".out", ".csv"))
    if os.path.exists(res):
        os.unlink(res)
    subprocess.run([sys.executable, "-I", os.path.join(ROOT, "tools", "ise2inp.py"), txt, inp,
                    "--mode", "eps", "--results", res], check=True, capture_output=True)
    r = subprocess.run([exe, inp, rpt, out], capture_output=True, text=True)
    if not os.path.exists(res):
        print("%-14s RUN FAILED: %s" % (base, (r.stdout + r.stderr).strip()[-300:]))
        return None
    R = read_csv(res)
    pid = [p["id"] for p in R["pipe"] if p["report"] == "0"]
    nid = [n["id"] for n in R["node"] if n["report"] == "0"]
    oid = [o["id"] for o in R["outlet"] if o["report"] == "0"]
    nrep = len(R["pipe"]) // max(1, len(pid))
    E = parse_out(ref, len(pid), len(nid), len(oid), start_hours(txt), nrep)
    worst, wkey, nbad, nval = 0.0, "", 0, 0
    def cmp(key, ident, rep, x, y):
        nonlocal worst, wkey, nbad, nval
        nval += 1
        d = abs(x - y)
        rel = d / abs(y) if abs(y) > 1e-12 else d
        if rel > 0:
            nbad += 1
            if verbose and nbad <= 40:
                print("   r%-3d %-9s %-22s swmm5_is %.17g engine %.17g" % (rep, key, ident, x, y))
        if rel > worst:
            worst, wkey = rel, "r%d %s %s (%.12g vs %.12g)" % (rep, key, ident, x, y)
    if nrep != len(E):
        print("%-14s report count %d vs engine %d" % (base, nrep, len(E)))
    for rep in range(min(nrep, len(E))):
        P = R["pipe"][rep * len(pid):(rep + 1) * len(pid)]
        Nn = R["node"][rep * len(nid):(rep + 1) * len(nid)]
        Oo = R["outlet"][rep * len(oid):(rep + 1) * len(oid)]
        for key in PIPE_ARRAYS:
            for i, p in enumerate(P):
                cmp(key, p["id"], rep, float(p[key]), E[rep][key][i])
        for key in NODE_ARRAYS:
            for i, n in enumerate(Nn):
                cmp(key, n["id"], rep, float(n[key]), E[rep][key][i])
        for i, o in enumerate(Oo):
            cmp("outflow", o["id"], rep, float(o["outflow"]), E[rep]["outflow"][i])
            cmp("outhead", o["id"], rep, float(o["outhead"]), E[rep]["outhead"][i])
    erpt = os.path.splitext(ref)[0] + ".rpt"
    mb = ""
    if os.path.exists(erpt):
        em = parse_massbal(erpt)
        mm = R["massbal"][0]
        mi, mo, ms = float(mm["in"]), float(mm["out"]), float(mm["storage"])
        ck = abs(100.0 * (mi - (mo + ms)) / mi) if mi else 0.0
        mb = "  massbal in %.3f/%.3f out %.3f/%.3f stor %.3f/%.3f check %.3f/%.3f" % (
            mi, em["inflow"], mo, em["outflow"], ms, em["storage"], ck, em["check"])
    print("%-14s %7d values %6d differ  max rel diff %.3e  %s" % (base, nval, nbad, worst, wkey))
    if mb:
        print(mb)
    return worst


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
    work = os.path.join(HERE, "work_eps")
    os.makedirs(work, exist_ok=True)
    worst, n = 0.0, 0
    for t in txts:
        d, b = os.path.split(t)
        b = os.path.splitext(b)[0]
        ref = next((p for p in (os.path.join(d, b + ".out"), os.path.join(d, "ref_" + b + ".out"),
                                os.path.join(d, "ref_" + b.split("_")[-1] + ".out")) if os.path.exists(p)), None)
        if not ref:
            continue
        w = check(t, ref, exe, work, verbose)
        worst = float("inf") if w is None else max(worst, w)
        n += 1
    print("%d models, overall max relative difference %.3e" % (n, worst))
    return 0 if worst == 0.0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
