#!/usr/bin/env python3
"""gen_press.py - write steady-state pressure-network variants of the active
Thornton model (tha_steady.txt) that exercise engine paths Thornton itself
does not: Pump Flow Continuity 1, parallel pumps, pump speeds and an off pump,
standard siphons (top manhole, mid-network, chained through a load manhole and
through a new chamber), multi-force-main station siphons, a lift the pumps
cannot make, and a low wet well.

  python3 -I gen_press.py <tha_steady.txt> <outdir>
Written by Claude (AI) for SWMM Enablement, 2026-11-10.
"""
import sys, os


def load(path):
    secs, order, cur = {}, [], None
    for line in open(path, encoding="latin-1").read().split("\n"):
        if line.startswith("["):
            cur = line.strip()
            order.append(cur)
            secs[cur] = []
        elif cur is not None:
            secs[cur].append(line)
    return secs, order


def save(secs, order, path):
    with open(path, "w", encoding="latin-1", newline="\n") as f:
        f.write("\n".join(line for s in order for line in [s] + secs[s]))


def rows(secs, s):
    return [l for l in secs[s] if l.strip()]


def set_rows(secs, s, rs):
    secs[s] = rs + [""]


def edit(secs, s, rid, col, val):
    out, hit = [], False
    for l in secs[s]:
        t = l.split("\t")
        if l.strip() and t[0] == rid:
            t[col] = val
            hit = True
        out.append("\t".join(t))
    assert hit, (s, rid)
    secs[s] = out


def option(secs, name, val):
    secs["[Options]"] = [("%s : %s" % (name, val)) if l.split(" : ")[0] == name else l
                         for l in secs["[Options]"]]


def pipe_to_fm(secs, pid, c="120.000000"):
    rs = rows(secs, "[Pipes]")
    t = next(r for r in rs if r.split("\t")[0] == pid).split("\t")
    set_rows(secs, "[Pipes]", [r for r in rs if r.split("\t")[0] != pid])
    t[7] = c
    fms = rows(secs, "[Force Mains]") + ["\t".join(t)]
    set_rows(secs, "[Force Mains]", fms)
    return t


def split_fm(secs, fid, nparts, rim_above=8.0):
    """Split force main fid into nparts force mains joined by new chambers."""
    rs = rows(secs, "[Force Mains]")
    i = next(k for k, r in enumerate(rs) if r.split("\t")[0] == fid)
    t = rs[i].split("\t")
    fi, ti, L = float(t[3]), float(t[4]), float(t[5])
    new, ch = [], []
    prev = t[1]
    for p in range(nparts):
        a = fi + (ti - fi) * p / nparts
        b = fi + (ti - fi) * (p + 1) / nparts
        to = t[2] if p == nparts - 1 else "%s_C%d" % (fid, p + 1)
        if p < nparts - 1:
            ch.append("%s\t%.6f\t4.000000\t0.250000" % (to, b + rim_above))
        new.append("\t".join([fid if p == 0 else "%s_%d" % (fid, p + 1), prev, to,
                              "%.6f" % a, "%.6f" % b, "%.6f" % (L / nparts)] + t[6:]))
        prev = to
    set_rows(secs, "[Force Mains]", rs[:i] + new + rs[i + 1:])
    set_rows(secs, "[Junction Chambers]", rows(secs, "[Junction Chambers]") + ch)


def main(src, outdir):
    os.makedirs(outdir, exist_ok=True)
    V = {}

    def variant(name, fn):
        secs, order = load(src)
        fn(secs)
        save(secs, order, os.path.join(outdir, name + ".txt"))
        V[name] = 1

    # pump columns: 0 ID 1 wetwell 2 to 3 type 4 extra 5 cap 6 h0 7 h1 8 q1 9 h2 10 q2 11 speed
    variant("pv1_cont", lambda s: option(s, "Pump Flow Continuity", "1"))

    def parallel(s):
        for pid, n in (("GHC-1", "1"), ("5054", "1"), ("BDC-1", "1"), ("REM_2", "2")):
            edit(s, "[Pumps]", pid, 4, n)
    variant("pv2_parallel", parallel)

    def speed(s):
        for pid, v in (("REM_3", "0.000000"), ("GHC-2", "0.900000"), ("BDC-1", "1.100000")):
            edit(s, "[Pumps]", pid, 11, v)
    variant("pv3_speed", speed)
    # engine stops: Newton and bisection both fail for 5042/RIV-2 (Error 3900)
    variant("xf_pumpflow", lambda s: edit(s, "[Pumps]", "5042", 11, "0.800000"))

    def stdsiphon(s):
        pipe_to_fm(s, "2313")          # top manhole 3643 -> 3613
        pipe_to_fm(s, "1000")          # mid-network 2264 -> 2256
    variant("pv4_stdsiphon", stdsiphon)

    def chain(s):
        pipe_to_fm(s, "678")           # top 1963 -> 1964 (load) ...
        pipe_to_fm(s, "679")           # ... -> next: chained through a load manhole
        pipe_to_fm(s, "862")           # top 2136 -> 2137, then split through a chamber
        split_fm(s, "862", 2)
        split_fm(s, "P-7958", 2)       # station siphon 5054 -> MH-6347 -> 5119
        split_fm(s, "959", 3)          # station siphon BDC-1 -> 2513 -> 1758
    variant("pv5_chain", chain)

    # engine stops: REM pumps cannot make the lift (Error 1441)
    variant("xf_lift", lambda s: edit(s, "[Force Mains]", "1675", 4, "5112.000000"))
    # engine stops: initial level below the minimum level (Error 1222)
    variant("xf_welllevel", lambda s: edit(s, "[Wet Wells]", "9004", 5, "0.000000"))

    def low(s):
        edit(s, "[Wet Wells]", "9004", 3, "0.000000")        # wet well at its bottom
        edit(s, "[Wet Wells]", "9004", 5, "0.000000")
        edit(s, "[Wet Wells]", "9022", 5, "0.500000")
    variant("pv7_lowwell", low)

    def combo(s):
        option(s, "Pump Flow Continuity", "1")
        parallel(s)
        speed(s)
        stdsiphon(s)
    variant("pv8_combo", combo)
    # flow units: same network in MGD (US) and in LPS and MLD (SI: m, mm, m3)
    FF = {"GPM": 448.831169, "MGD": 0.64632, "LPS": 28.316847, "MLD": 2.446575}

    def units(unit):
        def fn(s):
            r = FF[unit] / FF["GPM"]
            si = unit in ("LPS", "MLD")
            m = 0.3048 if si else 1.0
            mm = 25.4 if si else 1.0
            v3 = 0.028317 if si else 1.0
            f = lambda x, k: "%.9f" % (float(x) * k)
            option(s, "Flow Unit", unit)
            option(s, "Flow Tolerance", f(dict(l.split(" : ") for l in rows(s, "[Options]")
                                                if " : " in l)["Flow Tolerance"], r))
            out = []
            for l in rows(s, "[Loading Manholes]"):
                t = l.split("\t")
                t[1], t[2] = f(t[1], m), f(t[2], m)
                for c in range(18, len(t), 4):
                    t[c] = f(t[c], r)
                out.append("\t".join(t))
            set_rows(s, "[Loading Manholes]", out)
            for sec, cols in (("[Junction Chambers]", (1, 2)), ("[Outlets]", (1, 2)),
                              ("[Wet Wells]", (2, 3, 4, 5, 6))):
                out = []
                for l in rows(s, sec):
                    t = l.split("\t")
                    for c in cols:
                        if t[c] != "-":
                            t[c] = f(t[c], m)
                    out.append("\t".join(t))
                set_rows(s, sec, out)
            for sec in ("[Pipes]", "[Force Mains]"):
                out = []
                for l in rows(s, sec):
                    t = l.split("\t")
                    t[3], t[4], t[5], t[6] = f(t[3], m), f(t[4], m), f(t[5], m), f(t[6], mm)
                    out.append("\t".join(t))
                set_rows(s, sec, out)
            out = []
            for l in rows(s, "[Pumps]"):
                t = l.split("\t")
                for c in (5, 8, 10):
                    t[c] = f(t[c], r)
                for c in (6, 7, 9):
                    t[c] = f(t[c], m)
                out.append("\t".join(t))
            set_rows(s, "[Pumps]", out)
            wwc = {l.split("\t")[7] for l in rows(s, "[Wet Wells]")} - {"-"}
            out, cur = [], None
            for l in rows(s, "[Curves]"):
                t = l.split("\t")
                if t[0]:
                    cur = t[0].split(" ")[0]
                if cur in wwc:
                    t[-2], t[-1] = f(t[-2], m), f(t[-1], v3)
                out.append("\t".join(t))
            set_rows(s, "[Curves]", out)
        return fn
    for u in ("MGD", "LPS", "MLD"):
        variant("pv9_" + u.lower(), units(u))
    print("\n".join(sorted(V)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
