#!/usr/bin/env python3
"""gen_variants.py - write InfoSewer steady-state input variants of base.txt
(the 125_21A model) that exercise branches the base model does not:
peakable/coverage loads, global multiplier, infiltration types, percent splits,
minor losses, parallel barrels, locked heads, submerged outlet, flow units,
peaking curves, zero-flow split, heavy surcharge.

Usage: python3 -I gen_variants.py <base.txt> <out folder>
Written by Claude (AI) for SWMM Enablement, 2026-11-09.
"""
import sys, os

base, outdir = sys.argv[1], sys.argv[2]
LINES = open(base, encoding="latin-1").read().split("\n")


def sections(lines):
    idx, cur = {}, None
    for i, l in enumerate(lines):
        if l.startswith("["):
            cur = l.strip()
            idx[cur] = i
    return idx


def edit(name, fn):
    lines = list(LINES)
    fn(lines)
    with open(os.path.join(outdir, name + ".txt"), "w", encoding="latin-1") as f:
        f.write("\n".join(lines))


def set_option(lines, key, value):
    for i, l in enumerate(lines):
        if l.split(":")[0].strip() == key:
            lines[i] = "%s : %s" % (key, value)
            return
    s = sections(lines)["[Options]"]
    lines.insert(s + 1, "%s : %s" % (key, value))


def add_to_section(lines, sec, rows):
    s = sections(lines)[sec]
    for r in reversed(rows):
        lines.insert(s + 1, r)


def lm_rows(lines):
    s = sections(lines)
    a, b = s["[Loading Manholes]"] + 1, s["[Junction Chambers]"]
    return [i for i in range(a, b) if lines[i].strip()]


def pipe_rows(lines):
    s = sections(lines)
    a, b = s["[Pipes]"] + 1, s["[Channel Shape]"]
    return [i for i in range(a, b) if lines[i].strip()]


def map_loads(lines, fn):
    """fn(k, base, type, cover) -> (base, type, cover); k = manhole ordinal"""
    for k, i in enumerate(lm_rows(lines)):
        t = lines[i].split("\t")
        b, ty, cv = float(t[18]), int(t[19]), float(t[21])
        b, ty, cv = fn(k, b, ty, cv)
        t[18], t[19], t[21] = "%.6f" % b, str(ty), "%.6f" % cv
        lines[i] = "\t".join(t)


def set_field(lines, rows, col, fn):
    for k, i in enumerate(rows):
        t = lines[i].split("\t")
        t[col] = fn(k, t[col])
        lines[i] = "\t".join(t)


# a) all loads peakable
edit("peak", lambda L: map_loads(L, lambda k, b, ty, cv: (b, 1, 0.0)))
# b) all loads coverage-type, cover count = 2000 x base
edit("cover", lambda L: map_loads(L, lambda k, b, ty, cv: (b, 2, 2000.0 * b if b else 0.0)))


# c) mixed types, multiplier, non-default peaking parameters
def mixed(L):
    map_loads(L, lambda k, b, ty, cv: (b, k % 3, (500.0 + 37 * k) if k % 3 == 2 else 0.0))
    set_option(L, "Global Loading Multiplier", "1.700000")
    set_option(L, "Peaking Factor k", "1.900000")
    set_option(L, "Peaking Factor p", "0.920000")
    set_option(L, "Peaking Parameter a", "18.000000")
    set_option(L, "Peaking Parameter b", "4.000000")
    set_option(L, "Peaking Parameter c", "0.500000")
    set_option(L, "Peaking Parameter d", "0.300000")
    set_option(L, "Peaking Parameter e", "2.000000")
    add_to_section(L, "[Extra Loads]", ["CTS-500\t0.050000\t1\t-\t0", "CTS-700\t0.020000\t2\t-\t150",
                                        "CTS-300\t0.000050\t0\t-\t0"])
edit("mixed", mixed)


# d) infiltration of every steady type
def infil(L):
    rows = []
    ids = [L[i].split("\t")[0] for i in pipe_rows(L)]
    for k, pid in enumerate(ids):
        ty = 1 + k % 4
        rate = {1: "0.000020", 2: "0.0000015", 3: "0.0000004", 4: "0.003000"}[ty]
        rows.append("%s\t%d\t%s\t%s\t-" % (pid, ty, rate, "3" if ty == 4 else "0"))
    add_to_section(L, "[Infiltrate Rates]", rows)
edit("infil", infil)
# e) percent splits at CTS-460
edit("pct1", lambda L: add_to_section(L, "[Flow Splits]", ["CTS-460:CTS-450\t1\t30.000000\t-"]))
edit("pct2", lambda L: add_to_section(L, "[Flow Splits]", ["CTS-460:CTS-450\t1\t30.000000\t-",
                                                           "CTS-460:CTS-453\t1\t50.000000\t-"]))
edit("pct0", lambda L: add_to_section(L, "[Flow Splits]", ["CTS-460:CTS-450\t1\t0.000000\t-"]))
# f) minor loss coefficient at every manhole
edit("minor", lambda L: set_field(L, lm_rows(L), 3, lambda k, v: "%.6f" % (0.25 + 0.05 * (k % 7))))
# g) parallel barrels on some pipes, incl. both split legs
edit("parallel", lambda L: set_field(L, pipe_rows(L), 8, lambda k, v: str(k % 3)))


# h) locked heads + high submerged outlet + one unlocked node
def locked(L):
    set_option(L, "Manhole Head Type", "1")
    add_to_section(L, "[Outlet Heads]", ["10\t58.000000"])
    add_to_section(L, "[Load Node Head Type]", ["CTS-030\t0"])
edit("locked", locked)
edit("submerged", lambda L: add_to_section(L, "[Outlet Heads]", ["10\t57.500000"]))


# i) flow units (loads converted so internal flows ~ base)
def units(unit, f):
    def fn(L):
        set_option(L, "Flow Unit", unit)
        map_loads(L, lambda k, b, ty, cv: (b * f, ty, cv))
    return fn
edit("mgd", units("MGD", 0.64632))
edit("gpm", units("GPM", 448.831169))


def si(L):
    set_option(L, "Flow Unit", "LPS")
    map_loads(L, lambda k, b, ty, cv: (b * 28.316847, ty, cv))
    set_field(L, lm_rows(L), 1, lambda k, v: "%.6f" % (float(v) * 0.3048))
    set_field(L, lm_rows(L), 2, lambda k, v: "%.6f" % (float(v) * 0.3048))
    for c in (3, 4, 5):
        set_field(L, pipe_rows(L), c, lambda k, v: "%.6f" % (float(v) * 0.3048))
    set_field(L, pipe_rows(L), 6, lambda k, v: "%.6f" % (float(v) * 25.4))
    s = sections(L)["[Outlets]"]
    t = L[s + 1].split("\t")
    t[1], t[2] = "%.6f" % (float(t[1]) * 0.3048), "%.6f" % (float(t[2]) * 0.3048)
    L[s + 1] = "\t".join(t)
    for i, l in enumerate(L):   # criteria curves are in diameter units (inch -> mm)
        if l.startswith("AC ") or l.startswith("DC ") or l.startswith("\t\t"):
            pass
edit("lps", si)


# j) flow-based peaking curve (type 4) and population curve (type 7)
def pcurve(L):
    map_loads(L, lambda k, b, ty, cv: (b, 1 if k % 2 else 2, 0.0 if k % 2 == 1 else 900.0 + k))
    add_to_section(L, "[Curves]", ["PK 4\t0.500000\t1.500000", "\t\t2.000000\t5.000000",
                                   "\t\t10.000000\t18.000000",
                                   "POP 7\t1.000000\t4.000000", "\t\t50.000000\t2.500000",
                                   "\t\t500.000000\t1.800000"])
    set_option(L, "Peaking Factor Curve ID", "PK")
    set_option(L, "Coverage Based Peaking Curve", "POP")
edit("pcurve", pcurve)


# k) heavy loads -> widespread surcharge
edit("heavy", lambda L: map_loads(L, lambda k, b, ty, cv: (b * 6.0 + 0.4, 1 if k % 2 else 0, 0.0)))


# l) zero flow above the split (CTS-460 upward) and light flow elsewhere
def zero(L):
    def fn(k, b, ty, cv):
        return (b, ty, cv)
    rows = lm_rows(L)
    for i in rows:
        t = L[i].split("\t")
        num = int(t[0].split("-")[1])
        if num >= 460:
            t[18] = "0.000000"
        lines_ = t
        L[i] = "\t".join(lines_)
edit("zero", zero)


# m) in-range population curve and a peak curve whose first point is above most flows
def pcurve2(L):
    map_loads(L, lambda k, b, ty, cv: (b, 1 if k % 2 else 2, 0.0 if k % 2 == 1 else 40.0 + k))
    add_to_section(L, "[Curves]", ["PK 4\t3.000000\t6.000000", "\t\t8.000000\t13.000000",
                                   "\t\t30.000000\t40.000000",
                                   "POP 7\t10.000000\t4.000000", "\t\t500.000000\t2.500000",
                                   "\t\t5000000.000000\t1.500000"])
    set_option(L, "Peaking Factor Curve ID", "PK")
    set_option(L, "Coverage Based Peaking Curve", "POP")
edit("pcurve2", pcurve2)


# n) larger flow tolerance
def tol(L):
    set_option(L, "Flow Tolerance", "0.050000")
    map_loads(L, lambda k, b, ty, cv: (b if b > 0.04 or k % 4 else 0.03, ty, cv))
edit("tol", tol)
# o) 100 % to one split leg -> other (hydraulic) leg carries zero
edit("pct100", lambda L: add_to_section(L, "[Flow Splits]", ["CTS-460:CTS-453\t1\t100.000000\t-"]))


# p) second outlet + two-level split (CTS-480 -> CTS-460 bypass, CTS-600 -> OUT2)
def topo(L):
    s = sections(L)
    L.insert(s["[Outlets]"] + 1, "OUT2\t110.000000\t4.000000")
    s = sections(L)
    L.insert(s["[Pipes]"] + 1, "BYP-480\tCTS-480\tCTS-460\t91.400000\t90.700000\t540.000000\t15.000000\t0.013000\t0")
    L.insert(s["[Pipes]"] + 1, "REL-600\tCTS-600\tOUT2\t95.900000\t95.000000\t120.000000\t10.000000\t0.013000\t0")
    add_to_section(L, "[Flow Splits]", ["REL-600\t1\t25.000000\t-"])
    add_to_section(L, "[Outlet Heads]", ["OUT2\t95.600000"])
    map_loads(L, lambda k, b, ty, cv: (b * 1.5, ty, cv))
edit("topo", topo)


# q) split-ordering probe: like topo, but OUT2 listed BEFORE outlet 10 and a high
#    tailwater at 10, so CTS-600's head via CTS-590 arrives after its upstream was traversed
def topo2(L, first):
    topo(L)
    s = sections(L)
    a = s["[Outlets]"] + 1
    rows = [L[a], L[a + 1]]           # OUT2, 10 (topo inserted OUT2 first)
    if first == "10":
        rows = rows[::-1]
    L[a], L[a + 1] = rows
    add_to_section(L, "[Outlet Heads]", ["10\t104.000000"])
edit("topo2a", lambda L: topo2(L, "OUT2"))
edit("topo2b", lambda L: topo2(L, "10"))
