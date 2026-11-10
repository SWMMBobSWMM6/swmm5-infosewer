#!/usr/bin/env python3
"""ise2inp.py - turn an InfoSewer engine input file (the text file HeEngPkg reads,
written by iedb2ise.py from an .IEDB) into a SWMM5 .inp for swmm5-infosewer, or for
stock EPA SWMM 5.2.4.

  python3 -I ise2inp.py model.txt model.inp [--mode steady|eps]
                        [--stock KINWAVE|DYNWAVE|STEADY] [--iedb MODEL.IEDB]
                        [--results results.csv]

Without --stock the file is for swmm5-infosewer: FLOW_ROUTING INFOSEWER, the
InfoSewer options, and the InfoSewer-only data in five extra sections:
  [INFOSEWER_LOADS]         Node Base Type(UNPEAK|PEAK|COVER) Pattern|* CoverUnits
  [INFOSEWER_PATTERNS]      Name STEPWISE|LINEAR f1 f2 ...
  [INFOSEWER_CURVES]        Name X Y   (peaking curves)
  [INFOSEWER_SPLITS]        Link Percent
  [INFOSEWER_INFILTRATION]  Link LENGTH|DIAMETER|SURFACE|COUNT|PATTERN Rate [Count|Pattern]
With --stock the same network gets its loads as plain [DWF] and [INFLOWS] so stock
EPA SWMM 5.2.4 can run it (no peaking: steady-state peaking has no stock SWMM5
equivalent).

Network mapping (InfoSewer pipe diameters are inches, or mm with LPS/MLD):
  loading manhole -> [JUNCTIONS]; invert = lowest connected pipe invert,
                     max depth = rim - invert, surcharge depth 1000 if the head
                     type is Locked (InfoSewer does not cap a locked head at the rim)
  manhole minor loss K -> [LOSSES] inlet loss of each conduit leaving the manhole
  outlet          -> [OUTFALLS] FREE, or FIXED at the [Outlet Heads] head
  gravity pipe    -> [CONDUITS] + [XSECTIONS] CIRCULAR, barrels = count + 1,
                     LINK_OFFSETS ELEVATION with the InfoSewer inverts
Numbers are written with Python's shortest round-trip repr, so SWMM reads back
exactly the doubles InfoSewer uses.
With --iedb, NODE.DBF X/Y and VERTEX.DBF give [COORDINATES] and [VERTICES].

Written by Claude (AI) for SWMM Enablement, 2026-11-09.
"""
import sys, os, argparse, collections

SWMM_FLOW_UNITS = {"CFS": "CFS", "GPM": "GPM", "MGD": "MGD", "LPS": "LPS", "MLD": "MLD"}
LOAD_TYPES = {0: "UNPEAK", 1: "PEAK", 2: "COVER"}
INFIL_TYPES = {1: "LENGTH", 2: "DIAMETER", 3: "SURFACE", 4: "COUNT", 5: "PATTERN"}
ZERO_TOL = 0.0001            # HMAP_ZeroTolerance: manhole loads below this are dropped


def num(x):
    """shortest string that reads back as the same double"""
    return repr(float(x))


def toks(line):
    return line.split(";", 1)[0].split()


def read_sections(path):
    secs = collections.OrderedDict()
    cur = None
    for raw in open(path, encoding="latin-1"):
        line = raw.rstrip("\r\n")
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            cur = s[1:-1]
            secs[cur] = []
            continue
        if cur is None or not s:
            continue
        secs[cur].append(line)
    return secs


def kv(lines):
    d = {}
    for l in lines:
        if ":" in l:
            k, v = l.split(":", 1)
            d[k.strip()] = v.strip()
    return d


class Model:
    def __init__(self, path):
        S = self.S = read_sections(path)
        self.title = " ".join(l.strip() for l in S.get("Titles", []))
        self.opt = kv(S.get("Options", []))
        self.time = kv(S.get("Time", []))
        for unsupported in ("Junction Chambers", "Wet Wells", "Force Mains", "Pumps", "Channel Shape"):
            if any(toks(l) for l in S.get(unsupported, [])):
                sys.exit("model has [%s] data: this version handles gravity-only networks" % unsupported)
        if self.opt.get("Storm Runoff Analysis", "0").strip() not in ("0", "0.0", ""):
            sys.exit("Storm Runoff Analysis is on: storm loads are not supported yet")

        self.manholes = collections.OrderedDict()
        for l in S.get("Loading Manholes", []):
            t = toks(l)
            if not t:
                continue
            loads = []
            for i in range(10):
                k = 18 + 4 * i
                if k >= len(t) or t[k] == "-":
                    continue
                base = float(t[k])
                if abs(base) < ZERO_TOL:          # he_input.cpp:4252 (input flow units)
                    continue
                ltype = int(float(t[k + 1]))
                pat = None if t[k + 2] == "-" else t[k + 2]
                cover = float(t[k + 3]) if (ltype == 2 and k + 3 < len(t)) else 0.0
                loads.append((base, ltype, pat, cover))
            self.manholes[t[0]] = dict(rim=float(t[1]), diam=float(t[2]), loss=float(t[3]),
                                       loads=loads, head=None)
        for l in S.get("Extra Loads", []):            # no zero filter here (he_input.cpp:8488)
            t = toks(l)
            if len(t) != 5 or t[1] == "-" or t[0] not in self.manholes:
                continue
            ltype = int(float(t[2]))
            pat = None if t[3] == "-" else t[3]
            cover = float(t[4]) if ltype == 2 else 0.0
            self.manholes[t[0]]["loads"].append((float(t[1]), ltype, pat, cover))
        for l in S.get("Load Node Head Type", []):
            t = toks(l)
            if len(t) >= 2 and t[0] in self.manholes and int(float(t[1])) in (0, 1):
                self.manholes[t[0]]["head"] = int(float(t[1]))
        ht = int(float(self.opt.get("Manhole Head Type", "0") or 0))
        for mh in self.manholes.values():
            if mh["head"] is None:
                mh["head"] = ht if ht in (0, 1) else 0

        self.outlets = collections.OrderedDict()
        for l in S.get("Outlets", []):
            t = toks(l)
            if t:
                self.outlets[t[0]] = dict(rim=float(t[1]), diam=float(t[2]))
        self.outlet_heads = {}
        for l in S.get("Outlet Heads", []):
            t = toks(l)
            if len(t) >= 2 and t[0] in self.outlets:
                self.outlet_heads[t[0]] = float(t[1])

        self.pipes = collections.OrderedDict()
        for l in S.get("Pipes", []):
            t = toks(l)
            if not t:
                continue
            self.pipes[t[0]] = dict(n1=t[1], n2=t[2], z1=float(t[3]), z2=float(t[4]), length=float(t[5]),
                                    diam_in=float(t[6]), n=float(t[7]), extra=int(float(t[8])))

        self.splits = collections.OrderedDict()
        for l in S.get("Flow Splits", []):
            t = toks(l)
            if len(t) != 4 or (t[2][0] == "-" and t[3][0] == "-") or t[0] not in self.pipes:
                continue
            st = int(float(t[1]))
            if st == 1 and t[2][0] != "-":
                pct = float(t[2])
                if 0.0 <= pct <= 100.0:
                    self.splits[t[0]] = pct
            elif st in (2, 3):
                sys.exit("curve-based [Flow Splits] are not supported yet (pipe %s)" % t[0])

        self.infil = collections.OrderedDict()
        for l in S.get("Infiltrate Rates", []):
            t = toks(l)
            if len(t) < 3 or t[0] not in self.pipes:
                continue
            it = int(float(t[1]))
            if it == 0:
                continue
            if it not in INFIL_TYPES:
                sys.exit("infiltration type %d on pipe %s is not supported" % (it, t[0]))
            count = float(t[3]) if (it == 4 and len(t) > 3) else 0.0
            pat = t[4] if (it == 5 and len(t) > 4) else None
            self.infil[t[0]] = (it, float(t[2]), count, pat)

        self.patterns = collections.OrderedDict()   # id -> [type, values]
        for l in S.get("Patterns", []):
            t = toks(l)
            if len(t) < 2:
                continue
            pid, ptype, vals = t[0], int(float(t[1])), [float(x) for x in t[2:]]
            if pid in self.patterns:
                self.patterns[pid][1].extend(vals)
            else:
                self.patterns[pid] = [ptype, vals]

        self.curves = collections.OrderedDict()     # id -> (type, [(x, y)])
        cid = None
        for l in S.get("Curves", []):
            t = toks(l)
            if not t:
                continue
            if len(t) == 4:
                cid = t[0]
                self.curves[cid] = (int(float(t[1])), [(float(t[2]), float(t[3]))])
            elif cid is not None and len(t) >= 2:
                self.curves[cid][1].append((float(t[0]), float(t[1])))

        # node inverts: lowest connected pipe invert
        inv = {}
        for p in self.pipes.values():
            inv[p["n1"]] = min(inv.get(p["n1"], 1e30), p["z1"])
            inv[p["n2"]] = min(inv.get(p["n2"], 1e30), p["z2"])
        self.invert = inv
        for nid in list(self.manholes) + list(self.outlets):
            if nid not in inv:
                sys.exit("node %s has no connected pipe" % nid)

    def hours(self, key, default):
        v = self.time.get(key)
        if v is None or not v.strip():
            return default
        if ":" in v:
            hh, mm, *ss = (float(x) for x in v.split(":"))
            return hh + mm / 60.0 + (ss[0] if ss else 0.0) / 3600.0
        num_, *unit = v.split()
        x = float(num_)
        u = (unit[0].lower() if unit else "hours")
        return x / 60.0 if u.startswith("min") else x * 24.0 if u.startswith("day") else x / 3600.0 if u.startswith("sec") else x

    def curve_option(self, key):
        v = self.opt.get(key, "").strip()
        return v if v in self.curves else None


def hhmm(h):
    s = int(round(h * 3600))
    return "%d:%02d:%02d" % (s // 3600, (s % 3600) // 60, s % 60)


def date_time(h):
    """hours from 01/01/2000 00:00 -> (date, time) strings"""
    d, rem = divmod(h, 24.0)
    return "01/%02d/2000" % (1 + int(d)), hhmm(rem)


def write_inp(m, out, mode, stock, iedb, results):
    W = []
    w = W.append
    fu = m.opt.get("Flow Unit", "CFS").upper()
    if fu not in SWMM_FLOW_UNITS:
        sys.exit("flow unit %s not supported" % fu)
    us = fu in ("CFS", "GPM", "MGD")          # InfoSewer: US lengths in ft, diameters in inches
    start = m.hours("Simulation Start Time", 0.0)
    dur = m.hours("Simulation Duration", 24.0)
    rstep = m.hours("Report Time Step", 1.0)
    pstep = m.hours("Pattern Time Step", 1.0)
    if mode == "steady" and not stock:
        dur, rstep = 1.0, 1.0                     # one period: the solution does not change
    if dur <= 0:
        dur = rstep
    w("[TITLE]\n" + (m.title or "converted from InfoSewer") + "\n")
    w("[OPTIONS]")
    w("FLOW_UNITS           %s" % SWMM_FLOW_UNITS[fu])
    w("INFILTRATION         HORTON")
    w("FLOW_ROUTING         %s" % (stock or "INFOSEWER"))
    w("LINK_OFFSETS         ELEVATION")
    sd, st = date_time(start)
    ed, et = date_time(start + dur)
    w("START_DATE           %s\nSTART_TIME           %s" % (sd, st))
    w("REPORT_START_DATE    %s\nREPORT_START_TIME    %s" % (sd, st))
    w("END_DATE             %s\nEND_TIME             %s" % (ed, et))
    w("DRY_STEP             01:00:00\nWET_STEP             00:15:00")
    w("REPORT_STEP          %s" % hhmm(rstep))
    w("ROUTING_STEP         60")
    w("ALLOW_PONDING        NO")
    if stock == "DYNWAVE":
        w("VARIABLE_STEP        0.75\nMIN_SURFAREA         12.566\nHEAD_TOLERANCE       0.005")
    if not stock:
        a = m.opt.get
        w("INFOSEWER_ANALYSIS   %s" % mode.upper())
        w("INFOSEWER_LOAD_MULT  %s" % a("Global Loading Multiplier", "1.0"))
        for key, name in (("Peaking Factor k", "INFOSEWER_PEAK_K"), ("Peaking Factor p", "INFOSEWER_PEAK_P"),
                          ("Peaking Parameter a", "INFOSEWER_PEAK_A"), ("Peaking Parameter b", "INFOSEWER_PEAK_B"),
                          ("Peaking Parameter c", "INFOSEWER_PEAK_C"), ("Peaking Parameter d", "INFOSEWER_PEAK_D"),
                          ("Peaking Parameter e", "INFOSEWER_PEAK_E"),
                          ("Minimum Pipe Segment", "INFOSEWER_MIN_SEGMENT"),
                          ("Maximum Number of Segments", "INFOSEWER_MAX_SEGMENTS"),
                          ("Minimum Travel Time", "INFOSEWER_MIN_TRAVEL_TIME"),
                          ("Flow Attenuation", "INFOSEWER_FLOW_ATTENUATION"),
                          ("Flow Tolerance", "INFOSEWER_FLOW_TOLERANCE")):
            if a(key):
                w("%-20s %s" % (name, a(key)))
        pk = m.curve_option("Peaking Factor Curve ID")
        cv = m.curve_option("Coverage Based Peaking Curve")
        if pk:
            w("INFOSEWER_PEAK_CURVE %s" % pk)
        if cv:
            w("INFOSEWER_COVER_CURVE %s" % cv)
        w("INFOSEWER_PATTERN_STEP %s" % hhmm(pstep))
        if results:
            w("INFOSEWER_RESULTS_FILE \"%s\"" % results)
    w("")

    w("[JUNCTIONS]\n;;Name           Elevation          MaxDepth           InitDepth  SurDepth   Aponded")
    for nid, mh in m.manholes.items():
        z = m.invert[nid]
        w("%-16s %-18s %-18s 0          %-10s 0" % (nid, num(z), num(max(mh["rim"] - z, 0.0)),
                                                    "1000" if mh["head"] == 1 else "0"))
    w("")
    w("[OUTFALLS]\n;;Name           Elevation          Type       Stage              Gated")
    for nid in m.outlets:
        if nid in m.outlet_heads:
            w("%-16s %-18s FIXED      %-18s NO" % (nid, num(m.invert[nid]), num(m.outlet_heads[nid])))
        else:
            w("%-16s %-18s FREE                          NO" % (nid, num(m.invert[nid])))
    w("")
    w("[CONDUITS]\n;;Name           From             To               Length     Roughness  InOffset           OutOffset          InitFlow   MaxFlow")
    for pid, p in m.pipes.items():
        w("%-16s %-16s %-16s %-10s %-10s %-18s %-18s 0          0" % (
            pid, p["n1"], p["n2"], num(p["length"]), num(p["n"]), num(p["z1"]), num(p["z2"])))
    w("")
    w("[XSECTIONS]\n;;Link           Shape        Geom1                Geom2  Geom3  Geom4  Barrels")
    for pid, p in m.pipes.items():
        geom1 = p["diam_in"] / 12.0 if us else p["diam_in"] / 1000.0     # in -> ft, mm -> m
        w("%-16s CIRCULAR     %-20s 0      0      0      %d" % (pid, num(geom1), p["extra"] + 1))
    w("")
    losses = [(pid, m.manholes[p["n1"]]["loss"]) for pid, p in m.pipes.items()
              if p["n1"] in m.manholes and m.manholes[p["n1"]]["loss"] != 0.0]
    if losses:
        w("[LOSSES]\n;;Link           Kentry     Kexit      Kavg       FlapGate   Seepage")
        for pid, k in losses:
            w("%-16s %-10s 0          0          NO         0" % (pid, num(k)))
        w("")

    mult = float(m.opt.get("Global Loading Multiplier", "1.0"))
    if not stock:
        w("[INFOSEWER_LOADS]\n;;Node           Base               Type    Pattern          CoverUnits")
        for nid, mh in m.manholes.items():
            for base, ltype, pat, cover in mh["loads"]:
                w("%-16s %-18s %-7s %-16s %s" % (nid, num(base), LOAD_TYPES[ltype], pat or "*", num(cover)))
        w("")
        used = {pat for mh in m.manholes.values() for _, _, pat, _ in mh["loads"] if pat}
        used |= {v[3] for v in m.infil.values() if v[3]}
        if used:
            w("[INFOSEWER_PATTERNS]\n;;Name           Type      Factors...")
            for pid, (ptype, vals) in m.patterns.items():
                if pid in used:
                    kind = "LINEAR" if ptype == 1 else "STEPWISE"
                    for i in range(0, len(vals), 6):
                        w("%-16s %-9s %s" % (pid, kind if i == 0 else "", " ".join(num(v) for v in vals[i:i + 6])))
            w("")
        cids = [c for c in (m.curve_option("Peaking Factor Curve ID"),
                            m.curve_option("Coverage Based Peaking Curve")) if c]
        if cids:
            w("[INFOSEWER_CURVES]\n;;Name           X                  Y")
            for cid in dict.fromkeys(cids):
                for x, y in m.curves[cid][1]:
                    w("%-16s %-18s %s" % (cid, num(x), num(y)))
            w("")
        if m.splits:
            w("[INFOSEWER_SPLITS]\n;;Link           Percent")
            for pid, pct in m.splits.items():
                w("%-16s %s" % (pid, num(pct)))
            w("")
        if m.infil:
            w("[INFOSEWER_INFILTRATION]\n;;Link           Type      Rate               Count|Pattern")
            for pid, (it, rate, count, pat) in m.infil.items():
                extra = num(count) if it == 4 else (pat or "") if it == 5 else ""
                w("%-16s %-9s %-18s %s" % (pid, INFIL_TYPES[it], num(rate), extra))
            w("")
    else:
        # stock SWMM5: constant loads as DWF baselines; patterned loads as a step
        # time series per node in [INFLOWS] (base x pattern x multiplier)
        dwf, ts = {}, {}
        for nid, mh in m.manholes.items():
            for base, ltype, pat, cover in mh["loads"]:
                if pat is None:
                    dwf[nid] = dwf.get(nid, 0.0) + base * mult
                else:
                    ts.setdefault(nid, []).append((base * mult, m.patterns[pat][1]))
        if dwf:
            w("[DWF]\n;;Node           Constituent      Baseline")
            for nid, q in dwf.items():
                w("%-16s FLOW             %s" % (nid, num(q)))
            w("")
        if ts:
            w("[INFLOWS]\n;;Node           Constituent      Time Series      Type     Mfactor  Sfactor  Baseline Pattern")
            for nid in ts:
                w("%-16s FLOW             TS_%-13s FLOW     1.0      1.0" % (nid, nid))
            w("")
            w("[TIMESERIES]\n;;Name           Time       Value")
            nsteps = int(round(dur / pstep))
            for nid, items in ts.items():
                for k in range(nsteps + 1):
                    q = sum(b * vals[k % len(vals)] for b, vals in items)
                    t0 = k * pstep
                    w("TS_%-13s %-10s %s" % (nid, hhmm(t0), num(q)))
                    if k < nsteps:   # hold the value to the end of the step (InfoSewer patterns are steps)
                        w("TS_%-13s %-10s %s" % (nid, hhmm(t0 + pstep - 1.0 / 3600.0), num(q)))
            w("")

    w("[REPORT]\nINPUT      NO\nCONTROLS   NO\nSUBCATCHMENTS NONE\nNODES ALL\nLINKS ALL\n")
    if iedb:
        try:
            import dbfread
        except ImportError:
            sys.exit("--iedb needs dbfread")
        ids = set(m.manholes) | set(m.outlets)
        w("[COORDINATES]\n;;Node           X-Coord            Y-Coord")
        for r in dbfread.DBF(os.path.join(iedb, "NODE.DBF")):
            if r["ID"] in ids:
                w("%-16s %-18.3f %-18.3f" % (r["ID"], r["X"], r["Y"]))
        w("")
        vp = os.path.join(iedb, "VERTEX.DBF")
        if os.path.exists(vp):
            w("[VERTICES]\n;;Link           X-Coord            Y-Coord")
            for r in dbfread.DBF(vp):
                if r["ID"] in m.pipes:
                    w("%-16s %-18.3f %-18.3f" % (r["ID"], r["X"], r["Y"]))
            w("")
    open(out, "w").write("\n".join(W) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ise_txt")
    ap.add_argument("out_inp")
    ap.add_argument("--mode", default="steady", choices=["steady", "eps"])
    ap.add_argument("--stock", choices=["KINWAVE", "DYNWAVE", "STEADY"])
    ap.add_argument("--iedb")
    ap.add_argument("--results", help="write INFOSEWER_RESULTS_FILE (full-precision CSV)")
    a = ap.parse_args()
    m = Model(a.ise_txt)
    write_inp(m, a.out_inp, a.mode, a.stock, a.iedb, a.results)
    print("%d junctions, %d outfalls, %d conduits -> %s" % (len(m.manholes), len(m.outlets), len(m.pipes), a.out_inp))


if __name__ == "__main__":
    main()
