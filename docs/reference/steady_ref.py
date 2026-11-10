#!/usr/bin/env python3
"""steady_ref.py - standalone re-implementation of the InfoSewer engine
(Innovyze/MWH Soft HeEngPkg) STEADY-STATE run for gravity-only sanitary
networks: loading manholes, circular gravity pipes, outlets.  It follows
/home/claude/ise/spec/steady.md, which cites the C++ (file:line) for every
rule; the function names below are the C++ ones.

Usage:
  python3 -I steady_ref.py <input.txt>                      print results
  python3 -I steady_ref.py <input.txt> <engine.out>         compare with the engine .out
  python3 -I steady_ref.py <input.txt> <engine.out> -v      ... and list every value

Scope: no junction chambers, wet wells, pumps, force mains, channels
(non-circular shapes), storm (rational method) loads, design mode or water
quality.  The script stops if the input uses any of these.

Written by Claude (AI) for SWMM Enablement, 2026-11-09, from a reading of the
engine source; it reproduces the engine's quirks on purpose (see steady.md 5).
"""
import sys, math, struct

# ---------------------------------------------------------------- constants
# he_define.h:93-108, stdafx.h:124-125, he_define.h:35, he_define.h:40-51
HMAP_ZeroTolerance        = 0.0001
HMAP_DepthTolerance       = 0.0000001
HMAP_IterFlowTolerance    = 0.00001
HMAP_MaxManningCoeff      = 1.0
HMAP_WetAreaTolerance     = 0.00001
HMAP_DepthDiamTolerance   = 0.01
HMAP_FlowTolerance        = 0.0001
SYS_Zero                  = 1.0e-10
SYS_ZeroFloat             = 1.0e-7
PI_DEF                    = 3.1415926          # 'pi' macro, used by infiltration BySurface
FLOW_UNITS = ["CFS", "GPM", "MGD", "IMGD", "AFD", "LPS", "LPM", "MLD", "CMH", "CMD"]
FLOW_FACTOR = {"CFS": 1.0, "GPM": 448.831169, "MGD": 0.64632, "IMGD": 0.5382, "AFD": 1.9837,
               "LPS": 28.316847, "LPM": 1699.010796, "MLD": 2.446575, "CMH": 101.940648,
               "CMD": 2446.575546}
MperFT = 0.3048
BIG = 1.0e37

LOAD_UNPEAK, LOAD_PEAK, LOAD_COVER = 0, 1, 2
SPLIT_HYD, SPLIT_PCT, SPLIT_CURVE, SPLIT_INOUT = 0, 1, 2, 3
INF_NONE, INF_LEN, INF_LENDIA, INF_SURF, INF_COUNT, INF_PATTERN = 0, 1, 2, 3, 4, 5
CURVE_PEAKFLOW, CURVE_PEAKPOP = 4, 7


class Model:
    pass


# ---------------------------------------------------------------- input
def toks(line):
    # Fv_GetTokenData he_input.cpp:13163: ';' starts a comment; separators " \t\n\r"
    line = line.split(";", 1)[0]
    for c in "\t\r\n":
        line = line.replace(c, " ")
    return [t for t in line.split(" ") if t]


def getdouble(s):
    # Fv_GetDouble he_input.cpp:13233: first ',' becomes '.', strtod must consume all
    s = s.replace(",", ".", 1)
    return float(s)


def read_sections(path):
    sec, cur = {}, None
    with open(path, encoding="latin-1") as f:
        for ln in f:
            s = ln.rstrip("\n")
            if s.startswith("["):
                cur = s.strip().lower()
                sec[cur] = []
            elif cur is not None:
                sec[cur].append(s)
    return sec


def read_model(path):
    sec = read_sections(path)
    for name in ("[junction chambers]", "[wet wells]", "[force mains]", "[pumps]", "[channel shape]"):
        if any(toks(l) for l in sec.get(name, [])):
            sys.exit("steady_ref.py: %s present - outside this reference's scope" % name)
    M = Model()
    # defaults he_input.cpp:457-535
    M.flowUnit = "CFS"
    M.loadMult = 1.0
    M.K, M.P = 2.4, 0.89
    M.A, M.B, M.C, M.D, M.E = 5.0, 0.0, 0.2, 0.0, 1.0
    M.peakCurveID = None
    M.popCurveID = None
    M.flowTol = HMAP_FlowTolerance
    M.iterTol = HMAP_IterFlowTolerance
    M.headTypeOpt = 0                      # EN_HeadType_Unlocked
    M.storm = 0
    # [Options] Fv_ReadOption he_input.cpp:9583
    for ln in sec.get("[options]", []):
        t = toks(ln)
        if not t:
            continue
        item = t[-1]
        itemIndex = len(t) - 1
        try:
            val = getdouble(item); ok = True
        except ValueError:
            val = 0.0; ok = False
        ps = {i: x for i, x in enumerate(t) if x != ":"}
        p = lambda i: ps.get(i, "").lower()
        if p(0) == "manhole" and p(1) == "head" and p(2) == "type":
            o = int(float(item)) if ok else 0
            if o in (0, 1):
                M.headTypeOpt = o
        elif p(0) == "storm" and p(1) == "runoff" and p(2) == "analysis":
            M.storm = int(val) if val >= 0 else 0
        elif t[0].lower().startswith("global"):
            M.loadMult = val
        elif p(0) == "peaking":
            c = p(2)[:1]
            if c == "k":
                M.K = val
            elif c == "p":
                M.P = val
            elif p(2) == "curve":
                if itemIndex > 3:
                    M.peakCurveID = item
            elif c == "a":
                M.A = val
            elif c == "b":
                M.B = val
            elif c == "c":
                M.C = val
            elif c == "d":
                M.D = val
            elif c == "e":
                M.E = val
        elif p(0) == "flow" and p(1) == "unit":
            M.flowUnit = item.upper()
        elif p(0) == "flow" and p(1) == "tolerance":
            M.flowTol = val
        elif p(0) == "coverage" and p(1) == "based":
            M.popCurveID = item
    if M.storm == 1:
        sys.exit("steady_ref.py: Storm Runoff Analysis = 1 (rational method) is outside scope")
    fF = FLOW_FACTOR[M.flowUnit]
    us = FLOW_UNITS.index(M.flowUnit) <= FLOW_UNITS.index("AFD")
    lenF, diaF = (1.0, 12.0) if us else (MperFT, 304.8)
    M.fF, M.lenF, M.diaF = fF, lenF, diaF

    # [Curves] Fv_ReadCurve he_input.cpp:1648
    M.curves = {}
    cid = None
    for ln in sec.get("[curves]", []):
        t = toks(ln)
        if not t:
            continue
        if len(t) == 4:
            cid = t[0]
            M.curves[cid] = [int(t[1]), [getdouble(t[2])], [getdouble(t[3])]]
        else:
            M.curves[cid][1].append(getdouble(t[0]))
            M.curves[cid][2].append(getdouble(t[1]))

    # [Loading Manholes] Fv_ReadLoadNode he_input.cpp:3733
    M.nodes = []          # node dicts, index = node index
    M.nidx = {}
    for ln in sec.get("[loading manholes]", []):
        t = toks(ln)
        if not t:
            continue
        if len(t) < 58:
            sys.exit("loading manhole line with %d < 58 fields" % len(t))
        n = dict(id=t[0], kind="load", rim=getdouble(t[1]), dia=getdouble(t[2]), K=getdouble(t[3]),
                 loads=[], bottomRef=BIG, headType=-1)
        for i in range(10):
            i0 = 4 * i + 18
            if t[i0] == "-":
                continue
            base = float(t[i0]) if _isnum(t[i0]) else _atof(t[i0])
            if math.fabs(base) < HMAP_ZeroTolerance:        # he_input.cpp:4252 (input units!)
                continue
            ltype = int(_atof(t[i0 + 1]))
            cover = _atof(t[i0 + 3]) if ltype == LOAD_COVER else 0.0
            n["loads"].append([base, ltype, cover])
        M.nidx[n["id"]] = len(M.nodes)
        M.nodes.append(n)
    M.nLoad = len(M.nodes)
    # [Outlets] Fv_ReadOutlet he_input.cpp:4516
    for ln in sec.get("[outlets]", []):
        t = toks(ln)
        if not t:
            continue
        n = dict(id=t[0], kind="outlet", rim=getdouble(t[1]), dia=getdouble(t[2]), submerge=-BIG)
        M.nidx[n["id"]] = len(M.nodes)
        M.nodes.append(n)
    # [Pipes] Fv_ReadPipe he_input.cpp:4985
    M.pipes = []
    M.pidx = {}
    for ln in sec.get("[pipes]", []):
        t = toks(ln)
        if not t:
            continue
        if len(t) != 9:
            sys.exit("pipe line with %d != 9 fields" % len(t))
        fi, ti, L, d, n_ = (getdouble(x) for x in t[3:8])
        p = dict(id=t[0], frm=M.nidx[t[1]], to=M.nidx[t[2]], fromInv=fi, toInv=ti, length=L,
                 dia=d, n=n_, count=int(_atof(t[8])) + 1, split=SPLIT_HYD, splitPct=0.0,
                 infType=INF_NONE, infRate=0.0, infCount=0.0)
        # Fv_UpdateSlopeData he_data.cpp:5600 - computed in INPUT units
        p["slope"] = (fi - ti) / L
        if math.fabs(p["slope"]) < SYS_ZeroFloat:
            p["slopeType"] = "flat"
        elif p["slope"] > 0.0:
            p["slopeType"] = "standard"
        else:
            p["slopeType"] = "adverse"
        p["adverse"] = fi < ti                         # Fv_CheckPipeInvert he_input.cpp:11800
        M.pidx[p["id"]] = len(M.pipes)
        M.pipes.append(p)
        # bottom reference level he_input.cpp:5256-5270, F_SetBottomRefLevel keeps the minimum
        fn, tn = M.nodes[p["frm"]], M.nodes[p["to"]]
        if fn["kind"] == "load":
            fn["bottomRef"] = min(fn["bottomRef"], fi)
        if tn["kind"] == "load":
            tn["bottomRef"] = min(tn["bottomRef"], ti)
    # [Flow Splits] Fv_ReadFlowSplit he_input.cpp:7875
    for ln in sec.get("[flow splits]", []):
        t = toks(ln)
        if len(t) != 4 or (t[2][0] == "-" and t[3][0] == "-") or t[0] not in M.pidx:
            continue
        p = M.pipes[M.pidx[t[0]]]
        st = int(_atof(t[1]))
        if st == SPLIT_PCT and t[2][0] != "-":
            pct = getdouble(t[2])
            if 0.0 <= pct <= 100.0:
                p["splitPct"] = pct
                if pct > 0.0:                        # F_SetSplitFlowPercent he_data.cpp
                    p["split"] = SPLIT_PCT
        elif st in (SPLIT_CURVE, SPLIT_INOUT):
            sys.exit("curve-based flow splits are outside this reference's scope")
    # [Infiltrate Rates] Fv_ReadInfiltrateRate he_input.cpp:8031
    for ln in sec.get("[infiltrate rates]", []):
        t = toks(ln)
        if not t or t[0] not in M.pidx:
            continue
        p = M.pipes[M.pidx[t[0]]]
        p["infType"] = int(_atof(t[1]))
        p["infRate"] = getdouble(t[2])
        if p["infType"] == INF_COUNT:
            p["infCount"] = getdouble(t[3])
        if p["infType"] == INF_PATTERN:
            sys.exit("pattern infiltration is outside this reference's scope")
    # [Extra Loads] Fv_ReadExtraLoad he_input.cpp:8488 (no zero-tolerance filter)
    for ln in sec.get("[extra loads]", []):
        t = toks(ln)
        if len(t) != 5 or t[1] == "-":
            continue
        ni = M.nidx[t[0]]
        ltype = int(_atof(t[2]))
        cover = _atof(t[4]) if ltype == LOAD_COVER else 0.0
        M.nodes[ni]["loads"].append([_atof(t[1]), ltype, cover])
    # [Load Node Head Type] Fv_ReadLoadNodeHeadType he_input.cpp:8676
    for ln in sec.get("[load node head type]", []):
        t = toks(ln)
        if len(t) >= 2 and t[0] in M.nidx and M.nodes[M.nidx[t[0]]]["kind"] == "load":
            ht = int(_atof(t[1]))
            if ht in (0, 1):
                M.nodes[M.nidx[t[0]]]["headType"] = ht
    # [Outlet Heads] Fv_ReadOutletHead he_input.cpp:8742
    for ln in sec.get("[outlet heads]", []):
        t = toks(ln)
        if len(t) >= 2 and t[0] in M.nidx and M.nodes[M.nidx[t[0]]]["kind"] == "outlet":
            M.nodes[M.nidx[t[0]]]["submerge"] = getdouble(t[1])

    # ---- Fv_FillCheckDefaultData he_input.cpp:11131: default head type
    for n in M.nodes:
        if n["kind"] == "load" and n["headType"] < 0:
            n["headType"] = M.headTypeOpt
    # ---- Fv_UnitConversion he_input.cpp:12143
    for n in M.nodes:
        n["rim"] /= lenF
        n["dia"] /= lenF
        if n["kind"] == "load":
            n["bottomRef"] /= lenF
            for ld in n["loads"]:
                ld[0] /= fF                           # HD_LoadData::F_ConvertUnit he_data.cpp:312
        else:
            n["submerge"] /= lenF
    for p in M.pipes:                                 # HD_PipeData::F_ConvertUnit he_data.cpp:5091
        p["fromInv"] /= lenF
        p["toInv"] /= lenF
        p["length"] /= lenF
        p["dia"] /= diaF
        p["infRate"] /= fF
        if p["infType"] == INF_LEN:
            p["infRate"] *= lenF
        elif p["infType"] in (INF_LENDIA, INF_SURF):
            p["infRate"] *= lenF
            p["infRate"] *= diaF
    for cid, c in M.curves.items():                   # he_input.cpp:12269 (peak-flow curves only)
        if c[0] == CURVE_PEAKFLOW:
            s = 1.0 / fF
            c[1] = [x * s for x in c[1]]
            c[2] = [y * s for y in c[2]]
    M.flowTol /= fF                                   # he_input.cpp:12321
    if M.iterTol > M.flowTol:                         # he_input.cpp:12339
        M.iterTol = M.flowTol
    M.peakCurve = M.curves.get(M.peakCurveID) if M.peakCurveID else None
    M.popCurve = M.curves.get(M.popCurveID) if M.popCurveID else None
    return M


def _isnum(s):
    try:
        float(s); return True
    except ValueError:
        return False


def _atof(s):
    # C atof: longest numeric prefix, 0.0 if none
    best = 0.0
    for k in range(len(s), 0, -1):
        try:
            return float(s[:k])
        except ValueError:
            continue
    return best


def curve_y(c, x):
    # HD_CurveData::F_GetInterpolatedY he_data.cpp:10224
    X, Y = c[1], c[2]
    n = len(X)
    if n <= 0 or x < 0.0:
        return 0.0
    last = n - 1
    if x <= X[0]:
        if X[0] == 0.0:
            return 0.0
        return x * Y[0] / X[0]
    elif x >= X[last]:
        if n == 1:
            return Y[0]
        dx = X[last] - X[last - 1]; dy = Y[last] - Y[last - 1]
        if math.fabs(dx) < HMAP_FlowTolerance:
            return Y[last]
        return Y[last] + (x - X[last]) * dy / dx
    if n == 1:
        return Y[0]
    for i in range(1, n):
        if X[i] >= x:
            dx = X[i] - X[i - 1]; dy = Y[i] - Y[i - 1]
            if math.fabs(dx) < HMAP_FlowTolerance:
                return Y[i]
            return Y[i] - (X[i] - x) * dy / dx
    return 0.0


# ---------------------------------------------------------------- topology
def topology(M):
    N, L = len(M.nodes), len(M.pipes)
    # Fv_SetLinkConnection he_topology.cpp:308 - connecting links in link input order
    M.conn = [[] for _ in range(N)]
    for i, p in enumerate(M.pipes):
        M.conn[p["frm"]].append(i)
        M.conn[p["to"]].append(i)
    # Flow-split kernel type: Fv_SetFlowSplitLoadNode he_input.cpp:5680
    for i, n in enumerate(M.nodes):
        nout = sum(1 for l in M.conn[i] if M.pipes[l]["frm"] == i)
        n["split"] = nout > 1
        if n["split"] and n["kind"] != "load":
            sys.exit("flow split at a non-loading node")
    # Fv_ManholeOrder he_topology.cpp:1377
    indeg = [sum(1 for l in M.conn[i] if M.pipes[l]["to"] == i) for i in range(N)]
    vn, vl = [False] * N, [False] * L
    order = []
    while True:
        done = False
        for i in range(N):
            if vn[i]:
                continue
            k = i
            while k >= 0 and indeg[k] == 0:
                done = True
                vn[k] = True
                order.append(k)
                nvisit, nto = 0, -1
                for l in M.conn[k]:
                    if vl[l]:
                        continue
                    nto = M.pipes[l]["to"]
                    if vn[nto]:
                        sys.exit("loop found")
                    vl[l] = True
                    if indeg[nto] > 0:
                        nvisit += 1
                        indeg[nto] -= 1
                if nvisit == 1:
                    k = nto
                else:
                    break
        if len(order) == N or not done:
            break
    if len(order) < N:
        sys.exit("multiple loops")
    M.order = order
    # Fv_LinkOrder he_topology.cpp:1525 (after HU_ResetVisit)
    vl = [False] * L
    M.down = [[] for _ in range(N)]
    M.up = [[] for _ in range(N)]
    for k in order:
        n = M.nodes[k]
        edge = len(M.conn[k])
        if n["kind"] == "outlet":
            continue
        if edge == 1:
            l = M.conn[k][0]
            if not vl[l]:
                M.down[k] = [l]
                vl[l] = True
        elif edge > 1:
            if not n["split"]:
                dl = -1
                for l in M.conn[k]:
                    if not vl[l] and M.pipes[l]["frm"] == k:
                        dl = l
                        break
                M.down[k] = [dl]
                M.up[k] = [l for l in M.conn[k] if l != dl]
                vl[dl] = True
            else:
                for l in M.conn[k]:
                    if not vl[l]:
                        M.down[k].append(l)
                        vl[l] = True
                    else:
                        M.up[k].append(l)
    # outlet upstream links (used only by He_RunHead::Fv_SetUpstreamNode, which uses conn - down)
    M.upHead = [[(M.pipes[l]["frm"], l) for l in M.conn[k] if l not in M.down[k]] for k in range(N)]


# ---------------------------------------------------------------- circular geometry
def wet_area(d, r):
    # HHU_GetWetArea he_hydrun.cpp:2995-2999
    angle = (d - r) / r
    angle2 = (2.0 * 3.1416) - (2.0 * math.acos(angle))
    dp = 2.0 * r
    return 0.125 * (angle2 - math.sin(angle2)) * (dp * dp)


def wet_perimeter(d, r):
    # HHU_GetWetPerimeter he_hydrun.cpp:3945-3947
    angle = (d - r) / r
    angle2 = (2.0 * 3.1416) - (2.0 * math.acos(angle))
    return (0.5 * angle2) * (2.0 * r)


def wet_width(d, r):
    # HHU_GetWetWidth he_hydrun.cpp:2873
    return 2.0 * math.sqrt(d * (2 * r - d))


def get_flow(d, r, slope, n):
    # HHU_GetFlow he_hydrun.cpp:3767-3824
    if slope == 0.0:
        slope = 0.00001
    A = wet_area(d, r)
    P = wet_perimeter(d, r)
    R = A / P
    return (1.49 / n) * A * math.pow(R, 0.6667) * math.sqrt(slope)


def get_full_flow(p):
    # HHU_GetFullFlow he_hydrun.cpp:3588-3621 (per barrel)
    D = p["dia"]; r = 0.5 * p["dia"]; s = p["slope"]
    if s == 0.0:
        s = 0.00001
    if s < 0:
        return -1.0
    return get_flow(D, r, s, p["n"])


def get_link_flow(M, i, d, r=-1.0, n=-1.0):
    # HHU_GetLinkFlow he_hydrun.cpp:3312 (circular branch 3373-3396)
    p = M.pipes[i]
    if p["adverse"]:
        return M.fullFlow[i]
    s = p["slope"]
    if s == 0.0:
        s = 0.00001
    if r < 0.0:
        r = 0.5 * p["dia"]
    if n <= 0.0:
        n = p["n"]
    return get_flow(d, r, s, n)


def get_flow_depth(M, i, q):
    # HHU_GetFlowDepth he_hydrun.cpp:2164 (circular)
    if math.fabs(q) < HMAP_FlowTolerance:
        return 0.0
    p = M.pipes[i]
    D = p["dia"]; n = p["n"]; r = 0.5 * D
    dL, qL = 0.0, 0.0
    dU = D
    qU = get_link_flow(M, i, dU, r, n)
    if qU <= q:
        return D
    depth = 0.5 * (dU + dL)
    qN = get_link_flow(M, i, depth, r, n)
    d1 = q - qL
    d2 = q - qU
    av = q - qN
    while True:
        if av * d2 <= 0.0:
            dL = depth
            qL = get_link_flow(M, i, dL, r, n)
            d1 = q - qL
        elif av * d1 <= 0.0:
            dU = depth
            qU = get_link_flow(M, i, dU, r, n)
            d2 = q - qU
        depth = 0.5 * (dL + dU)
        av = q - get_link_flow(M, i, depth, r, n)
        if not ((math.fabs(av) / q) >= 0.0001):
            break
    if p["slopeType"] != "standard":                  # he_hydrun.cpp:2299-2305
        depth = depth + 0.5 * p["length"] * (-1.0 * p["slope"])
        depth = min(depth, D)
    return depth


def get_critical_depth(M, i, q, D):
    # HHU_GetCriticalDepth he_hydrun.cpp:2586 (circular 2642-2679) - bisection bug kept
    if math.fabs(q) < M.flowTol:
        return 0.0
    r = 0.5 * D
    d0, d1 = 0.0, D
    depth = 0.5 * (d1 - d0)
    ref = 0.0
    A = wet_area(depth, r)
    if A == 0.0:
        return 0.0
    Fr = (q / A) / math.sqrt(A / wet_width(depth, r) * 9.81 * 1 / 0.3048)
    while math.fabs(depth - ref) > HMAP_DepthTolerance:
        ref = depth
        if Fr < 1:
            d1 = depth
            depth = 0.5 * (d1 - d0)
        else:
            d0 = depth
            depth = d0 + 0.5 * (d1 - d0)
        A = wet_area(depth, r)
        if A == 0.0:
            return 0.0
        T = wet_width(depth, r)
        Fr = (q / A) / math.sqrt(A / T * 9.81 * 1 / 0.3048)
    return depth


def pipe_velocity(d, r, q, checkCutoff=False):
    # HHU_GetPipeVelocity he_hydrun.cpp:3138-3148
    A = wet_area(d, r)
    if A < HMAP_WetAreaTolerance:
        return 0.0
    if q < HMAP_ZeroTolerance:
        return 0.0
    if checkCutoff and d / r < 0.02:
        return 0.0
    return q / A


def headloss(p, q):
    # HD_PipeData::F_GetHeadloss he_data.cpp:5480 (circular): uses TOTAL link flow
    if p["n"] < HMAP_MaxManningCoeff:
        varN = 4.66 * p["n"] * p["n"] * p["length"] * q * q
        varD = math.pow(p["dia"], 5.33)
        return varN / varD if math.fabs(varD) > SYS_Zero else 0.0
    fp = math.pow(p["n"], 1.852)
    dp = math.pow(p["dia"], 4.871)
    qp = math.pow(q, 1.852)
    varN = 4.727 * p["length"] * qp
    varD = fp * dp
    return varN / varD if math.fabs(varD) > SYS_Zero else 0.0


# ---------------------------------------------------------------- loads
def infiltration(p):
    # HD_PipeData::F_GetInfiltrateFlow he_data.cpp:5420
    t = p["infType"]
    if t == INF_NONE:
        return 0.0
    if t == INF_LEN:
        f = p["infRate"] * p["length"]
    elif t == INF_LENDIA:
        f = p["infRate"] * p["dia"] * p["length"]
    elif t == INF_SURF:
        f = p["infRate"] * PI_DEF * p["dia"] * p["length"]
    elif t == INF_COUNT:
        f = p["infRate"] * p["infCount"]
    else:
        f = 0.0
    return f * p["count"]


def static_get_loading(M, k):
    # HHU_StaticGetLoading he_hydrun.cpp:984 -> (unpeak, peak, cover, coverCount, infil)
    unpeak = peak = cover = ccount = 0.0
    for base, ltype, cu in M.nodes[k]["loads"]:
        if ltype == LOAD_UNPEAK:
            unpeak += M.loadMult * base
        elif ltype == LOAD_PEAK:
            peak += M.loadMult * base
        elif ltype == LOAD_COVER:
            ccount += 0.001 * cu                      # F_GetPerThousandUnit he_data.h:149
            cover += M.loadMult * base
    infil = 0.0
    for l in M.down[k]:
        infil += infiltration(M.pipes[l])
    unpeak += infil
    return unpeak, peak, cover, ccount, infil


def peak_flow(M, base):
    # HHU_StaticGetPeakFlow he_hydrun.cpp:1423
    if M.peakCurve is None:
        return M.K * math.pow(base, M.P)
    return curve_y(M.peakCurve, base)


def cover_flow(M, cflow, cunit):
    # HHU_StaticGetCoverFlow he_hydrun.cpp:1458
    if cunit < HMAP_ZeroTolerance:
        return 0.0
    if M.popCurve is None:
        varC = math.pow(cunit * M.E, M.C)
        factor = (M.A / (M.B + varC)) + M.D
        return factor * cflow
    return cflow * curve_y(M.popCurve, cunit * 1000)


# ---------------------------------------------------------------- static run
def static_run(M):
    N, L = len(M.nodes), len(M.pipes)
    # H_SetRunData / HH_PipeCapacity he_hydrun.cpp:1546
    M.fullFlow = [get_full_flow(p) * p["count"] for p in M.pipes]
    Q = [0.0] * L; UP = [0.0] * L; PK = [0.0] * L; CV = [0.0] * L; CC = [0.0] * L; INF = [0.0] * L
    visit = [False] * L
    static = [None] * N
    for k in M.order:                                 # Fv_StaticRun he_staticrun.cpp:344
        n = M.nodes[k]
        edge = len(M.conn[k])
        if edge == 1 and n["kind"] == "load":
            # Fv_GetSingleSourceNodeFlow he_staticrun.cpp:423
            l = M.conn[k][0]
            static[k] = static_total(M, k)
            UP[l], PK[l], CV[l], CC[l], INF[l] = static_get_loading(M, k)
            Q[l] = UP[l] + peak_flow(M, PK[l]) + cover_flow(M, CV[l], CC[l]) + 0.0
            visit[l] = True
            continue
        # Fv_GetConnectNodeFlow he_staticrun.cpp:518
        up = pk = cv = cc = inf = 0.0
        if n["kind"] == "load":
            static[k] = static_total(M, k)
            up, pk, cv, cc, inf = static_get_loading(M, k)
        for l in M.conn[k]:                           # HHU_StaticAddPrecedeFlow he_hydrun.cpp:354
            if visit[l]:
                up += UP[l]; pk += PK[l]; cv += CV[l]; cc += CC[l]; inf += INF[l]
        dn = M.down[k]
        inflow = up + peak_flow(M, pk) + cover_flow(M, cv, cc) + 0.0
        if math.fabs(inflow) < M.flowTol:
            for l in dn[1:]:                          # BUG kept: loop starts at i = 1
                Q[l] = UP[l] = PK[l] = CV[l] = INF[l] = 0.0
                visit[l] = True
            continue
        if len(dn) == 1:                              # HHU_StaticSingleDownLink he_hydrun.cpp:396
            l = dn[0]
            Q[l], UP[l], PK[l], CV[l], CC[l], INF[l] = inflow, up, pk, cv, cc, inf
            visit[l] = True
        elif len(dn) > 1:
            split_links(M, k, inflow, Q, visit)       # HHU_StaticFlowSplitLinks he_hydrun.cpp:477
            for l in dn:
                ratio = Q[l] / inflow
                UP[l] = ratio * up; PK[l] = ratio * pk; CV[l] = ratio * cv
                CC[l] = ratio * cc; INF[l] = ratio * inf
                visit[l] = True
    M.Q, M.UP, M.PK, M.CV, M.CC, M.INF, M.static = Q, UP, PK, CV, CC, INF, static

    # Fv_GetPipeHydraulic he_staticrun.cpp:681
    full = [False] * L
    for i, p in enumerate(M.pipes):
        if Q[i] > M.flowTol and Q[i] >= M.fullFlow[i]:
            full[i] = True
    # Fv_GetCriticalDepth he_staticrun.cpp:732
    crit = [0.0] * L
    for i, p in enumerate(M.pipes):
        if math.fabs(Q[i]) < M.flowTol:
            crit[i] = 0.0
        else:
            q = min(Q[i], M.fullFlow[i]) / p["count"]
            crit[i] = get_critical_depth(M, i, q, p["dia"])
    # Fv_GetPipeFlowDepth he_staticrun.cpp:776
    depth = [0.0] * L; vel = [0.0] * L
    for i, p in enumerate(M.pipes):
        q = Q[i] / p["count"]
        depth[i] = get_flow_depth(M, i, q) if not full[i] else p["dia"]
        if depth[i] > 0.0:
            vel[i] = pipe_velocity(depth[i], 0.5 * p["dia"], min(q, Q[i] / p["count"]))
        else:
            vel[i] = 0.0
    # Fv_StaticFilterResult he_staticrun.cpp:878 / HHU_FilterDepthVelocity he_hydrun.cpp:3262
    for i, p in enumerate(M.pipes):
        if not (depth[i] / p["dia"] >= HMAP_DepthDiamTolerance):
            depth[i] = 0.0; vel[i] = 0.0
    M.full, M.crit, M.depth, M.vel = full, crit, depth, vel


def static_total(M, k):
    # HHU_StaticGetTotalLoading he_hydrun.cpp:933
    up, pk, cv, cc, inf = static_get_loading(M, k)
    return up + peak_flow(M, pk) + cover_flow(M, cv, cc) + 0.0


def split_links(M, k, inflow, Q, visit):
    # HHU_StaticFlowSplitLinks he_hydrun.cpp:477 + Fc_GetDownstreamSplitData :589
    hyd, pct = [], []
    for l in M.down[k]:
        if visit[l]:
            continue
        Q[l] = 0.0
        (hyd if M.pipes[l]["split"] == SPLIT_HYD else pct).append(l)
    s = 0.0
    if pct:
        # HHU_CompPercentSplitFlow he_hydrun.cpp:1744 (percent type only)
        w = [M.pipes[l]["splitPct"] for l in pct]
        tot = sum(w)
        ratio = 1.0
        if tot == 0.0:
            ratio = 0.0
        elif tot < 100.0:
            ratio = 1.0 if hyd else 100.0 / tot
        elif tot > 100.0:
            ratio = 100.0 / tot
        if not (tot < 100.0 and hyd):
            w = [x * ratio for x in w]
        for l, x in zip(pct, w):
            Q[l] = 0.01 * x * inflow
            s += Q[l]
            visit[l] = True
    if not hyd:
        return
    if len(hyd) == 1:
        Q[hyd[0]] = inflow - s
        visit[hyd[0]] = True
        return
    multi_hyd_split(M, hyd, inflow - s, Q, visit)


def multi_hyd_split(M, hyd, qh, Q, visit):
    # Fc_StaticMultiHydSplitFlow he_hydrun.cpp:690
    if qh < M.flowTol:
        for l in hyd:
            Q[l] = 0.0; visit[l] = True
        return
    p0 = M.pipes[hyd[0]]
    minInv = p0["fromInv"]; maxCrown = minInv + p0["dia"]
    for l in hyd[1:]:
        p = M.pipes[l]
        minInv = min(minInv, p["fromInv"]); maxCrown = max(maxCrown, p["fromInv"] + p["dia"])
    wl = 0.5 * (minInv + maxCrown)
    out = 0.0
    for l in hyd:
        p = M.pipes[l]
        d = wl - p["fromInv"]
        if d <= 0.0:
            Q[l] = 0.0
        else:
            d = min(d, p["dia"])
            Q[l] = p["count"] * get_link_flow(M, l, d)
        out += Q[l]
    lo, hi = minInv, maxCrown
    if out > qh:
        hi = wl
    if out < qh:
        lo = wl
    wl = 0.5 * (hi + lo)
    while math.fabs(out - qh) / out > M.iterTol:
        out = 0.0
        allFull = True
        for l in hyd:
            p = M.pipes[l]
            D = p["dia"]
            d = wl - p["fromInv"]
            if d <= 0.0:
                Q[l] = 0.0; allFull = False
            else:
                Q[l] = p["count"] * get_link_flow(M, l, min(d, D))
                if d < (D - HMAP_ZeroTolerance):
                    allFull = False
            out += Q[l]
        if allFull:
            tf = sum(M.fullFlow[l] for l in hyd)
            if tf <= 0.0:
                break
            for l in hyd:
                Q[l] = qh * M.fullFlow[l] / tf
            break
        if out > qh:
            hi = wl
        if out < qh:
            lo = wl
        wl = 0.5 * (hi + lo)
    for l in hyd:
        visit[l] = True


# ---------------------------------------------------------------- heads (He_RunHead)
def run_head(M):
    N, L = len(M.nodes), len(M.pipes)
    hyd = [-1.0e37] * N
    out = [0.0] * N            # mv_pOutHead is never initialised in C++; every reached node is set
    level = [0] * N
    adjFull = [False] * L
    partial = [False] * L
    ftype = [0] * L
    splitNodes = [k for k in range(N) if M.nodes[k].get("split")]
    splitOrder = {k: -1 for k in splitNodes}
    counter = [0]
    depth, Q = M.depth, M.Q

    def fully(k):
        return level[k] == len(M.upHead[k])

    def adjust_node(k, h):
        # Fv_AdjustNodeHead he_runhead.cpp:1910 (load-manhole / flow-split branches)
        n = M.nodes[k]
        # Fv_AdjustLoadNodeSurcharge: Dynamic only -> no-op
        # Fv_AdjustLoadNodeDownstream he_runhead.cpp:2152
        dh = -1.0e37
        for l in M.down[k]:
            if Q[l] > HMAP_FlowTolerance:
                dh = depth[l] + M.pipes[l]["fromInv"]
            h = max(h, dh)
        # Fv_AdjustNodeMinorLoss he_runhead.cpp:2258 (velocity uses loop index i as pipe index)
        for i, l in enumerate(M.down[k]):
            if Q[l] > HMAP_FlowTolerance:
                D = M.pipes[l]["dia"]
                v = pipe_velocity(D, 0.5 * D, Q[l])
                h += math.fabs(n["K"]) * ((v * v) / (2.0 * 32.174))
        if n["split"]:
            hyd[k] = max(hyd[k], h)
        else:
            hyd[k] = h
        if n["headType"] == 0:
            out[k] = min(hyd[k], n["rim"])
        else:
            out[k] = hyd[k]
        hyd[k] = out[k]

    def adjust_pipe(l, dh, uh):
        # Fv_AdjustPipeHead he_runhead.cpp:1591
        p = M.pipes[l]
        if Q[l] >= (M.fullFlow[l] - HMAP_FlowTolerance):
            return uh
        L_, D, s, ti = p["length"], p["dia"], p["slope"], p["toInv"]
        downRef = ti + D
        if Q[l] > HMAP_FlowTolerance:
            if dh > uh and dh < downRef:
                upH = dh + headloss(p, Q[l])
                if upH > uh:
                    uh = upH
                    partial[l] = True
                    if uh >= p["fromInv"] + D:
                        adjFull[l] = True
            if dh <= downRef + HMAP_DepthTolerance:
                return uh
            deltaLevel = dh - (ti + depth[l])
            angle = math.atan(s)
            refLength = deltaLevel / math.sin(angle)
            dH = 0.0
            if refLength > L_:
                dH = headloss(p, Q[l])
            crossing = dh + dH
            if dh > ti:
                partial[l] = True
            if crossing > uh:
                uh = crossing
                partial[l] = True
                if uh >= p["fromInv"] + D:
                    adjFull[l] = True
        return uh

    def pipe_head(l, uk, dh):
        # Fv_GetPipeHead he_runhead.cpp:1404 (circular)
        p = M.pipes[l]
        if p["slopeType"] == "standard":
            if Q[l] < M.fullFlow[l]:
                if Q[l] <= HMAP_FlowTolerance:        # HM_ZeroLinkFlow_ByInvert (fixed default)
                    uh = max(dh, M.nodes[uk]["bottomRef"])
                else:
                    uh = depth[l] + p["fromInv"]
                uh = adjust_pipe(l, dh, uh)
            else:
                dH = headloss(p, Q[l])
                ref = p["dia"] + p["toInv"]
                uh = max(dh, ref) + dH
        else:
            ui, di = p["fromInv"], p["toInv"]
            if Q[l] > HMAP_FlowTolerance:
                D = p["dia"]
                dH = headloss(p, Q[l]) if depth[l] >= D else 0.0
                ref = depth[l] + di
                dd = max(dh, ref)
                uh = max(ui + depth[l], dd + dH)
                dh = dd
                uh = adjust_pipe(l, dh, uh)
            else:
                uh = ui
                dh = max(dh, di)
                if dh > di:
                    partial[l] = True
                if dh > di + p["dia"]:
                    adjFull[l] = True
        adjust_node(uk, uh)
        uh = hyd[uk]
        ftype[l] = 1 if uh >= (p["fromInv"] + p["dia"] + HMAP_DepthTolerance) else 0

    def node_head(uk, dk, l):
        # Fv_GetNodeHead he_runhead.cpp:1353
        pipe_head(l, uk, hyd[dk])

    def dfs(k):
        # Fv_DepthSearch he_runhead.cpp:1162
        if fully(k):
            return
        if M.nodes[k].get("split"):
            if splitOrder[k] < 0:
                splitOrder[k] = counter[0]; counter[0] += 1
            return
        for (nk, l) in M.upHead[k]:
            node_head(nk, k, l)
            level[k] += 1
            dfs(nk)

    def split_dfs(k):
        # Fv_FlowSplitDepthSearch he_runhead.cpp:1239
        if fully(k):
            return
        for (nk, l) in M.upHead[k]:
            node_head(nk, k, l)
            level[k] += 1
            dfs(nk)

    # F_Run he_runhead.cpp:554
    for k, n in enumerate(M.nodes):
        if n["kind"] != "outlet":
            continue
        # Fv_GetOutletHead he_runhead.cpp:1296
        ff = -1.0e37
        for (nk, l) in M.upHead[k]:
            ff = max(ff, depth[l] + M.pipes[l]["toInv"])
        hyd[k] = max(n["submerge"], ff)
        out[k] = hyd[k]
        dfs(k)
    remain = len(splitNodes)
    start = 0
    while remain > 0:
        tmpStart = counter[0]
        end = counter[0]
        if start >= end:
            break
        for i in range(start, end):
            k = next(s for s in splitNodes if splitOrder[s] == i)
            split_dfs(k)
        remain -= 1
        start = tmpStart
    M.head, M.ftype, M.partial, M.adjFull = out, ftype, partial, adjFull


# ---------------------------------------------------------------- outputs (Fv_PipeData etc.)
def outputs(M):
    P = M.pipes
    L = len(P)
    R = {}
    Q = [max(q, 0.0) for q in M.Q]                    # he_staticrun.cpp:1244-1248
    R["flow"] = Q
    R["unpeak"] = [M.UP[i] - M.INF[i] for i in range(L)]
    R["peak"] = list(M.PK)
    R["cover"] = list(M.CV)
    R["infil"] = list(M.INF)
    R["storm"] = [0.0] * L
    R["flowtype"] = [1.0 if t == 1 else -1.0 for t in M.ftype]
    R["full"] = list(M.fullFlow)
    R["rescap"] = [max(M.fullFlow[i] - Q[i], 0.0) for i in range(L)]
    fr = []
    for i, p in enumerate(P):                         # he_staticrun.cpp:1309-1321
        D = p["dia"]; d = M.depth[i]
        T = wet_width(d, D / 2.0)
        A = wet_area(d, D / 2.0)
        if T <= 0.0 and d > 0.0:
            fr.append(M.vel[i] / math.sqrt(32.2 * d))
        elif A > 0.0:
            with_inf = (A / T) if T != 0.0 else math.inf
            fr.append(M.vel[i] / math.sqrt(32.2 * with_inf))
        else:
            fr.append(0.0)
    R["vel"] = list(M.vel)
    R["dD"] = [M.depth[i] / P[i]["dia"] for i in range(L)]
    R["qQ"] = [Q[i] / get_link_flow(M, i, P[i]["dia"]) for i in range(L)]
    R["depth"] = list(M.depth)
    R["crit"] = list(M.crit)
    R["froude"] = fr
    R["covcnt"] = [1000 * c for c in M.CC]
    bw, ad, av = [], [], []
    for i, p in enumerate(P):                         # he_staticrun.cpp:1429-1460
        uh = M.head[p["frm"]]; dh = M.head[p["to"]]
        dse = M.depth[i] + p["toInv"]
        if dh > dse:
            a = 0.5 * ((dh - p["toInv"]) + (uh - p["fromInv"]))
            a = min(p["dia"], a)
            ad.append(a)
            av.append(pipe_velocity(a, p["dia"] / 2, Q[i], False))
            bw.append(1.0)
        else:
            ad.append(M.depth[i]); av.append(M.vel[i]); bw.append(0.0)
    R["backwater"], R["adjdepth"], R["adjvel"] = bw, ad, av
    R["partial"] = [1.0 if f else -1.0 for f in M.partial]
    loads = [k for k, n in enumerate(M.nodes) if n["kind"] == "load"]
    R["base"] = [sum(ld[0] for ld in M.nodes[k]["loads"]) for k in loads]
    R["static"] = [M.static[k] for k in loads]
    R["nstorm"] = [0.0 for k in loads]
    R["head"] = [M.head[k] for k in loads]
    R["surch"] = [1.0 if M.head[k] >= M.nodes[k]["rim"] else 0.0 for k in loads]
    jump = []
    for k in loads:                                   # he_staticrun.cpp:1326-1369
        j = 0
        if any(fr[l] > 1 for l in M.up[k]):
            if any(fr[l] < 1 for l in M.down[k]):
                j = 1
        jump.append(float(j))
    R["jump"] = jump
    outs = [k for k, n in enumerate(M.nodes) if n["kind"] == "outlet"]
    # Fv_GetOutletFlow he_staticrun.cpp:842 runs inside Fv_StaticRun, i.e. BEFORE the
    # negative-flow clip of Fv_PipeData, so it sums the unclipped link flows
    R["oflow"] = [sum(M.Q[i] for i, p in enumerate(P) if p["to"] == k) for k in outs]
    R["ohead"] = [M.head[k] for k in outs]
    return R


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


def run(path):
    M = read_model(path)
    topology(M)
    static_run(M)
    run_head(M)
    return M, outputs(M)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    M, R = run(sys.argv[1])
    loads = [n["id"] for n in M.nodes if n["kind"] == "load"]
    if len(sys.argv) < 3:
        print("%-24s %14s %10s %10s %8s %10s" % ("pipe", "flow", "depth", "vel", "d/D", "q/Q"))
        for i, p in enumerate(M.pipes):
            print("%-24s %14.8f %10.6f %10.6f %8.5f %10.6f" %
                  (p["id"], R["flow"][i], R["depth"][i], R["vel"][i], R["dD"][i], R["qQ"][i]))
        for i, k in enumerate(loads):
            print("%-24s head %12.6f" % (k, R["head"][i]))
        return
    E = read_out(sys.argv[2], len(M.pipes), len(loads), len(R["ohead"]))
    verbose = "-v" in sys.argv[3:]
    worst_all = 0.0
    print("%-10s %4s %12s %12s %12s  %s" % ("array", "n", "max abs", "max rel", "rel>1e-6", "worst item"))
    for key in PIPE_ARRAYS + NODE_ARRAYS + ["oflow", "ohead"]:
        a, e = R[key], E[key]
        ids = [p["id"] for p in M.pipes] if key in PIPE_ARRAYS else (loads if key in NODE_ARRAYS else
                                                                     [n["id"] for n in M.nodes if n["kind"] == "outlet"])
        mabs, mrel, bad, worst = 0.0, 0.0, 0, ""
        for i, (x, y) in enumerate(zip(a, e)):
            d = abs(x - y)
            r = d / abs(y) if abs(y) > 1e-12 else d
            if d > mabs:
                mabs = d
            if r > mrel:
                mrel, worst = r, ids[i]
            if r > 1e-6:
                bad += 1
            if verbose and r > 1e-9:
                print("   %-8s %-24s ref %.12g engine %.12g" % (key, ids[i], x, y))
        worst_all = max(worst_all, mrel)
        print("%-10s %4d %12.3e %12.3e %12d  %s" % (key, len(a), mabs, mrel, bad, worst))
    print("overall max relative difference: %.3e" % worst_all)


if __name__ == "__main__":
    main()
