#!/usr/bin/env python3
"""eps_ref.py - standalone Python re-implementation of the InfoSewer (HeEngPkg)
EPS / "dynamic" gravity routing (runOption 2), written from the C++ source to
verify the specification in eps.md.  Gravity-only networks: loading manholes,
circular gravity pipes, outlets.  No wet wells, pumps, force mains, siphons,
junction chambers, open channels, storm runoff, water quality.

Usage:
  python3 -I eps_ref.py <input.txt> [--out engine.out] [--trace trace.csv]
                         [--rpt engine.rpt] [--dump ref.csv]

  --out    compare every pipe/node/time value with the engine's binary .out
  --rpt    compare with the 2-decimal [Pipe Results] tables of the engine .rpt
  --dump   write this script's per-report-step pipe results to a CSV

Every formula below carries the he_*.cpp line it reproduces.  Floating-point
operation order follows the C++ so that results are bit-for-bit where possible.

Written by Claude (AI) for SWMM Enablement, 2026-11-09.
"""
import math
import struct
import sys

# ---------------------------------------------------------------- constants
PI = 3.1415926                  # he_define.h:35 (#define pi)
PI_AREA = 3.1416                # literal in HHU_GetWetArea/Perimeter/Width
HMAP_DoubleTimeTolerance = 0.0001   # he_define.h:71
HMAP_FlowTolerance = 0.0001         # he_define.h:108 (compile-time constant)
HMAP_ZeroTolerance = 0.0001         # he_define.h:93
HMAP_DepthTolerance = 0.0000001     # he_define.h:95
HMAP_WetAreaTolerance = 0.00001     # he_define.h:101
HMAP_DepthDiamTolerance = 0.01      # he_define.h:102
HMAP_IterFlowTolerance = 0.00001    # he_define.h:96
ZERO_TIME_VALUE = 0.001             # he_dynamicrun.cpp:30
SYS_ZeroFloat = 1.0e-7              # stdafx.h:125
BIG = 1.0e37

FLOW_UNIT_FACTOR = {   # he_input.cpp:10847-10883 with he_define.h:40-48; internal unit is cfs
    "CFS": 1.0, "GPM": 448.831169, "MGD": 0.64632, "IMGD": 0.5382, "AFD": 1.9837,
}                      # SI units (LPS, LPM, MLD, CMH, CMD) also convert lengths: not covered here


def cint(x):
    """C (int) cast: truncation toward zero."""
    return int(x)


# ---------------------------------------------------------------- input
class Pipe:
    pass


class Node:
    pass


def read_input(path):
    sec, cur = {}, None
    for ln in open(path, encoding="latin-1"):
        s = ln.rstrip("\r\n")
        if s.startswith("["):
            cur = s.strip()
            sec.setdefault(cur, [])
        elif cur is not None and s.strip():
            sec[cur].append(s)
    m = type("Model", (), {})()
    # --- options (he_input.cpp Fv_ReadOption) ---
    opt = {}
    for s in sec.get("[Options]", []):
        if ":" in s:
            k, v = s.split(":", 1)
            opt[k.strip().lower()] = v.strip()
    def fopt(name, dflt):
        v = opt.get(name.lower())
        try:
            return float(v) if v is not None else dflt
        except ValueError:
            return dflt
    m.flowUnit = opt.get("flow unit", "CFS").upper()
    m.fuf = FLOW_UNIT_FACTOR[m.flowUnit]
    m.loadMult = fopt("Global Loading Multiplier", 1.0)
    m.flowTol = fopt("Flow Tolerance", HMAP_FlowTolerance) / m.fuf   # Fv_UnitConversion
    m.iterTol = min(HMAP_IterFlowTolerance, m.flowTol)               # he_input.cpp:12340
    m.attenuation = fopt("Flow Attenuation", 0.0) >= HMAP_ZeroTolerance
    m.minDeltaX = fopt("Minimum Pipe Segment", 100.0)                # ft (US units)
    v = fopt("Maximum Number of Segments", 100.0)
    m.maxSeg = int(v) if v >= 1.0 else 1000                          # he_input.cpp:10597
    v = fopt("Minimum Travel Time", 1.0)
    m.minTravel = v if v >= 1.0 else 1.0                             # he_input.cpp:10554
    m.headLocked = int(fopt("Manhole Head Type", 0.0)) == 1
    if fopt("Storm Runoff Analysis", 0) == 1:
        raise SystemExit("storm runoff analysis not covered by this reference")
    # --- time (he_input.cpp Fv_ReadTime; Fv_GetDecimalTime -> seconds) ---
    def tval(s):                                  # Fv_GetDecimalTime: number + unit prefix
        p = s.split(":", 1)[1].split()
        x = [float(v) for v in p[0].split(":")]
        u = p[1].lower() if len(p) > 1 else "hour"
        f = 1.0 if u.startswith("sec") else 60.0 if u.startswith("min") else \
            86400.0 if u.startswith("day") else 3600.0
        return x[0] * f
    t = {}
    for s in sec.get("[Time]", []):
        if ":" in s:
            k = s.split(":")[0].strip().lower()
            t[k] = tval(s)
    m.startTime = t.get("simulation start time", 0.0)
    m.simDuration = t.get("simulation duration", 24 * 3600.0)
    m.patternTimeStep = t.get("pattern time step", 3600.0)
    m.reportTimeStep = t.get("report time step", 3600.0)
    # he_input.cpp:957-985
    m.patternTimeInc = 1
    if m.patternTimeStep > m.reportTimeStep + HMAP_DoubleTimeTolerance:
        m.patternTimeInc = max(1, cint(m.patternTimeStep / m.reportTimeStep + 0.9995))
        m.patternTimeStep /= m.patternTimeInc
    m.endTime = m.startTime + m.simDuration
    m.reportStepCount = cint(m.simDuration / m.reportTimeStep + 1.0 + HMAP_DoubleTimeTolerance)
    m.H = m.patternTimeStep                      # hm_hydraulicDeltaT (he_input.cpp:985)
    m.nHyd = cint(m.simDuration / m.H + 1.0 + HMAP_DoubleTimeTolerance)   # he_dynamicrun.cpp:396
    # --- patterns (Fv_ReadPattern): values repeated patternTimeInc times ---
    m.pat = {}
    m.patType = {}
    for s in sec.get("[Patterns]", []):
        tok = s.split()
        pid = tok[0].upper()
        m.patType[pid] = int(tok[1])
        vals = m.pat.setdefault(pid, [])
        for x in tok[2:]:
            vals.extend([float(x)] * m.patternTimeInc)
    # --- nodes ---
    m.nodes, m.nodeList = {}, []
    for s in sec.get("[Loading Manholes]", []):
        tok = s.split()
        n = Node()
        n.id, n.kind = tok[0], "load"
        n.rim, n.diam, n.lossCoef = float(tok[1]), float(tok[2]), float(tok[3])
        n.loads = []
        for i in range(10):                       # HMAP_MAX_ManholeLoadEntry
            i0 = 4 * i + 18
            if i0 + 3 >= len(tok) or tok[i0] == "-":
                continue
            base = float(tok[i0])                 # atof
            if abs(base) < HMAP_ZeroTolerance:    # he_input.cpp:4252
                continue
            pat = None if tok[i0 + 2] == "-" else tok[i0 + 2].upper()
            n.loads.append([base / m.fuf, pat])   # HD_LoadData::F_ConvertUnit
        n.bottomRef = BIG
        n.headLocked = m.headLocked
        m.nodes[n.id.upper()] = n
        m.nodeList.append(n)
    for s in sec.get("[Outlets]", []):
        tok = s.split()
        n = Node()
        n.id, n.kind = tok[0], "outlet"
        n.rim, n.diam = float(tok[1]), float(tok[2])
        n.submergeHead = -BIG
        n.loads = []
        m.nodes[n.id.upper()] = n
        m.nodeList.append(n)
    for s in sec.get("[Junction Chambers]", []) + sec.get("[Wet Wells]", []) + \
            sec.get("[Force Mains]", []) + sec.get("[Pumps]", []) + sec.get("[Channel Shape]", []):
        raise SystemExit("model has chambers/wet wells/force mains/pumps/channels: not covered")
    for s in sec.get("[Extra Loads]", []):
        tok = s.split()
        n = m.nodes[tok[0].upper()]
        if tok[1] != "-" and abs(float(tok[1])) >= HMAP_ZeroTolerance:
            pat = None if len(tok) < 4 or tok[3] == "-" else tok[3].upper()
            n.loads.append([float(tok[1]) / m.fuf, pat])
    for s in sec.get("[Outlet Heads]", []):
        tok = s.split()
        m.nodes[tok[0].upper()].submergeHead = float(tok[1])
    for s in sec.get("[Load Node Head Type]", []):
        tok = s.split()
        m.nodes[tok[0].upper()].headLocked = int(tok[1]) == 1
    # default load pattern -> flat pattern (he_input.cpp:11124, 11259)
    for n in m.nodeList:
        for ld in n.loads:
            if ld[1] is not None and ld[1] not in m.pat:
                raise SystemExit("undefined pattern " + ld[1])
    # --- pipes (Fv_ReadPipe) ---
    m.pipes = []
    for s in sec.get("[Pipes]", []):
        tok = s.split()
        p = Pipe()
        p.id = tok[0]
        p.frm, p.to = m.nodes[tok[1].upper()], m.nodes[tok[2].upper()]
        p.fromInv, p.toInv = float(tok[3]), float(tok[4])
        p.length = float(tok[5])
        p.diam = float(tok[6]) / 12.0              # inches -> ft (F_ConvertUnit)
        p.n = float(tok[7])
        p.pc = int(tok[8]) + 1                      # F_SetParallelCount(parallelCount+1)
        p.slope = (p.fromInv - p.toInv) / p.length  # Fv_UpdateSlopeData
        if abs(p.slope) < SYS_ZeroFloat:
            p.slopeType = "flat"
        elif p.slope > 0.0:
            p.slopeType = "standard"
        else:
            p.slopeType = "adverse"
        p.adverse = p.fromInv < p.toInv             # he_input.cpp:11797
        p.outletPipe = p.to.kind == "outlet"
        p.infil = (0, 0.0, 0.0, None)
        p.index = len(m.pipes)
        if p.frm.kind == "load":                    # bottom reference level = min connecting invert
            p.frm.bottomRef = min(p.frm.bottomRef, p.fromInv)
        if p.to.kind == "load":
            p.to.bottomRef = min(p.to.bottomRef, p.toInv)
        m.pipes.append(p)
    m.pipeById = {p.id.upper(): p for p in m.pipes}
    for s in sec.get("[Infiltrate Rates]", []):
        tok = s.split()
        p = m.pipeById[tok[0].upper()]
        typ, rate = int(tok[1]), float(tok[2]) / m.fuf    # HD_PipeData::F_ConvertUnit
        if typ in (2, 3):
            rate *= 12.0                                  # x diameterUnitFactor (in -> ft)
        cnt = float(tok[3]) if typ == 4 else 0.0
        pat = tok[4].upper() if typ == 5 else None
        p.infil = (typ, rate, cnt, pat)
    for s in sec.get("[Flow Splits]", []):
        raise SystemExit("[Flow Splits] percentage/curve splits are not covered by this reference")
    # topology
    for n in m.nodeList:
        n.down = [p for p in m.pipes if p.frm is n]
        n.up = [p for p in m.pipes if p.to is n]
        n.flowSplit = n.kind == "load" and len(n.down) > 1     # Fv_SetFlowSplitLoadNode
    for p in m.pipes:
        p.fullFlow = full_flow(p) * p.pc                         # HH_PipeCapacity
    return m


# ---------------------------------------------------------------- patterns
def pattern_value(m, pid, simTime):
    """HU_GetPatternItemValue for loading patterns (he_util.cpp:693)."""
    if pid is None:                      # flat default pattern -> 1.0
        return 1.0
    v = m.pat[pid]
    size = len(v)
    if simTime < 0.0:
        return 0.0
    item = cint(simTime / m.patternTimeStep + HMAP_DoubleTimeTolerance)
    while item >= size:
        item -= size
    if item < 0:
        item = 0
    if m.patType[pid] == 0:              # HM_PatDataType_LoadStepwise
        return v[item]
    # linear (HM_PatDataType_LoadLinear)
    inc, ps = m.patternTimeInc, m.patternTimeStep
    index1 = cint(simTime / (inc * ps)) * inc
    while index1 >= size:
        index1 -= size
    index2 = index1 + 1
    if index2 >= size:
        index2 = 0
    value1, value2 = v[index1], v[index2]
    if index2 * ps < simTime:
        ref = min(index1, index2) + size
        while (simTime - ref * ps) > (ps * inc):
            ref += size
        if (size - index1) <= 1:
            ref -= 1
    else:
        ref = index1
    deltaTime = simTime - ref * ps
    deltaValue = value2 - value1
    return value1 + deltaTime * deltaValue / (ps * inc)


def node_load(m, n, t):
    """HD_LoadNodeData::F_GetDynamicTimeFlow / F_GetDynamicTotalFlow (he_data.cpp:2071,2119)."""
    q = 0.0
    for base, pid in n.loads:
        q += base * (m.loadMult * pattern_value(m, pid, t))
    return q


def node_load_total(m, n, t):
    q = 0.0
    for base, pid in n.loads:
        q += pattern_value(m, pid, t) * m.loadMult * base
    return q


def infil_flow(m, p, t):
    """HD_PipeData::F_GetInfiltrateFlow (he_data.cpp:5420)."""
    typ, rate, cnt, pat = p.infil
    if typ == 0:
        return 0.0
    if typ == 1:
        f = rate * p.length
    elif typ == 2:
        f = rate * p.diam * p.length
    elif typ == 3:
        f = rate * PI * p.diam * p.length
    elif typ == 4:
        f = rate * cnt
    else:
        f = rate * pattern_value(m, pat, t)
    return f * p.pc


# ---------------------------------------------------------------- hydraulics (he_hydrun.cpp)
def wet_area(depth, radius):           # HHU_GetWetArea circular (2995)
    angle = (depth - radius) / radius
    angle2 = (2.0 * PI_AREA) - (2.0 * math.acos(angle))
    d = 2.0 * radius
    return 0.125 * (angle2 - math.sin(angle2)) * (d * d)


def wet_perimeter(depth, radius):      # HHU_GetWetPerimeter circular (3944)
    angle = (depth - radius) / radius
    angle2 = (2.0 * PI_AREA) - (2.0 * math.acos(angle))
    return (0.5 * angle2) * (2.0 * radius)


def wet_width(depth, radius):          # HHU_GetWetWidth circular
    return 2.0 * math.sqrt(depth * (2 * radius - depth))


def manning_flow(p, depth, radius, slope, n):   # HHU_GetFlow circular (3816)
    if slope == 0.0:
        slope = 0.00001
    a = wet_area(depth, radius)
    per = wet_perimeter(depth, radius)
    r = a / per
    return (1.49 / n) * a * math.pow(r, 0.6667) * math.sqrt(slope)


def full_flow(p):                      # HHU_GetFullFlow (3558)
    slope = p.slope
    if slope == 0.0:
        slope = 0.00001
    if slope < 0:
        return -1.0
    return manning_flow(p, p.diam, 0.5 * p.diam, slope, p.n)


def link_flow(p, depth, radius=-1.0, friction=-1.0):   # HHU_GetLinkFlow (3312)
    if p.adverse:
        return p.fullFlow
    slope = p.slope
    if slope == 0.0:
        slope = 0.00001
    if radius < 0.0:
        radius = 0.5 * p.diam
    if friction <= 0.0:
        friction = p.n
    return manning_flow(p, depth, radius, slope, friction)


def flow_depth(p, actualFlow):         # HHU_GetFlowDepth circular (2164)
    if abs(actualFlow) < HMAP_FlowTolerance:
        return 0.0
    diameter = p.diam
    friction = p.n
    radius = 0.5 * diameter
    dLower, lowerFlow = 0.0, 0.0
    dUpper = diameter
    upperFlow = link_flow(p, dUpper, radius, friction)
    if upperFlow <= actualFlow:
        return diameter
    depth = 0.5 * (dUpper + dLower)
    newFlow = link_flow(p, depth, radius, friction)
    difference1 = actualFlow - lowerFlow
    difference2 = actualFlow - upperFlow
    avdifference = actualFlow - newFlow
    while True:
        if (avdifference * difference2) <= 0.0:
            dLower = depth
            lowerFlow = link_flow(p, dLower, radius, friction)
            difference1 = actualFlow - lowerFlow
        elif (avdifference * difference1) <= 0.0:
            dUpper = depth
            upperFlow = link_flow(p, dUpper, radius, friction)
            difference2 = actualFlow - upperFlow
        depth = 0.5 * (dLower + dUpper)
        avdifference = actualFlow - link_flow(p, depth, radius, friction)
        if not ((abs(avdifference) / actualFlow) >= 0.0001):
            break
    if p.slopeType != "standard":      # he_hydrun.cpp:2299
        depth = depth + 0.5 * p.length * (-1.0 * p.slope)
        depth = min(depth, diameter)
    return depth


def pipe_velocity(p, depth, radius, flow):     # HHU_GetPipeVelocity circular (checkCutoff=false)
    a = wet_area(depth, radius)
    if a < HMAP_WetAreaTolerance:
        return 0.0
    if flow < HMAP_ZeroTolerance:
        return 0.0
    return flow / a


def headloss(p, flow):                 # HD_PipeData::F_GetHeadloss Manning (5506)
    varN = 4.66 * p.n * p.n * p.length * flow * flow
    varD = math.pow(p.diam, 5.33)
    return varN / varD if abs(varD) > 1.0e-10 else 0.0


def is_flat_gravity_main(p):           # HU_IsFlatGravityMain (he_util.cpp:254)
    return not (abs(p.slope) > 1.0e-5)


# ---------------------------------------------------------------- routing
def topo_order(m):
    """Any upstream-first order; the engine uses he_topology's hm_pNodeOrder."""
    done = set()
    order = []
    pending = [n for n in m.nodeList if n.kind != "outlet"]
    while pending:
        nxt = []
        for n in pending:
            if all(p.index in done for p in n.up):
                order.append(n)
                for p in n.down:
                    done.add(p.index)
            else:
                nxt.append(n)
        if len(nxt) == len(pending):
            raise SystemExit("loop in network")
        pending = nxt
    return order


def interp_hydro(m, arr, hydTime):
    """Fv_SinglePrecedeFlow / Fv_SplitPrecedeFlow interpolation (5831-5851)."""
    H = m.H
    item = cint(hydTime / H + ZERO_TIME_VALUE)
    if item > m.nHyd - 2:
        item = m.nHyd - 2
    deltaTime = hydTime - item * H
    if abs(deltaTime) <= ZERO_TIME_VALUE:
        deltaTime = 0.0
    deltaFlow = arr[item + 1] - arr[item]
    return arr[item] + deltaTime * deltaFlow / H


def split_flows(m, n, hydro, mb):
    """Fv_SplitFlow (he_dynamicrun.cpp:2010): hydraulic split at every hydraulic step."""
    pipes = n.down
    out = [[0.0] * (m.nHyd + 1) for _ in pipes]
    simTime = m.startTime
    hydTime = 0.0
    ti = 0
    while True:
        total = node_load_total(m, n, simTime)
        mb["in"] += total * m.H
        for p in n.up:
            total += hydro[p.index][ti]
        inf = 0.0
        for p in pipes:                     # all split pipes are hydraulic here
            inf += infil_flow(m, p, simTime)
        total += inf
        mb["in"] += inf * m.H
        q = multi_hyd_split(m, n, pipes, total)
        for i in range(len(pipes)):
            out[i][ti] = q[i]
        simTime += m.H
        hydTime += m.H
        ti += 1
        if not (simTime <= m.endTime):
            break
    return out


def multi_hyd_split(m, n, pipes, hydCalcFlow):
    """Fv_MultiHydSplitFlow (he_dynamicrun.cpp:2343)."""
    k = len(pipes)
    if k == 1:
        return [hydCalcFlow]
    q = [0.0] * k
    if hydCalcFlow < m.flowTol:
        return q
    capacity = 0.0
    for p in pipes:
        capacity += p.fullFlow
    if hydCalcFlow >= capacity:
        return [(p.fullFlow / capacity) * hydCalcFlow for p in pipes]
    upInv = [p.fromInv for p in pipes]
    minInvert = upInv[0]
    maxCrown = minInvert + pipes[0].diam
    for i in range(1, k):
        minInvert = min(minInvert, upInv[i])
        maxCrown = max(maxCrown, upInv[i] + pipes[i].diam)
    waterLevel = 0.5 * (minInvert + maxCrown)
    outflow = 0.0
    for i, p in enumerate(pipes):
        depth = waterLevel - upInv[i]
        if depth <= 0.0:
            q[i] = 0.0
        else:
            depth = min(depth, p.diam)
            q[i] = p.pc * link_flow(p, depth)
        outflow += q[i]
    minLevel, maxLevel = minInvert, maxCrown
    if outflow > hydCalcFlow:
        maxLevel = waterLevel
    if outflow < hydCalcFlow:
        minLevel = waterLevel
    waterLevel = 0.5 * (maxLevel + minLevel)
    while abs(outflow - hydCalcFlow) / outflow > m.iterTol:
        outflow = 0.0
        allFull = True
        for i, p in enumerate(pipes):
            depth = waterLevel - p.fromInv
            if depth <= 0.0:
                q[i] = 0.0
                allFull = False
            else:
                fd = min(depth, p.diam)
                q[i] = p.pc * link_flow(p, fd)
                if depth < (p.diam - HMAP_ZeroTolerance):
                    allFull = False
            outflow += q[i]
        if allFull:
            return [p.fullFlow for p in pipes]
        if outflow > hydCalcFlow:
            maxLevel = waterLevel
        if outflow < hydCalcFlow:
            minLevel = waterLevel
        waterLevel = 0.5 * (maxLevel + minLevel)
    return q


def pipe_init(m, p):
    """Fv_ManholePipeInit (he_dynamicrun.cpp:4946)."""
    length = p.length
    segCount = cint(length / m.minDeltaX - 0.0001) + 1
    if segCount > m.maxSeg:
        segCount = m.maxSeg
    radius = 0.5 * p.diam
    linkArea = p.pc * PI * radius * radius
    if length > m.minDeltaX:
        deltaX = length / segCount
    else:
        segCount = 1
        deltaX = length
    if is_flat_gravity_main(p):
        dt = m.H * 1.0
    else:
        dt = float(cint((linkArea * deltaX) / p.fullFlow))
    if dt > m.H * 1.0:
        dt = m.H * 1.0
    if dt < m.minTravel:
        dt = m.minTravel
    multiple = float(cint(m.H / dt) + 1)
    dt = m.H / multiple
    if dt < 1.0:
        dt = 1.0
    return segCount, deltaX, dt


def report_data(m, p, qUp):
    """Fv_CalPipeReportData (he_dynamicrun.cpp:10419): returns velocity, depth."""
    if abs(qUp) < m.flowTol:
        return 0.0, 0.0
    pipeFlow = qUp / p.pc
    if qUp < p.fullFlow:
        depth = flow_depth(p, pipeFlow)
    else:
        depth = p.diam
    vel = pipe_velocity(p, depth, 0.5 * p.diam, pipeFlow)
    if depth / p.diam < HMAP_DepthDiamTolerance:      # HHU_FilterDepthVelocity
        depth, vel = 0.0, 0.0
    return vel, depth


def mc_coeffs(m, p, flow0, deltaX):
    """Fv_SegmentFlowAttenuate coefficient block (he_dynamicrun.cpp:6030-6091)."""
    D, n = p.diam, p.n
    flow = 0.001 if flow0 < 0.001 else flow0
    flow = 3.5 * flow
    flow = 0.5 * (flow + 2.5 * flow)
    fullFlow = full_flow(p)
    if fullFlow < m.flowTol:
        fullFlow = flow
    if (flow >= fullFlow) or p.slopeType != "standard":
        flow = 0.25 * fullFlow
        depth = flow_depth(p, flow)   # HHU_GetFlowDepth(pipe, flow, D, n)
        if depth >= D:
            depth = 0.25 * D
            flow = manning_flow(p, depth, 0.5 * D, p.slope, n)
    else:
        depth = flow_depth(p, flow)   # HHU_GetFlowDepth(pipe, flow, D, n)
    radius = 0.5 * D
    area = wet_area(depth, radius)
    top = wet_width(depth, radius)
    celerity = 1.25 * (flow / area)
    unitq = flow / top
    dT = deltaX / (flow / area)
    C = celerity * (dT / deltaX)
    if C > 1.25:
        C = 1.25
        celerity = C / (dT / deltaX)
    elif C < 0.9:
        C = 0.9
        celerity = C / (dT / deltaX)
    Dn = (unitq / p.slope) / (celerity * deltaX)
    c1 = (1 + C - Dn) / (1 + C + Dn)
    c2 = (-1 + C + Dn) / (1 + C + Dn)
    c3 = (1 - C + Dn) / (1 + C + Dn)
    return c1, c2, c3


def route_pipe(m, n, p, downIdx, hydro, splitArr, mb, trace=None):
    """Fv_ManholePipeInit + Fv_TimeLoopAttenuate for one pipe (whole simulation)."""
    segCount, deltaX, dt = pipe_init(m, p)
    p.segCount, p.deltaX, p.dt = segCount, deltaX, dt
    # number of sub-steps K (identical loop in every segment pass)
    times = []
    simTime = m.startTime
    hydTime = 0.0
    while simTime <= (m.endTime + dt):
        times.append((simTime, hydTime))
        simTime = simTime + dt
        hydTime += dt
    K = len(times)
    Former = [0.0] * (K + 3)
    Current = [0.0] * (K + 3)
    Inflow = [0.0] * (K + 3)
    Outflow = [0.0] * (K + 3)
    DF = [0.0] * (m.nHyd + 1)
    DF[0] = 0.0
    hydroIndex = 1
    linkRes = []          # list of (qUp, qDn, vel, depth) per report index
    reportIndex = 0
    rng = {"max": None, "maxTime": None}
    single = len(n.down) == 1
    c1 = c2 = c3 = 0.0
    for seg in range(segCount):
        lf0 = 0.0
        for k in range(1, K + 1):
            simTime, hydTime = times[k - 1]
            nextSimTime = simTime + dt
            if seg == 0:
                if k <= 1:
                    lf0 = 0.0
                lf1 = 0.0
                if single:
                    nextHyd = nextSimTime - m.startTime
                    if k <= 1:
                        lf1 = node_load(m, n, nextSimTime)
                        lf0 = node_load(m, n, simTime)
                        lf0 += infil_flow(m, p, simTime)
                        lf1 += infil_flow(m, p, nextSimTime)
                        mb["in"] += lf0 * dt
                    else:
                        lf1 = node_load(m, n, nextSimTime)
                        lf1 += infil_flow(m, p, nextSimTime)
                        mb["in"] += lf1 * dt
                    # Fv_SinglePrecedeFlow
                    hT = simTime - m.startTime
                    for up in n.up:
                        if k <= 1:
                            lf0 += interp_hydro(m, hydro[up.index], hT)
                        lf1 += interp_hydro(m, hydro[up.index], nextHyd)
                else:
                    arr = splitArr[downIdx]
                    hT = simTime - m.startTime
                    if k <= 1:
                        lf0 = interp_hydro(m, arr, hT)
                    lf1 = interp_hydro(m, arr, nextSimTime - m.startTime)
                if k <= 1:
                    Inflow[k] = lf0
                    Former[k] = Inflow[k]
                Inflow[k + 1] = lf1
                Former[k + 1] = Inflow[k + 1]
            # -- segment outflow --
            if m.attenuation and p.slope > 0.00001:
                if seg <= 0 and k <= 10:
                    c1, c2, c3 = mc_coeffs(m, p, Former[k], deltaX)
                if k <= 1:
                    Current[k] = 0.0
                Current[k + 1] = c1 * Former[k] + c2 * Former[k + 1] + c3 * Current[k]
                if Current[k + 1] <= m.flowTol:
                    Current[k + 1] = 0.0
            elif m.attenuation:
                if k <= 1:
                    Current[k] = 0.0
                Current[k + 1] = Former[k]
                if Current[k + 1] <= m.flowTol:
                    Current[k + 1] = 0.0
            else:                                  # Fv_SegmentOutflow
                Current[k + 1] = Former[k]
                if Current[k + 1] <= m.flowTol:
                    Current[k + 1] = 0.0
            Inflow[k] = Current[k]
            if seg == segCount - 1:
                Outflow[k] = Current[k]
                Outflow[k + 1] = Current[k + 1]
            if seg <= 0:
                Inflow[k] = Former[k]
                Inflow[k + 1] = Former[k + 1]
            if seg == segCount - 1:
                # Fv_DownstreamFlow
                if hydroIndex <= m.nHyd:
                    if hydTime > (hydroIndex * m.H + HMAP_DoubleTimeTolerance):
                        dF = Outflow[k] - DF[hydroIndex - 1]
                        dT = hydTime - (hydroIndex - 1) * m.H
                        DF[hydroIndex] = max(0.0, DF[hydroIndex - 1] + m.H * dF / dT)
                        hydroIndex += 1
                    elif abs(hydTime - hydroIndex * m.H) < HMAP_DoubleTimeTolerance:
                        DF[hydroIndex] = Outflow[k]
                        hydroIndex += 1
                # Fv_ReportPipeData
                if p.outletPipe:
                    if k == 1:
                        mb["out"] += Outflow[k] * dt
                    else:
                        mb["out"] += Outflow[k + 1] * dt
                if k == 1:
                    qUp, qDn = Inflow[k], Outflow[k]
                    linkRes.append((qUp, qDn) + report_data(m, p, qUp))
                    reportIndex += 1
                else:
                    if hydTime >= reportIndex * m.reportTimeStep and reportIndex <= m.reportStepCount:
                        deltaTime = reportIndex * m.reportTimeStep + dt - hydTime
                        if abs(deltaTime) <= ZERO_TIME_VALUE:
                            deltaTime = 0.0
                        qUp = Inflow[k - 1] + deltaTime * (Inflow[k] - Inflow[k - 1]) / dt
                        qDn = Outflow[k - 1] + deltaTime * (Outflow[k] - Outflow[k - 1]) / dt
                        linkRes.append((qUp, qDn) + report_data(m, p, qUp))
                        reportIndex += 1
                # Fv_OutputPipeFlowRange
                q = Current[k]
                if abs(hydTime) < HMAP_DoubleTimeTolerance:
                    rng["max"], rng["maxTime"] = q, simTime
                elif q > rng["max"] + m.flowTol:
                    rng["max"], rng["maxTime"] = q, simTime
                if trace is not None:
                    trace(p, k, hydTime, Former[k] if segCount == 1 else p._in0[k], Current[k])
        if seg == 0:
            p._in0 = list(Former)
        if seg < segCount - 1:
            # shift: Former[k] = Current[k-1] for k = 1..K
            for k in range(1, K + 1):
                Former[k] = Current[k - 1]
    p.linkRes = linkRes
    p.hydro = DF[:m.nHyd]
    p.maxFlow, p.maxTime = rng["max"], rng["maxTime"]
    return DF[:m.nHyd]


def run(m, trace=None):
    mb = {"in": 0.0, "out": 0.0}
    hydro = {}
    for n in topo_order(m):
        splitArr = None
        if len(n.down) > 1 and n.flowSplit:
            splitArr = split_flows(m, n, hydro, mb)
        for i, p in enumerate(n.down):
            hydro[p.index] = route_pipe(m, n, p, i, hydro, splitArr, mb, trace)
    m.mb = mb
    return m


# ---------------------------------------------------------------- report-time outputs and heads
def outputs(m):
    """Per report step: emulate Fv_ReadPipeMainResultFromFile, Fv_DynamicGetHead,
    Fv_PipeFlow/Velocity/Depth/DoverD/QoverQ and node/outlet outputs."""
    res = []
    pipeStorage = 0.0
    for r in range(m.reportStepCount):
        up = [p.linkRes[r][0] for p in m.pipes]
        dn = [p.linkRes[r][1] for p in m.pipes]
        vel = [p.linkRes[r][2] for p in m.pipes]
        dep = [p.linkRes[r][3] for p in m.pipes]
        froude = []
        for i, p in enumerate(m.pipes):                 # Fv_ReadPipeMainResultFromFile
            radius = p.diam / 2.0
            w = wet_width(dep[i], radius)
            a = wet_area(dep[i], radius)
            if (w <= 0.0) and dep[i] > 0.0:
                fr = vel[i] / math.sqrt(32.2 * dep[i])
            elif a > 0.0 and w > 0.0:
                fr = vel[i] / math.sqrt(32.2 * (a / w))
            else:
                fr = 0.0
            froude.append(fr)
            if r == m.reportStepCount - 1 and r > 0:
                pipeStorage += a * p.pc * p.length
        head = node_heads(m, up, dep)
        # Fv_PipeFlow (adverse pipes: sum of upstream pipes' upstream flows)
        flow = []
        upAdj = list(up)
        for i, p in enumerate(m.pipes):
            if p.adverse:
                f = 0.0
                for q in p.frm.up:
                    f += upAdj[q.index]
                upAdj[i] = f
                flow.append(max(0.0, f))
            else:
                flow.append(max(0.0, up[i]))
        velo = []
        for i, p in enumerate(m.pipes):
            if p.adverse:
                d = min(p.diam, dep[i])
                velo.append(pipe_velocity(p, d, 0.5 * p.diam, upAdj[i]))
            else:
                velo.append(vel[i])
        depth = [min(p.diam, dep[i]) for i, p in enumerate(m.pipes)]
        dD = [depth[i] / p.diam for i, p in enumerate(m.pipes)]
        qQ = []
        for i, p in enumerate(m.pipes):
            ff = link_flow(p, p.diam)
            v = max(0.0, flow[i] / ff)
            if dD[i] >= 1.0:
                v = max(v, 1.0)
            qQ.append(v)
        # Fv_PipeBackwaterAdjustment (12760): uses the adverse-modified flows and velocities
        bw, adjD, adjV = [], [], []
        for i, p in enumerate(m.pipes):
            uh, dh = head[p.frm.id], head[p.to.id]
            if dh > dep[i] + p.toInv:
                ad = 0.5 * ((dh - p.toInv) + (uh - p.fromInv))
                ad = min(p.diam, ad)
                ad = max(ad, dep[i])
                adjD.append(ad)
                adjV.append(pipe_velocity(p, p.diam, 0.5 * p.diam, max(0.0, upAdj[i])))
                bw.append(1.0)
            else:
                adjD.append(dep[i])
                adjV.append(velo[i])
                bw.append(0.0)
        status = [1.0 if res_h >= n.rim else 0.0 for n, res_h in
                  ((n, head[n.id]) for n in m.nodeList if n.kind == "load")]
        jump = []                                   # Fv_ReadPipeMainResultFromFile (14747)
        for n in m.nodeList:
            if n.kind != "load":
                continue
            f = 0.0
            if any(froude[q.index] > 1 for q in n.up):
                if any(froude[q.index] < 1 for q in n.down):
                    f = 1.0
            jump.append(f)
        m._jumpHist = getattr(m, "_jumpHist", [])
        loads = [node_load_total(m, n, m.startTime + r * m.reportTimeStep) for n in m.nodeList if n.kind == "load"]
        outflow = {}
        for p in m.pipes:
            if p.outletPipe:
                outflow[p.to.id] = outflow.get(p.to.id, 0.0) + dn[p.index]
        res.append(dict(up=up, dn=dn, flow=flow, vel=velo, depth=depth, dD=dD, qQ=qQ, froude=froude,
                        head=head, load=loads, outflow=outflow, backwater=bw, adjdepth=adjD, adjvel=adjV,
                        status=status, jump=jump,
                        flowtype=[1.0 if t == 1 else -1.0 for t in m._flowType],
                        partial=[1.0 if x else -1.0 for x in m._partial]))
    m.pipeStorage2 = 2.0 * pipeStorage     # accumulated by Fv_DynamicData and again by Fv_LinkReport
    return res


def node_heads(m, linkFlow, pipeDepth):
    """He_RunHead::F_Run for gravity-only networks (he_runhead.cpp:554)."""
    hyd = {n.id: -BIG for n in m.nodeList}
    out = {}
    partial = [False] * len(m.pipes)       # mv_pAdjustedPartial (reset every F_Run)
    flowType = getattr(m, "_flowType", None) or [0] * len(m.pipes)   # mv_pFlowType: never reset
    m._flowType = flowType
    visit = {n.id: 0 for n in m.nodeList}
    splitOrder = {}
    counter = [0]

    def adjust_node(n, h):
        # Fv_AdjustNodeHead (1910): LoadManhole and FlowSplit cases
        # Fv_AdjustLoadNodeSurcharge: surcharge volume is always 0 (V.4) -> no-op
        downHead = -BIG                              # Fv_AdjustLoadNodeDownstream (2152)
        for q in n.down:
            if linkFlow[q.index] > HMAP_FlowTolerance:
                downHead = pipeDepth[q.index] + q.fromInv
            h = max(h, downHead)
        for q in n.down:                             # Fv_AdjustNodeMinorLoss (2258)
            if linkFlow[q.index] > HMAP_FlowTolerance:
                # BUG reproduced: velocity uses loop counter i as pipe index (only area matters
                # for circular pipes of equal diameter); coefficient usually 0
                pass
        if n.lossCoef != 0.0:
            for i, q in enumerate(n.down):
                if linkFlow[q.index] > HMAP_FlowTolerance:
                    pi_ = m.pipes[i]
                    v = pipe_velocity(pi_, q.diam, 0.5 * q.diam, linkFlow[q.index])
                    h += abs(n.lossCoef) * ((v * v) / (2.0 * 32.174))
        if n.flowSplit:
            hyd[n.id] = max(hyd[n.id], h)
        else:
            hyd[n.id] = h
        if not n.headLocked:
            out[n.id] = min(hyd[n.id], n.rim)
        else:
            out[n.id] = hyd[n.id]
        hyd[n.id] = out[n.id]

    def adjust_pipe(p, downHead, upHead):
        # Fv_AdjustPipeHead (1591)
        Q = linkFlow[p.index]
        if Q >= (p.fullFlow - HMAP_FlowTolerance):
            return upHead
        downRef = p.toInv + p.diam
        if Q > HMAP_FlowTolerance:
            if downHead > upHead and downHead < downRef:
                uh = downHead + headloss(p, Q)
                if uh > upHead:
                    upHead = uh
                    partial[p.index] = True
            if downHead <= downRef + HMAP_DepthTolerance:
                return upHead
            deltaLevel = downHead - (p.toInv + pipeDepth[p.index])
            angle = math.atan(p.slope)
            refLength = deltaLevel / math.sin(angle)
            dh = 0.0
            if refLength > p.length:
                dh = headloss(p, Q)
            crossing = downHead + dh
            if downHead > p.toInv:
                partial[p.index] = True
            if crossing > upHead:
                upHead = crossing
                partial[p.index] = True
        return upHead

    def pipe_head(p, upNode, downHead):
        # Fv_GetPipeHead (1404) for circular pipes, hm_optHeadZeroLinkFlow = ByInvert
        Q = linkFlow[p.index]
        if p.slopeType == "standard":
            if Q < p.fullFlow:
                if Q <= HMAP_FlowTolerance:
                    upHead = max(downHead, upNode.bottomRef)
                else:
                    upHead = pipeDepth[p.index] + p.fromInv
                upHead = adjust_pipe(p, downHead, upHead)
            else:
                dh = headloss(p, Q)
                ref = p.diam + p.toInv
                upHead = max(downHead, ref) + dh
        else:
            if Q > HMAP_FlowTolerance:
                dh = 0.0
                if pipeDepth[p.index] >= p.diam:
                    dh = headloss(p, Q)
                ref = pipeDepth[p.index] + p.toInv
                dd = max(downHead, ref)
                upHead = max(p.fromInv + pipeDepth[p.index], dd + dh)
                upHead = adjust_pipe(p, dd, upHead)
            else:
                upHead = p.fromInv
                dd = max(downHead, p.toInv)
                if dd > p.toInv:
                    partial[p.index] = True
        adjust_node(upNode, upHead)
        flowType[p.index] = 1 if hyd[upNode.id] >= (p.fromInv + p.diam + HMAP_DepthTolerance) else 0

    def depth_search(n, splitPhase=False):
        nup = len(n.up)
        if visit[n.id] == nup:
            return
        if (not splitPhase) and n.flowSplit:
            if n.id not in splitOrder:
                splitOrder[n.id] = counter[0]
                counter[0] += 1
            return
        for q in n.up:
            pipe_head(q, q.frm, hyd[n.id])
            visit[n.id] += 1
            depth_search(q.frm)

    for o in m.nodeList:
        if o.kind != "outlet":
            continue
        ff = -BIG
        for q in o.up:
            ff = max(ff, pipeDepth[q.index] + q.toInv)
        hyd[o.id] = max(o.submergeHead, ff)
        out[o.id] = hyd[o.id]
        depth_search(o)
        visit[o.id] = len(o.up)
    remain = sum(1 for n in m.nodeList if n.flowSplit)
    start = 0
    while remain > 0:
        tmpStart = counter[0]
        end = counter[0]
        if start >= end:
            break
        for i in range(start, end):
            nid = [k for k, v in splitOrder.items() if v == i][0]
            depth_search(m.nodes[nid.upper()], splitPhase=True)
        remain -= 1
        start = tmpStart
    m._partial = partial
    return out


# ---------------------------------------------------------------- comparison helpers
def parse_out(path, m):
    P = len(m.pipes)
    N = sum(1 for n in m.nodeList if n.kind == "load")
    O = sum(1 for n in m.nodeList if n.kind == "outlet")
    b = open(path, "rb").read()
    for off in range(0, len(b) - 12, 4):
        n = struct.unpack_from("<i", b, off)[0]
        if not (2 <= n <= 100000) or off + 4 + 20 * n > len(b):
            continue
        t = struct.unpack_from("<%dd" % n, b, off + 4)
        if t[0] != m.startTime / 3600.0:
            continue
        if b[off + 4 + 8 * n:off + 4 + 8 * n + 2].isdigit():
            break
    else:
        raise SystemExit("time block not found in .out")
    p = off + 4 + 8 * n + 12 * n
    full = struct.unpack_from("<%dd" % P, b, p); p += 8 * P
    steps = []
    for k in range(n):
        d = {}
        for name in ["load", "storm", "totstorm", "head", "status", "jump"]:
            d[name] = struct.unpack_from("<%dd" % N, b, p); p += 8 * N
        d["outflow"] = struct.unpack_from("<%dd" % O, b, p); p += 8 * O
        d["outhead"] = struct.unpack_from("<%dd" % O, b, p); p += 8 * O
        for name in ["flow", "vel", "depth", "froude", "dD", "overflow", "backwater",
                     "adjdepth", "adjvel", "qQ", "flowtype", "partial"]:
            d[name] = struct.unpack_from("<%dd" % P, b, p); p += 8 * P
        steps.append(d)
    return full, steps


def parse_rpt(path):
    tabs = []
    cur = None
    for ln in open(path, encoding="latin-1"):
        if ln.startswith("[Pipe Results]"):
            cur = {}
            tabs.append(cur)
            continue
        if ln.startswith("["):
            cur = None
            continue
        if cur is not None:
            tok = ln.split()
            if len(tok) == 7:
                try:
                    cur[tok[0]] = [float(x) for x in tok[1:]]
                except ValueError:
                    pass
    return tabs


def stat(label, pairs):
    ad = [abs(a - b) for a, b in pairs]
    rd = [abs(a - b) / abs(b) for a, b in pairs if abs(b) > 1e-6]
    big = max(range(len(ad)), key=lambda i: ad[i]) if ad else None
    print("  %-22s n=%5d  max|d| %.3e  max rel %.3e" % (label, len(pairs), max(ad) if ad else 0,
                                                       max(rd) if rd else 0))
    return max(ad) if ad else 0.0


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    inp = argv[1]
    args = dict(zip(argv[2::2], argv[3::2]))
    m = read_input(inp)
    tr = None
    if "--trace" in args:
        tf = open(args["--trace"], "w")
        tf.write("pipe,k,time_s,inflow_cfs,outflow_cfs\n")
        def tr(p, k, t, qi, qo):
            tf.write("%s,%d,%.17g,%.17g,%.17g\n" % (p.id, k, t, qi, qo))
    run(m, tr)
    res = outputs(m)
    print("model %s: %d pipes, %d nodes, %d report steps, H=%gs, nHyd=%d" %
          (inp, len(m.pipes), len(m.nodeList), m.reportStepCount, m.H, m.nHyd))
    tflow = m.mb["out"] + m.pipeStorage2 / 2.0
    cerr = abs(100.0 * (m.mb["in"] - tflow) / m.mb["in"]) if m.mb["in"] else 0.0
    print("mass balance: in %.3f  out %.3f  final pipe storage %.3f  check %.3f %%" %
          (m.mb["in"], m.mb["out"], m.pipeStorage2 / 2.0, cerr))
    if "--dump" in args:
        with open(args["--dump"], "w") as f:
            f.write("time_h,pipe,up_flow,dn_flow,flow_out,velocity,depth,dD\n")
            for r, d in enumerate(res):
                for i, p in enumerate(m.pipes):
                    f.write("%g,%s,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n" %
                            ((m.startTime + r * m.reportTimeStep) / 3600.0, p.id, d["up"][i], d["dn"][i],
                             d["flow"][i], d["vel"][i], d["depth"][i], d["dD"][i]))
    if "--out" in args:
        full, steps = parse_out(args["--out"], m)
        print("vs engine .out (all pipes x all report times):")
        stat("full flow", [(p.fullFlow, full[i]) for i, p in enumerate(m.pipes)])
        for key in ["flow", "vel", "depth", "dD", "qQ", "froude"]:
            stat("pipe " + key, [(res[r][key][i], steps[r][key][i])
                                 for r in range(len(steps)) for i in range(len(m.pipes))])
        for key in ["backwater", "adjdepth", "adjvel", "flowtype", "partial"]:
            stat("pipe " + key, [(res[r][key][i], steps[r][key][i])
                                 for r in range(len(steps)) for i in range(len(m.pipes))])
        stat("pipe overflow (=0)", [(0.0, steps[r]["overflow"][i])
                                    for r in range(len(steps)) for i in range(len(m.pipes))])
        for key in ["status", "jump"]:
            stat("node " + key, [(res[r][key][j], steps[r][key][j])
                                 for r in range(len(steps)) for j in range(len(res[r][key]))])
        loads = [n for n in m.nodeList if n.kind == "load"]
        stat("node head", [(res[r]["head"][n.id], steps[r]["head"][j])
                           for r in range(len(steps)) for j, n in enumerate(loads)])
        stat("node load", [(res[r]["load"][j], steps[r]["load"][j])
                           for r in range(len(steps)) for j in range(len(loads))])
        outs = [n for n in m.nodeList if n.kind == "outlet"]
        stat("outlet flow", [(res[r]["outflow"].get(o.id, 0.0), steps[r]["outflow"][j])
                             for r in range(len(steps)) for j, o in enumerate(outs)])
        stat("outlet head", [(res[r]["head"][o.id], steps[r]["outhead"][j])
                             for r in range(len(steps)) for j, o in enumerate(outs)])
    if "--rpt" in args:
        tabs = parse_rpt(args["--rpt"])
        pairs = []
        worst = (0, None)
        for r, tab in enumerate(tabs[:m.reportStepCount]):
            for i, p in enumerate(m.pipes):
                e = tab[p.id][2]
                v = round(res[r]["up"][i] * m.fuf + 0.0, 2)     # rpt prints flow * hm_flowUnitFactor
                pairs.append((v, e))
                if abs(v - e) > worst[0]:
                    worst = (abs(v - e), (r, p.id, res[r]["up"][i], e))
        print("vs engine .rpt [Pipe Results] D/S Flow column (2 decimals):")
        stat("flow (rounded)", pairs)
        print("  worst:", worst)
        # [Maximum Gravity Pipe Flows]: Fv_OutputPipeFlowRange + Fv_FlowRange
        mx = {}
        on = False
        for ln in open(args["--rpt"], encoding="latin-1"):
            if ln.startswith("[Maximum Gravity Pipe Flows]"):
                on = True
                continue
            if on and ln.startswith("["):
                break
            tok = ln.split()
            if on and len(tok) == 11:
                try:
                    mx[tok[0]] = [float(x) for x in tok[3:9]] + [tok[9]]
                except ValueError:
                    pass
        bad = 0
        for p in m.pipes:
            if p.id not in mx:           # [Report] Link : No -> no range table
                bad = -1
                break
            e = mx[p.id]
            qmax = p.maxFlow
            if abs(qmax) <= m.flowTol:
                dep = vel = 0.0
            else:
                dep = flow_depth(p, qmax / p.pc)
                vel = pipe_velocity(p, dep, 0.5 * p.diam, qmax / p.pc)
                if dep / p.diam < HMAP_DepthDiamTolerance:
                    dep = vel = 0.0
            tt = p.maxTime / 3600.0 + 0.0002777          # HU_GetMilTime
            hh = cint(tt)
            mm = cint((tt - hh) * 60.0)
            mine = ["%.2f" % (qmax * m.fuf), "%.2f" % vel, "%.2f" % dep, "%.2f" % (dep / p.diam), "%02d:%02d" % (hh, mm)]
            eng = ["%.2f" % e[2], "%.2f" % e[3], "%.2f" % e[4], "%.2f" % e[5], e[6]]
            if mine != eng:
                bad += 1
                if bad <= 5:
                    print("  max-table mismatch", p.id, mine, eng)
        if bad >= 0:
            print("vs engine .rpt [Maximum Gravity Pipe Flows] (flow, vel, depth, d/D, time): %d of %d pipes differ"
                  % (bad, len(m.pipes)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
