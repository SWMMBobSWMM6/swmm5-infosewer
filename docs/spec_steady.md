# InfoSewer engine (HeEngPkg) steady-state run: reimplementation spec

Scope: the STEADY-STATE run (`runOption 0`, `HM_RunType_Static`) of the InfoSewer engine (Innovyze/MWH Soft HeEngPkg, C++) for a gravity-only sanitary network: loading manholes, circular gravity pipes and outlets. Wet wells, pumps, force mains, siphons, junction chambers, design mode, water quality and storm (rational-method) loads are out of scope. Where they touch the code path, this spec says so.

Source read: the lower-cased Linux build copy `/home/claude/ise/port/build/src/`. Every `file:line` below refers to that copy.

Target: a C implementation inside EPA SWMM 5.2.4 that matches InfoSewer to about 1e-6 relative.

Status: `steady_ref.py` (same folder) implements this spec. It reproduces the engine's binary `.out` **bit for bit** (max difference 0.0, all 28 result arrays) on the 125_21A model and on 23 variants built to exercise every branch described here (section 6).

Written by Claude (AI) for SWMM Enablement, 2026-11-09, from a reading of the engine source. Every rule was checked against the engine by running it; none is a summary from memory.

---------------------------------------------------------------------------

## 0. Conventions that the whole calculation depends on

### 0.1 Internal units

| quantity | internal unit | conversion from input |
|---|---|---|
| flow | cfs | `x /= hm_flowUnitFactor` (loads, infiltration rates, Flow Tolerance, peak-flow curves) |
| elevation, length, manhole diameter, depth, head | ft | `x /= hm_lengthUnitFactor` |
| pipe diameter | ft | `x /= hm_pipeDiamUnitFactor` (12.0 for in, 304.8 for mm) |
| velocity | ft/s | |

Unit factors (he_define.h:40-51; selection he_input.cpp:10844-10885 and 10988-11006). Flow factors are *units per cfs*: CFS 1.0, GPM 448.831169, MGD **0.64632**, IMGD 0.5382, AFD 1.9837, LPS 28.316847, LPM 1699.010796, MLD 2.446575, CMH 101.940648, CMD 2446.575546. CFS/GPM/MGD/IMGD/AFD are US (length factor 1.0, diameter factor 12.0); LPS..CMD are SI (length 0.3048, diameter 304.8).

The binary `.out` is written in **internal units (cfs, ft) whatever the Flow Unit**; this was verified with MGD, GPM and LPS variants. The text `.rpt` multiplies back by `hm_flowUnitFactor` and the other factors (he_staticrun.cpp:1831-1906). InfoSewer's database (IEDB `PIPEOUT.DBF`, `MHOUT.DBF`) stores gpd, using 646316.64 gpd per cfs. That factor belongs to the GUI/database layer, not the engine.

### 0.2 Hard-coded constants (reproduce exactly, never "improve")

| constant | where | used for |
|---|---|---|
| `1.49` | he_hydrun.cpp:3824 | Manning coefficient (US), in every unit system |
| `0.6667` | he_hydrun.cpp:3824 | exponent of R (not 2/3) |
| `3.1416` | he_hydrun.cpp:2872, 2996, 3946 | pi inside circular geometry (`2.0*3.1416`) |
| `0.00001` | he_hydrun.cpp:3387, 3596, 3768 | slope used when slope == 0.0 |
| `4.66`, `5.33` | he_data.cpp:5508-5509 | full-pipe friction loss (Manning form) |
| `4.727`, `1.852`, `4.871` | he_data.cpp:5518-5523 | Hazen-Williams form (only if n >= 1.0) |
| `32.2` | he_staticrun.cpp:1316, 1318 | g in the reported Froude number |
| `9.81*1/0.3048` | he_hydrun.cpp:2654, 2677 | g in the circular critical-depth Froude test |
| `32.174` | he_runhead.cpp:2317 | g in the manhole minor loss |
| `3.1415926` (`pi`) | he_define.h:35 | infiltration "by surface" only |
| `HMAP_ZeroTolerance = 0.0001` | he_define.h:93 | several cut-offs (see text) |
| `HMAP_DepthTolerance = 1e-7` | he_define.h:95 | critical-depth stop, flow-type test, pipe-head tests |
| `HMAP_IterFlowTolerance = 1e-5` | he_define.h:96 | default `hm_iterTolerance` |
| `HMAP_WetAreaTolerance = 1e-5` | he_define.h:101 | velocity cut-off |
| `HMAP_DepthDiamTolerance = 0.01` | he_define.h:102 | d/D result filter |
| `HMAP_FlowTolerance = 1e-4` | he_define.h:108 | fixed flow cut-off (cfs) in depth/velocity/head code |
| `HMAP_MaxManningCoeff = 1.0` | he_define.h:97 | n >= 1.0 means Hazen-Williams C |
| `SYS_Zero = 1e-10`, `SYS_ZeroFloat = 1e-7` | stdafx.h:124-125 | headloss guard, flat-slope test |

Two flow tolerances coexist and are **not** interchangeable:

* `HMAP_FlowTolerance` = 1e-4 cfs, fixed.
* `hm_flowTolerance` = the [Options] "Flow Tolerance" divided by `hm_flowUnitFactor` (he_input.cpp:10468, 12321). Its default is 1e-4 in input units.

`hm_iterTolerance` = min(1e-5, hm_flowTolerance) (he_input.cpp:535, 12339-12340).

### 0.3 Why the exact constants and operation order matter

Normal depth is a bisection that stops on a *relative flow* error of 1e-4 (2.5). The depth it returns is wherever the error first falls below 1e-4, so it can be off the true root by about 1e-4 relative. A reimplementation that changes Q(d) by even 1e-6 (2/3 instead of 0.6667, M_PI instead of 3.1416, 1.486 instead of 1.49) can stop one step earlier or later, which moves depth, velocity, d/D and heads by up to about 1e-4.

To match to 1e-6, use the same constants, the same bisection bounds, the same update rule and the same stop tests. With that, the bisection visits the same midpoints, and the C result equals the engine's to the last bit except in astronomically rare tie cases. Python with glibc libm reproduces the engine exactly, so C with the same libm functions (`pow`, `acos`, `sin`, `sqrt`, `atan`) will too. Write expressions in the same order as quoted here; the C++ is evaluated left to right with no FMA contraction on x86-64.

### 0.4 Index conventions

* Node index: loading manholes in file order, then junction chambers, then outlets, then wet wells (he_input.cpp:3793, 4534, 4619, 4936). A loading manhole's node index equals its load-node index, and the code relies on that (`hm_pLoadNode[nodeIndex]`).
* Link index: pipes in [Pipes] file order come first (pipe index == link index, he_input.cpp:5201-5214).
* Parallel count: barrels = the [Pipes] last field + 1 (he_input.cpp:5207).

---------------------------------------------------------------------------

## 1. Order of computation

### 1.1 Pipeline

`Hx_EngineInputRun` (main.cpp:44) calls `H_ModuleManager`. That runs input (first pass, second pass, `Fv_DataOrganization` at he_input.cpp:197/222/239), then topology (`H_TopologyManager`, main.cpp:188), then `H_RunManager` (main.cpp:202, he_hydrun.cpp:85). `H_RunManager` runs:

1. `H_SetRunData` (he_hydrun.cpp:188): `HU_ResetVisit()` clears all link and node visit flags (202); `HH_PipeCapacity` sets `hm_pFullFlow[i] = HHU_GetFullFlow(i) * barrels` (1556-1565).
2. `H_StaticRunManager`, which builds `He_StaticRun`. `F_Initialize` (he_staticrun.cpp:192) zeroes the per-link arrays: peak, cover flow, cover count, unpeak, infiltration, link flow, velocity, backwater, adjusted velocity and depth (241-250). `F_Run` (276) then does:
   1. `Fv_StaticRun` (324): the node loop for flows (1.3-1.6), then `Fv_GetPipeHydraulic`, `Fv_GetCriticalDepth`, `Fv_GetPipeFlowDepth`, `Fv_GetOutletFlow` (401-404).
   2. `Fv_StaticFilterResult` (878): the d/D < 0.01 cut-off (2.8).
   3. `Fv_StaticGetHead` (905): `He_RunHead`, the HGL (section 3).
   4. `Fv_StaticOut` (987): the binary results (2.9).
   5. `Fv_StaticReport` (1769): the text report.

`HU_ResetVisit` is called before the node loop, so during the flow calculation a link's `visit` flag means "this link's flow has been set".

### 1.2 Topology (he_topology.cpp)

**Connections** (`Fv_SetLinkConnection`, 308-354). For each link `i` in index order, append `i` to `conn[from]` and to `conn[to]`. So `conn[node]` lists the node's links in link-index (file) order, upstream and downstream mixed.

**Flow-split nodes** (`Fv_SetFlowSplitLoadNode`, he_input.cpp:5680-5715). A node that is the from-node of more than one link gets kernel type `FlowSplit`; it must be a loading manhole, or input error 5322-5330. `F_IsKernelLoadNode()` is true for both `LoadManhole` and `FlowSplit` (he_data.cpp:993-1009).

**Node order** (`Fv_ManholeOrder`, 1377-1484). This is a Kahn-style topological sort that walks single-branch chains:

```
indeg[k] = number of links with to == k
clear node/link visit flags; order = []
repeat
  progress = false
  for i = 0 .. nodeCount-1:
    if visited(i): continue
    k = i
    while indeg[k] == 0:
      progress = true; visit k; order.append(k)
      nVisit = 0
      for l in conn[k]:                       # file order
        if visited(l): continue
        to = l.to                             # (nToNodeIndex)
        if visited(to): error "loop"
        visit l
        if indeg[to] > 0: nVisit++; indeg[to]--
      if nVisit == 1: k = to                  # follow a single downstream branch
      else: break
until len(order) == nodeCount or not progress
```

The flow results do not depend on which topological order is used: every node needs only its upstream links' values, and those are final. Any upstream-before-downstream order (SWMM's sorted order, for example) gives identical numbers, because the summation order is fixed by `conn[]` (1.4), not by the node order.

**Downstream and upstream links** (`Fv_LinkOrder`, 1525-1648, after `HU_ResetVisit`). For each node `k` in `order`:

* outlet: `down[k] = []` (1557-1561);
* `edge = len(conn[k]) == 1` and the link is not yet visited: `down[k] = [that link]`, `up[k] = []`, mark the link visited (1566-1579);
* `edge > 1`, non-split node: `down[k] = [first unvisited link in conn[k] with from == k]` (`Fv_GetDownstreamLink`, 1674-1689); `up[k]` = every other link of `conn[k]`, in conn order (1593-1612);
* `edge > 1`, split node: `down[k]` = every unvisited link of `conn[k]` in conn order; `up[k]` = the visited ones (1616-1640).

`hm_ppDownStreamLinks[k][0]` holds `len(down[k])` and `[1..]` the links; `hm_ppUpstreamLinkOrder[k][0]` and `hm_ppUpstreamLinkIndex[k][1..]` hold `up[k]`.

### 1.3 Loads at one loading manhole: `HHU_StaticGetLoading` (he_hydrun.cpp:984-1082)

The node's load list holds the manhole-line entries (up to 10, he_input.cpp:4239-4321), then any [Extra Loads] rows (he_input.cpp:8635-8659), in that order. Base loads are already in cfs (1.0).

```
unpeak = peak = cover = coverCount = 0
for each load (base, type, coverUnits):
   type 0 (Unpeakable): unpeak     += hm_loadMultFactor*base          (1024-1026)
   type 1 (Peakable)  : peak       += hm_loadMultFactor*base          (1030-1032)
   type 2 (Coverage)  : coverCount += 0.001*coverUnits                (1037; he_data.h:149)
                        cover      += hm_loadMultFactor*base          (1038)
infil = sum over l in down[k], if l is a pipe, of F_GetInfiltrateFlow(l)   (1062-1075)
unpeak += infil                                                       (1077)
```

* **Global Loading Multiplier** multiplies only base loads, never infiltration (1025, 1031, 1038).
* **Infiltration (I/I)** belongs to each pipe. It is charged to the pipe's *from* node and enters the accumulation as unpeakable flow, so it is never peaked. A split node gets the infiltration of all its downstream pipes.
* `F_GetInfiltrateFlow` (he_data.cpp:5420-5470) is rate x (length | diameter·length | pi·diameter·length | count), by type 1 | 2 | 3 | 4. It uses `pi` = 3.1415926 for type 3 and multiplies by the barrel count. Type 0 is none; type 5 (pattern) is EPS-only and uses the pattern value at time 0. Internal rate = input rate / flowF, x lengthF for type 1, x lengthF x diamF for types 2 and 3 (he_data.cpp:5119-5145). So the input rate is per ft (type 1), per ft·in (types 2 and 3), or per unit (type 4), in the input flow unit.
* Storm: when `hm_stormRunoffAnalysis == 1`, `HHU_RationalMethod` runs here (1054-1057); out of scope. Otherwise `hm_pNodeStormPeakFlow[] = 0` (memset he_input.cpp:1251).

### 1.4 Peaking functions

`HHU_StaticGetPeakFlow(B)` (he_hydrun.cpp:1423-1450). This returns the **peaked flow**, not a factor:

```
if no peaking curve: Qpeak = hm_peakFactK * pow(B, hm_peakFactP)      (1438)
else               : Qpeak = curve_Y(peakCurve, B)                    (1444)
```

`B` and the result are cfs, and `K` is "cfs based" (he_input.cpp:494). Changing Flow Unit does not rescale K. The K/P defaults 2.4 and 0.89 are the classic `Q = 2.4·Q^0.89` with Q in cfs.

`HHU_StaticGetCoverFlow(F, U)` (1458-1496), with U = accumulated coverCount in thousands of units:

```
if U < 0.0001: return 0
if no population curve:
   varC   = pow(U*hm_peakParamE, hm_peakParamC)
   factor = (hm_peakParamA/(hm_peakParamB + varC)) + hm_peakParamD
   return factor*F
else: return F * curve_Y(popCurve, U*1000)        # x axis = raw unit count
```

With the defaults a=5, b=0, c=0.2, d=0, e=1 this is the Babbitt formula `PF = 5/P^0.2`, with P in thousands.

`curve_Y` = `HD_CurveData::F_GetInterpolatedY` (he_data.cpp:10224-10306):

```
if size<=0 or x<0: 0
if x <= X0: (X0 == 0 ? 0 : x*Y0/X0)                    # straight line through the origin
elif x >= Xlast: size==1 ? Y0 : (|dx|<1e-4 ? Ylast : Ylast + (x-Xlast)*dy/dx)   # last segment extended
else first i with X[i] >= x: (|dx|<1e-4 ? Y[i] : Y[i] - (X[i]-x)*dy/dx)
```

Peak-flow curves (type 4) have both X and Y scaled by `1/hm_flowUnitFactor` (multiply by the reciprocal, he_input.cpp:12269-12273, he_data.cpp:9838-9900). Population curves (type 7) are not scaled.

### 1.5 The node loop: `Fv_StaticRun` (he_staticrun.cpp:344-396)

Per link the engine keeps `Q` (`hm_pLinkFlow`), `UP` (`hm_pUnpeakFlow`, base unpeakable **including** infiltration), `PK` (accumulated peakable *base*), `CV` (accumulated coverage *base*), `CC` (accumulated cover units/1000), `INF` (accumulated infiltration) and `STORM`.

```
for k in order:
  if len(conn[k]) == 1 and k is a loading manhole:            (355)
      Fv_GetSingleSourceNodeFlow(k)                           (423-510)
  else:
      Fv_GetConnectNodeFlow(k)                                (518-673)
```

**Source node** (`Fv_GetSingleSourceNodeFlow`), with `l = conn[k][0]`:

```
static[k] = HHU_StaticGetTotalLoading(k)                      (449)
(UP,PK,CV,CC,INF)[l] = HHU_StaticGetLoading(k)                (450-455)
Q[l] = UP[l] + PeakFlow(PK[l]) + CoverFlow(CV[l],CC[l]) + storm[k]   (479-484)
STORM[l] = storm[k]; visit(l)
```

There is no zero-flow test here: a source node with no load gives `Q = 0 + K·0^P + 0`, which is 0 for P > 0 and K if P == 0.

**Other nodes** (`Fv_GetConnectNodeFlow`) are loading manholes with several links, outlets and chambers:

```
up = pk = cv = cc = inf = 0
if k is a loading manhole:
   static[k] = HHU_StaticGetTotalLoading(k)                   (547)
   (up,pk,cv,cc,inf) = HHU_StaticGetLoading(k)                (548-553)
for l in conn[k]:                                             HHU_StaticAddPrecedeFlow (he_hydrun.cpp:371-383)
   if visited(l): up += UP[l]; pk += PK[l]; cv += CV[l]; cc += CC[l]; inf += INF[l]
inflow = up + PeakFlow(pk) + CoverFlow(cv,cc) + storm[k]      (596-599)
nd = len(down[k])
if |inflow| < hm_flowTolerance (and not wet well/chamber):    (606-621)
   for i = 1 .. nd-1:  zero Q,UP,PK,CV,INF,STORM of down[k][i] and visit it   # i starts at 1: see 5.1
   return
if nd == 1:  HHU_StaticSingleDownLink (he_hydrun.cpp:422-433):
   Q,UP,PK,CV,CC,INF of down[k][0] = inflow,up,pk,cv,cc,inf;  STORM = storm[k];  visit
elif nd > 1: HHU_StaticFlowSplitLinks(k, inflow)  (1.6), then for each l in down[k]:   (646-666)
   ratio = Q[l]/inflow
   UP[l]=ratio*up; PK[l]=ratio*pk; CV[l]=ratio*cv; CC[l]=ratio*cc; INF[l]=ratio*inf; STORM[l]=ratio*storm[k]; visit
```

`HHU_StaticAddPrecedeFlow` adds every *visited* connecting link. In a topological order those are exactly the upstream links whose flow has been set; the downstream links are not yet visited. It adds the five components but not `Q` and not `STORM`.

`HHU_StaticGetTotalLoading` (he_hydrun.cpp:933-976) is `unpeak + PeakFlow(peak) + CoverFlow(cover, count) + storm[k]` over the node's *own* loads only, plus its own downstream infiltration. This is the per-manhole "static flow" result (`STAT_FLOW`).

**Where peaking is applied.** At *every* node, to the *accumulated* peakable base flow (and the accumulated coverage base and units) of all loads upstream. A pipe's flow is `ΣUnpeak + ΣInfil + K·(ΣPeakableBase)^P + PF(ΣUnits)·ΣCoverBase`, summed over everything upstream. Peaked flows are never summed, so flow is not conserved at junctions: downstream `K(B1+B2)^P` is less than `K·B1^P + K·B2^P` for P < 1. **The reimplementation must accumulate base components along the tree and evaluate the peak at each pipe. It must not route peaked flows through SWMM's node continuity.**

**Outlet flow** (`Fv_GetOutletFlow`, 842-870): `outFlow[to] += Q[i]` over pipes whose to-node is an outlet. This runs *before* the negative-flow clip of 2.9.

### 1.6 Flow splits (more than one downstream pipe)

`HHU_StaticFlowSplitLinks` (he_hydrun.cpp:477-581). It returns an error code that the caller ignores (he_staticrun.cpp:646).

1. `Fc_GetDownstreamSplitData` (589-682): for each `l` in `down[k]` not yet visited, set `Q[l] = 0` and classify it as *hydraulic* (`EN_SplitFlowType_Hydraulic`, the default for every pipe, he_input.cpp:5224) or *percent/curve*.
2. Percent pipes (509-541): `HHU_CompPercentSplitFlow` (1744-1878) collects percentages. If their sum is 0 the ratio is 0. If the sum is below 100 and there is at least one hydraulic pipe, they are used unscaled; with no hydraulic pipe they are scaled by 100/sum. If the sum is above 100 they are scaled by 100/sum. Then `Q[l] = 0.01*pct*inflow` (536) and `splitSum += Q[l]`.
   Note: `F_SetSplitFlowPercent` sets the type to Percent only when pct > 0 (he_data.cpp:4460-4490), so **a 0 % split stays hydraulic**.
3. No hydraulic pipe: done. Exactly one: `Q[h] = inflow - splitSum` (553-558).
4. Several hydraulic pipes: `Fc_StaticMultiHydSplitFlow(hyd, Qh = inflow - splitSum)` (690-855). This is a common water level at the node:

```
if Qh < hm_flowTolerance: Q[h]=0 for all; visit; return            (704-714)
minInv   = min over h of fromInv[h];  maxCrown = max over h of fromInv[h] + D[h]    (721-735)
wl = 0.5*(minInv + maxCrown)
out = sum over h of Qpipe(h, wl)                                    (743-769)
   Qpipe(h, wl): d = wl - fromInv[h]; d <= 0 ? 0 : barrels[h]*HHU_GetLinkFlow(h, min(d, D[h]))
lo = minInv; hi = maxCrown
if out > Qh: hi = wl;  if out < Qh: lo = wl;  wl = 0.5*(hi + lo)    (774-780)
while |out - Qh|/out > hm_iterTolerance:                            (785)
   out = 0; allFull = true
   for h: d = wl - fromInv[h]
          if d <= 0: Q[h]=0; allFull=false
          else: Q[h] = barrels[h]*HHU_GetLinkFlow(h, min(d,D[h])); if d < D[h]-0.0001: allFull=false
          out += Q[h]
   if allFull: tf = sum fullFlow[h]; if tf <= 0 break
               Q[h] = Qh*fullFlow[h]/tf for all h; break           (819-835)
   if out > Qh: hi = wl;  if out < Qh: lo = wl;  wl = 0.5*(hi + lo)
visit all h
```

The stop test uses `out` from the *previous* water level, and the `Q[h]` left behind belong to that level. So the legs sum to `out`, not exactly to `Qh` (relative error up to `hm_iterTolerance`), and the base components prorated by `ratio = Q[l]/inflow` (he_staticrun.cpp:652) inherit that small imbalance. `hm_iterTolerance` is a *flow* tolerance used here as a *relative* tolerance (0.2).

In the 125_21A model CTS-460 splits by this rule into CTS-460:CTS-450 = 1.4895370141997326 cfs and CTS-460:CTS-453 = 11.433979985800265 cfs, reproduced exactly.

Curve-based (type 2) and in/out curve (type 3) splits use `HU_GetFlowSplitPercentFromInFlow` and were not traced; they are out of scope for `steady_ref.py`.

---------------------------------------------------------------------------

## 2. Pipe hydraulics (circular)

All of this is per barrel: `q = Q/barrels`, `D = diameter (ft)`, `r = 0.5*D`.

### 2.1 Geometry (he_hydrun.cpp)

```
HHU_GetWetArea      (2995-2999): angle  = (d - r)/r
                                 angle2 = (2.0*3.1416) - (2.0*acos(angle))
                                 A      = 0.125*(angle2 - sin(angle2))*((2.0*r)*(2.0*r))
HHU_GetWetPerimeter (3945-3947): P = (0.5*angle2)*(2.0*r)          # same angle2
HHU_GetWetWidth     (2873)     : T = 2.0*sqrt(d*(2*r - d))
```

Because 3.1416 is not pi, `angle2(d=0) = 1.469e-5` (A(0) is about 5e-16·D², not 0) and `angle2(D) = 6.2832`.

### 2.2 Manning flow

`HHU_GetFlow` (3753-3831), circular branch:

```
if slope == 0.0: slope = 0.00001
Q = (1.49/n)*A*pow(A/P, 0.6667)*sqrt(slope)            # evaluated left to right
```

The coefficient is **1.49**, in all unit systems: the internal units are always ft and cfs.

`HHU_GetFullFlow` (3558-3628), per barrel: `s = slope (0 -> 0.00001); s < 0 ? -1.0 : HHU_GetFlow(D, r, s, n)`. **Adverse pipes have full flow -1 per barrel.** `hm_pFullFlow[i] = HHU_GetFullFlow(i)*barrels` (he_hydrun.cpp:1563-1564).

`HHU_GetLinkFlow(i, d, r=-1, n=-1)` (3312-3402), circular:

```
if hm_pAdversePipe[i]: return hm_pFullFlow[i]          # the TOTAL (-barrels), whatever d is  (3378)
slope = slope (0 -> 0.00001); r<0 -> 0.5*D; n<=0 -> pipe n
return HHU_GetFlow(i, d, r, slope, n)
```

`slope = (fromInvert - toInvert)/length`, computed from the raw input values (`Fv_UpdateSlopeData`, he_data.cpp:5600-5622, called from the setters 4033-4108). Slope type: Flat if `|slope| < 1e-7`, Standard if `slope > 0`, else Adverse. Separately, `hm_pAdversePipe[i] = fromInv < toInv` (he_input.cpp:11797-11806). A pipe with a tiny negative slope is therefore Flat *and* adverse.

### 2.3 Surcharge flag: `Fv_GetPipeHydraulic` (he_staticrun.cpp:681-724)

```
full[i] = (Q[i] > hm_flowTolerance) and (Q[i] >= hm_pFullFlow[i])
```

Flow is **not capped** at full flow. A surcharged pipe keeps its full routed flow, gets depth = D and velocity = q/A(D). `HHU_StaticCompSurcharge` (1504), which would cap it, is "no longer in use". Every adverse pipe with flow is full (full flow -1). Flat pipes use S = 1e-5 for capacity and are usually full.

### 2.4 Critical depth: `Fv_GetCriticalDepth` (732-768) and `HHU_GetCriticalDepth` (2586-2685)

```
if |Q| < hm_flowTolerance: yc = 0
else: q = min(Q, fullFlow)/barrels;  yc = Crit(q, D)

Crit(q, D):                                          # circular, he_hydrun.cpp:2647-2678
  if |q| < hm_flowTolerance: return 0
  d0 = 0; d1 = D; y = 0.5*(d1 - d0); yref = 0
  A = area(y); if A == 0.0: return 0
  Fr = (q/A)/sqrt(A/T(y)*9.81*1/0.3048)
  while |y - yref| > 1e-7:
     yref = y
     if Fr < 1: d1 = y; y = 0.5*(d1 - d0)            # BUG kept: should be (d1+d0)/2, see 5.4
     else:      d0 = y; y = d0 + 0.5*(d1 - d0)
     A = area(y); if A == 0.0: return 0
     Fr = (q/A)/sqrt(A/T(y)*9.81*1/0.3048)
  return y
```

Despite the bug, results are within about 4e-7 relative of the true critical depth on 125_21A, so the bug costs some iterations but no accuracy. Reproduce it anyway, since the stopping point depends on the sequence. For adverse pipes `min(Q, -barrels)/barrels = -1`, so Fr < 0 always and y halves toward zero: the reported critical depth is 6.7055e-08 ft for a 2.25 ft adverse pipe.

### 2.5 Normal depth: `Fv_GetPipeFlowDepth` (776-834) and `HHU_GetFlowDepth` (2164-2311)

```
q = Q/barrels
depth = full ? D : FlowDepth(q)

FlowDepth(q):                                         # circular branch 2231-2305
  if |q| < 0.0001 (HMAP_FlowTolerance): return 0      (2176, 2245)
  dL = 0; qL = 0; dU = D; qU = HHU_GetLinkFlow(i, dU, r, n)
  if qU <= q: return D                                (2264)  # adverse pipes always exit here
  y  = 0.5*(dU + dL); qN = HHU_GetLinkFlow(i, y, r, n)
  e1 = q - qL; e2 = q - qU; ea = q - qN
  do:
     if ea*e2 <= 0:   dL = y; qL = LinkFlow(dL); e1 = q - qL
     elif ea*e1 <= 0: dU = y; qU = LinkFlow(dU); e2 = q - qU
     y  = 0.5*(dL + dU)
     ea = q - LinkFlow(y)
  while |ea|/q >= 0.0001                              (2293)
  if slope type != Standard: y = min(y + 0.5*length*(-slope), D)    (2299-2305)
  return y
```

* The tolerance is relative *flow* error 1e-4, there is no iteration cap and the start value is D/2.
* The upper bracket is Q(D), the full-pipe flow, not the maximum at about 0.94D. For q just below Q(D) the root found depends on the midpoint sequence, so mirror it exactly.
* For q < 0 (only possible with negative curve extrapolation, 5.9) neither branch fires and the stop test `|ea|/q` is negative, so y = D/2 is returned at once.
* The flat/adverse correction `y + 0.5·L·(-slope)` is effectively dead for circular pipes: adverse pipes return D before it, and flat pipes have |slope| < 1e-7.

### 2.6 Velocity: `HHU_GetPipeVelocity` (3079-3155), `checkCutoff = false` by default (he_function.h:155)

```
if depth > 0: V = Vel(depth, r, min(q, Q/barrels))     (he_staticrun.cpp:813-821)
else V = 0
Vel(d, r, q): A = area(d, r); if A < 1e-5: 0; if q < 0.0001: 0; else q/A
```

A per-barrel flow below 1e-4 cfs (about 65 gpd) reports zero velocity even when depth > 0.

### 2.7 Head loss of a full pipe: `HD_PipeData::F_GetHeadloss` (he_data.cpp:5480-5530)

```
if n < 1.0: hf = 4.66*n*n*L*Q*Q / pow(D, 5.33)      (if |pow(D,5.33)| > 1e-10, else 0)
else      : hf = 4.727*L*pow(Q,1.852) / (pow(C,1.852)*pow(D,4.871))
```

`Q` is the **total** link flow (all barrels) and `D` is one barrel's diameter (5.6). `4.66` is the Manning full-pipe constant for k = 1.486 (exactly 4.66 with 1.486; 4.637 with 1.49), so it is inconsistent with the 1.49 used for capacity (5.5). n > 1.0 is rejected at input (he_input.cpp:11696-11735); n = 1.0 exactly would select Hazen-Williams.

### 2.8 Result filter: `Fv_StaticFilterResult` (878-897) and `HHU_FilterDepthVelocity` (3262-3304)

```
if not (depth/D >= 0.01): depth = 0; V = 0
```

This runs **after** the depth and velocity calculation and **before** the head calculation, so the heads use the filtered depth. Critical depth is not filtered.

### 2.9 Reported pipe results: `Fv_PipeData` (he_staticrun.cpp:1232-1506)

The output order is the order of the `.out` arrays.

| # | array | rule |
|---|---|---|
| 1 | flow | `Q[i] < 0 -> 0` first (1244-1248); then Q |
| 2 | unpeak | `UP - INF` (1254-1255): unpeakable base without infiltration |
| 3 | peak | PK, accumulated peakable **base**, not peaked |
| 4 | cover | CV, accumulated coverage base, not peaked |
| 5 | infil | INF |
| 6 | storm | STORM |
| 7 | flowtype | `+1` if `ftype == 1` (pressurised, 3.6) else `-1` (1279-1285) |
| 8 | full | `hm_pFullFlow` (all barrels; adverse = -barrels) |
| 9 | rescap | `max(full - Q, 0)` (1295-1299) |
| 10 | vel | V (filtered) |
| 11 | dD | `depth/D` (1380-1381) |
| 12 | qQ | `Q / HHU_GetLinkFlow(i, D)` (1388-1392). The denominator is ONE barrel's full flow, so q/Q is x barrels too large for parallel pipes. For adverse pipes `HHU_GetLinkFlow` returns `hm_pFullFlow = -barrels`, so q/Q = `-Q/barrels` (negative) |
| 13 | depth | filtered depth |
| 14 | crit | yc |
| 15 | froude | below |
| 16 | covcnt | `1000*CC` (1410-1413), units |
| 17 | backwater | below |
| 18 | adjdepth | below |
| 19 | adjvel | below |
| 20 | partial | `+1` if `partial[i]` (3.4) else `-1` (1492-1498) |

**Froude** (1309-1321), with `T = T(depth)` and `A = area(depth)`:

```
if T <= 0 and depth > 0: Fr = V/sqrt(32.2*depth)          # full pipes: T(D) == 0 exactly
elif A > 0:              Fr = V/sqrt(32.2*(A/T))          # depth 0: A/0 = inf -> Fr = 0
else: Fr = 0
```

**Backwater adjustment** (1429-1460), with `hu = head[from]`, `hd = head[to]` (3, output heads):

```
if hd > depth + toInv:
   adj  = min(D, 0.5*((hd - toInv) + (hu - fromInv)))
   Vadj = Vel(adj, D/2, Q_total)            # TOTAL flow, not per barrel (5.6)
   backwater = 1
else: adj = depth; Vadj = V; backwater = 0
```

The flow in `Vadj` is the clipped `Q`. A negative `adj` is not guarded (acos of a value below -1 gives NaN); it is possible only in contrived zero-flow cases.

**Loading-manhole arrays** (`Fv_LoadNodeData`, 1042-1108), in load-node order:

| array | rule |
|---|---|
| base | `F_GetTotalBaseLoad`: Σ base over all loads and extra loads, cfs, **no** multiplier, no peaking, no infiltration (he_data.cpp:2271-2290) |
| static | `HHU_StaticGetTotalLoading` (1.5) |
| storm | `hm_pNodeStormPeakFlow` (0) |
| head | output head (3) |
| surch | `head >= rim ? 1 : 0` (1085-1093) |
| jump | 1 if any *upstream* pipe has Fr > 1 **and** any *downstream* pipe has Fr < 1, else 0 (1326-1369) |

**Outlets** (`Fv_OutletData`, 1152-1187): outlet flow (1.5, unclipped) and output head.

**`.out` layout.** Header, then 20 pipe arrays x pipeCount doubles, 6 x loadNodeCount, [chamber head], 2 x outletCount, [wet-well, force-main, pump blocks], then a module marker. All arrays are little-endian float64 (see `/home/claude/ise/port/compare_m125.py`). Bytes 9-16 of the header hold a value that changes from run to run, so compare the data only.

IEDB fields (from `compare_m125.py`): PIPEOUT `PIPE_FLOW, VELOCITY, D_OVER_D, FLOW_DEPTH, Q_OVER_Q, CRIT_DEPTH, FROUDE_NUM, FULL_FLOW, ADJ_DEPTH, ADJ_VELOC`; MHOUT `BASE_FLOW, STAT_FLOW, HEAD`; OUTLEOUT `FLOW, HEAD`. Flows are x 646316.64 gpd per cfs.

### 2.10 Other shapes (not traced in detail)

`[Channel Shape]` assigns `HM_ChannelType` (he_enum.h:75-84): Circular 0, Rectangular closed/open 1/2, Trapezoidal closed/open 3/4, Triangular closed/open 5/6, Parabolic closed/open 7/8. Non-circular pipes take separate routines: `HHU_GetTrap*` and `HHU_GetParab*`, the `F_GetNonCircularHeadLoss` closed-channel head loss (he_data.cpp:5541), and open-channel branches in the head code. They differ from circular in ways a port would have to copy:

* the depth bisection stops at relative error 0.001, not 1e-4 (2396, 2489);
* the flat/adverse correction is `+0.5·L·slope` with **no** sign flip and no cap (2398-2403, 2492-2497);
* critical depth uses g = 32.2 (2717, 2741);
* the trapezoid wet perimeter uses `sqrt(1 + Z)` instead of `sqrt(1 + Z²)` (5.11);
* the parabolic area uses the full top width (5.11).

---------------------------------------------------------------------------

## 3. Head (HGL) calculation: `He_RunHead` (he_runhead.cpp)

Inputs (he_staticrun.cpp:934-952): the filtered depths, `Q` (before the negative clip) and wet-well/pump data (none here). Two options are **fixed** because no input key sets them (he_input.cpp:489-490): `hm_optHeadZeroLinkFlow = ByInvert` and `hm_optHeadUpstreamNode = Adjust`.

### 3.1 Structure

* `Fv_SetUpstreamNode` (955-1029): for each node, `U[k]` = the links of `conn[k]` that are not in `down[k]`, in conn order, each paired with its from-node.
* Every node starts with `hyd[k] = -1e37` (899-904) and `visitLevel = 0`. "Fully visited" means `visitLevel == len(U[k])` (he_data.cpp:14080), so nodes with no upstream links are always fully visited.
* `F_Run` (554-629):

```
for each outlet o in outlet order: OutletHead(o); DFS(o)
remain = number of flow-split nodes; start = 0
while remain > 0:
   tmpStart = end = splitCounter
   if start >= end: return
   for i = start .. end-1: SplitDFS(split node whose order == i)
   remain--; start = tmpStart

DFS(k):                                   # Fv_DepthSearch 1162-1231
   if fully visited(k): return
   if k is a flow-split node: if not yet ordered: order[k] = splitCounter++;  return     # suspend
   for (u, l) in U[k]: NodeHead(u, k, l); visitLevel[k]++; DFS(u)
SplitDFS(k):                              # Fv_FlowSplitDepthSearch 1239-1288
   same loop as DFS without the split-suspend test
```

The calculation therefore runs **upstream from each outlet**. A split node stops the descent until all outlet trees are done; then the split nodes are resumed level by level, in the order they were first reached. That order matters (5.10).

### 3.2 Outlet head: `Fv_GetOutletHead` (1296-1345)

```
freeFall = max over (u,l) in U[o] of depth[l] + toInv[l]
hyd[o] = out[o] = max(submergeHead[o], freeFall)
```

`submergeHead` defaults to -1e37 (he_data.cpp:2809), and [Outlet Heads] sets it. In 125_21A the outlet pipe "13" is flat and full, so the outlet head = 43.5 + 2.75 = 46.25.

### 3.3 Pipe head: `Fv_GetPipeHead(l, u, hd = hyd[d])` (1404-1583), circular pipes

Standard slope (1422-1485):

```
if Q[l] < fullFlow[l]:                                  # partly full, incl. zero flow
   if Q[l] <= 0.0001: hu = max(hd, bottomRef[u])        # ByInvert (1454-1460)
   else:              hu = depth[l] + fromInv[l]        # normal-depth water surface
   hu = AdjustPipeHead(l, hd, hu)                       # 3.4
else:                                                   # surcharged (1469-1478)
   hu = max(hd, D + toInv) + headloss(Q[l])
```

`bottomRef[u]` = the minimum over all pipes touching node u of (fromInv for pipes leaving u, toInv for pipes entering u) (he_input.cpp:5256-5270, `F_SetBottomRefLevel` keeps the minimum, he_data.cpp:1417-1430).

Flat or adverse slope (1490-1560):

```
if Q[l] > 0.0001:
   dH  = depth[l] >= D ? headloss(Q[l]) : 0
   dd  = max(hd, depth[l] + toInv)
   hu  = max(fromInv + depth[l], dd + dH);  hd = dd
   hu  = AdjustPipeHead(l, hd, hu)                      # returns at once: Q >= fullFlow - 1e-4
else:
   hu = fromInv;  hd = max(hd, toInv)
   if hd > toInv:     partial[l] = true
   if hd > toInv + D: adjFull[l] = true
```

Then, always (1568-1577):

```
AdjustNodeHead(u, hu)                                   # 3.5; sets hyd[u]
ftype[l] = (hyd[u] >= fromInv + D + 1e-7) ? 1 : 0       # pressurised flag, using u's head *now*
```

### 3.4 `Fv_AdjustPipeHead(l, hd, &hu)` (1591-1714)

```
if Q >= fullFlow - 0.0001: return                         # (1605)
if open channel: return; if to-node is a junction chamber: return
downRef = toInv + D
if Q > 0.0001:
   if hd > hu and hd < downRef:                           # tailwater above normal depth, below crown
      up = hd + headloss(Q)                               # full-pipe loss over the whole length (5.7)
      if up > hu: hu = up; partial = true; if hu >= fromInv + D: adjFull = true
   if hd <= downRef + 1e-7: return
   deltaLevel = hd - (toInv + depth)
   refLength  = deltaLevel/sin(atan(slope))
   dH = (refLength > length) ? headloss(Q) : 0
   crossing = hd + dH
   if hd > toInv: partial = true
   if crossing > hu: hu = crossing; partial = true; if hu >= fromInv + D: adjFull = true
```

`partial[l]` is the reported "partial full" flag (array 20). `adjFull` is not reported in the static run.

### 3.5 Node head: `Fv_AdjustNodeHead(u, h)` (1910-1991), loading manhole and flow-split

```
(Fv_AdjustLoadNodeSurcharge: Dynamic mode only -> no-op in steady state, 2011-2012)
# Fv_AdjustLoadNodeDownstream (2152-2250):
dh = -1e37
for l in down[u]:
   if Q[l] > 0.0001: dh = depth[l] + fromInv[l]
   h = max(h, dh)                          # circular / closed pipes
# Fv_AdjustNodeMinorLoss (2258-2340):
for i, l in enumerate(down[u]):
   if Q[l] > 0.0001:
      V = Vel(D_l, D_l/2, Q[l])            # full-bore area, TOTAL flow; shape looked up with pipe index i (5.8)
      h += |K_u| * (V*V)/(2.0*32.174)      # K_u = loading-manhole field 4; added once per flowing downstream pipe
if u is a flow-split node: hyd[u] = max(hyd[u], h)  else: hyd[u] = h
out[u] = (head type Unlocked) ? min(hyd[u], rim[u]) : hyd[u]
hyd[u] = out[u]                            # "2002-11-26": the capped head is carried upstream
```

* **Surcharge at manholes.** No surcharge-volume or flooding calculation is done in steady state. An Unlocked manhole (the default: Manhole Head Type 0) is simply capped at its rim, and the cap propagates upstream. A Locked manhole (1) keeps the HGL above the rim; the `surch` flag reports `head >= rim`.
* **Minor loss** is applied whether or not the downstream pipe is surcharged; the conditional is commented out (2310-2311). It uses the full-bore velocity even for part-full pipes. With K = 0, as in 125_21A, it vanishes.
* `Fv_AdjustLoadNodeUpstream` exists but is commented out of every call (1936, 1951).

### 3.6 Flow type (pressurised) flag

`ftype[l] = 1` when the upstream node's carried head, at the moment pipe `l` is processed, is at least `fromInv + D + 1e-7` (1574-1577). For a split node's legs this uses the head *so far* (5.10).

### 3.7 Hydraulic-jump flag

See 2.9: it is computed from the reported Froude numbers. There is no head jump in the HGL.

---------------------------------------------------------------------------

## 4. Inputs that affect the steady-state result

The format is tab- or space-separated tokens; `;` starts a comment (he_input.cpp:13185-13186); a token may be quoted. Numbers go through `strtod` after replacing the first `,` with `.` (`Fv_GetDouble`), except where `atof` is noted.

### 4.1 [Options] (`Fv_ReadOption`, he_input.cpp:9583-10808)

Each line is `words : value`. The value is the last token; the `:` token is skipped when the words are matched.

| key (words) | variable | default (he_input.cpp:457-547) | effect |
|---|---|---|---|
| `Flow Unit` | `hm_flowUnitFlag` -> `hm_flowUnitFactor`, length/diameter factors | CFS | 0.1 (10452-10460) |
| `Flow Tolerance` | `hm_flowTolerance` (then /flowF), `hm_iterTolerance = min(1e-5, it)` | 0.0001 | zero-inflow cut, full flag, crit depth, split solver (10464-10469) |
| `Global Loading Multiplier` (first token starts with "Global") | `hm_loadMultFactor` | 1.0 | x base loads (10327-10331) |
| `Peaking Factor k` | `hm_peakFactK` | 2.4 | 1.4 (10347-10352); the first letter of word 3 is tested |
| `Peaking Factor p` | `hm_peakFactP` | 0.89 | 1.4 |
| `Peaking Factor Curve ID` (word 3 = "Curve", more than 4 tokens) | `hm_peakCurveID` -> `hm_peakCurveIndex` | none (-1) | replaces K·B^P (10365-10374, 11446-11449) |
| `Peaking Parameter a`..`e` | `hm_peakParamA..E` | 5, 0, 0.2, 0, 1 | coverage formula 1.4 |
| `Coverage Based ... : ID` (words 1-2 = "Coverage Based") | `hm_peakPopCurveID` -> `hm_peakPopCurveIndex` | none | replaces the a-e formula (10484-10490, 11453-11456) |
| `Manhole Head Type` | `hm_optLoadNodeHeadType` | 0 = Unlocked | default per-manhole head type (9668-9681, 11131-11139) |
| `Storm Runoff Analysis` | `hm_stormRunoffAnalysis` | 0 | 1 = rational-method storm loads (out of scope) (9720-9737) |
| `Advanced Forcemain Network Support` | `hm_pressureNet` (`YES`=1) | 0 | must be 0 for this spec (9656-9662) |

Read by the engine but without effect on a gravity steady run: Loading Pattern Type, Default Analysis/Design Criteria Curve ID, the other design keys, Minimum Pipe Segment, Minimum Travel Time, quality keys and Pump Flow Continuity. The criteria-curve checks in `Fv_CheckDesignData` return at once outside Design mode (he_input.cpp:11888). A run with `Default Analysis Criteria Curve ID : -` gave identical results. There is **no** option for `hm_optHeadZeroLinkFlow` or `hm_optHeadUpstreamNode`.

### 4.2 Sections

| section | fields (0-based tokens) | variable / rule |
|---|---|---|
| `[Loading Manholes]` (`Fv_ReadLoadNode`, 3733-4371; at least 58 tokens, 3839) | 0 ID, 1 rim, 2 manhole diameter, **3 minor-loss K**, 4 area, 5 C, 6 slope, 7 length, 8 Tc, 9 hyetograph, 10 hydrograph, 11-17 runoff/infiltration parameters, then 10 load entries at 18+4i: **base, type (0/1/2), pattern, cover units** | base via `atof`; entries with \|base\| < 0.0001 *in input units* are dropped (4252); cover units used only for type 2 (4302-4304) |
| `[Outlets]` (4516-4631) | 0 ID, 1 rim, 2 diameter | |
| `[Pipes]` (4985-5342; exactly 9 tokens) | 0 ID, 1 from, 2 to, 3 fromInv, 4 toInv, 5 length, 6 diameter (in or mm), 7 n, 8 extra barrels | barrels = token 8 + 1; slope from the raw values; bottomRef of end manholes |
| `[Flow Splits]` (7875-8024; 4 tokens) | pipe ID, type (0 hydraulic, 1 percent, 2 curve, 3 in/out), percent, curve ID | percent > 0 makes the pipe a Percent split; rows with `-` in fields 2 and 3 are skipped |
| `[Infiltrate Rates]` (8031-8164; 5 tokens) | pipe ID, type (0-5), rate, count (type 4), pattern (type 5) | 1.3 |
| `[Extra Loads]` (8488-8668; 5 tokens) | manhole ID, base, type, pattern, cover units | appended to the manhole's loads; **no** 1e-4 filter; base via `atof` |
| `[Load Node Head Type]` (8676-8734) | manhole ID, 0 Unlocked / 1 Locked | overrides the option |
| `[Outlet Heads]` (8742-8796) | outlet ID, submerged head | `submergeHead` (3.2) |
| `[Curves]` (1648-1781) | `ID type x y` then `x y` continuation lines | type 4 = peak flow (X, Y in flow units), type 7 = population (X = units, Y = factor) |
| `[Patterns]`, `[Time]`, `[Controls]`, `[Status]` | | not used by the steady gravity calculation |

---------------------------------------------------------------------------

## 5. Odd or buggy behaviour that a faithful port must reproduce

Each item quotes the code as it stands in the build copy.

**5.1 Zero-inflow branch skips the first downstream link** (he_staticrun.cpp:608-611):
```
   for (i = 1;  i < downSplitCount;  i++)
    {
     int  linkIndex = hm_ppDownStreamLinks[nodeIndex][i+1];
```
`i` starts at 1, so `down[k][0]` is never zeroed or marked visited. This is harmless in practice: the link's arrays are still at their initial 0 and the unvisited link adds nothing downstream. The "zero" variant (all loads from CTS-460 up removed) matched.

**5.2 Split error ignored** (he_staticrun.cpp:646): `HHU_StaticFlowSplitLinks(nodeIndex,  inflow);`. The returned `HMAP_RunErrorSplitFlowRange` is dropped.

**5.3 Multi-hydraulic split continuity** (he_hydrun.cpp:785): `while (fabs(outFlow-hydCalcFlow)/outFlow > hm_iterTolerance)`. The legs keep the flows of the previous water level, so their sum is not exactly the node inflow; `hm_iterTolerance` is an absolute flow tolerance used as a relative one.

**5.4 Critical-depth bisection** (he_hydrun.cpp:2662-2666):
```
     if (Fr < 1)
      { 
       depth1 = depth;
       depth  = 0.5*(depth1 - depth0);
      }
```
It should be `0.5*(depth1 + depth0)`. The bracket stays valid, so the answer is still within about 4e-7 relative, but the iteration path differs from a correct bisection. The circular test uses `9.81*1/0.3048`; channels use `32.2`.

**5.5 Two Manning constants.** Capacity and depth use `(1.49/friction)*wetArea*pow((double) hydRadius,  (double) 0.6667)*sqrt(slope)` (he_hydrun.cpp:3824). Full-pipe head loss uses `4.66*mv_friction*mv_friction*mv_length*flow*flow` / `pow((double) mv_diameter,  (double) 5.33)` (he_data.cpp:5508-5509), which corresponds to 1.486 and 16/3 ≈ 5.333. Pi is `3.1416` in geometry and `3.1415926` in infiltration.

**5.6 Total flow where per-barrel flow is meant** (parallel pipes):
* `F_GetHeadloss(pipeIndex, linkFlow)` with `linkFlow = mv_pLinkFlow[linkIndex]` (he_runhead.cpp:1473-1474 and 1630-1637): total flow, single-barrel D.
* `hm_pLinkAdjustedVelocity[i] = HHU_GetPipeVelocity(i,  hm_pLinkAdjustedDepth[i],  pipeDiameter/2, hm_pLinkFlow[i], false);` (he_staticrun.cpp:1448).
* `mv_pWork[i]  = hm_pLinkFlow[i]/HHU_GetLinkFlow(i,  pipeDiameter);` (1391): q/Q against one barrel.
* The minor-loss velocity uses total flow (he_runhead.cpp:2313-2315).

**5.7 Backwater below the crown adds a full-pipe friction loss** (he_runhead.cpp:1635-1637):
```
   if ((downstreamHead > upstreamHead) && (downstreamHead < downRefHead))
    {
     double upHead = downstreamHead + hm_pPipe[pipeIndex].F_GetHeadloss(pipeIndex,  linkFlow);
```
This applies the full-length, full-pipe loss to a part-full pipe.

**5.8 Minor loss uses the loop counter as the pipe index** (he_runhead.cpp:2315):
```
       double  velocity             = HHU_GetPipeVelocity(i,  diameter,  radius,  linkFlow);
```
`i` is the 0-based position in the node's downstream list. Only the shape lookup is affected, which is harmless when all pipes are circular. The loss is also added for every flowing downstream pipe (summed, not max), and is applied regardless of surcharge, because the condition `//  if ((downHead >= (fromInvert + diameter))||...` is commented out (2310).

**5.9 Unbounded curve extrapolation and negative flows.** `F_GetInterpolatedY` extends the last segment without a floor (he_data.cpp:10276). A falling population curve (`POP 7`, 1.8 at 500 units) drove 70 of 87 pipe flows negative in the "pcurve" variant. Negative `Q` then flows through every routine: depth D/2, velocity 0. Pipe flows are clipped to 0 only at output (he_staticrun.cpp:1246-1247), while the outlet flow was summed before the clip and is reported negative (-7.016 cfs in that variant).

**5.10 Nested flow splits leave stale upstream heads** (he_runhead.cpp:583-623 with 1181-1195). A split node is resumed when its *level* comes up, even if one of its downstream legs runs through another split node that has not been resumed yet. Its upstream is then traversed with a head that is later raised. In the "topo2" variants (bypass CTS-480 -> CTS-460, relief CTS-600 -> OUT2, tailwater 104 ft at outlet 10) the instrumented engine logged:
```
HEADTRACE node CTS-600 via REL-600 -> 96.592850
HEADTRACE node CTS-610 via CTS-610:CTS-600 -> 96.983881
HEADTRACE node CTS-600 via CTS-600:CTS-590 -> 101.998917
```
So CTS-610's reported head (96.98 ft) sits below its downstream manhole's (102.00 ft) on a flowing, standard-slope pipe. The same mechanism makes `ftype` of a split node's first-processed leg use the partial head (3.6). The port must reproduce the DFS and level loop exactly, not compute heads in plain reverse-topological order.

**5.11 Channel geometry** (out of scope, noted):
* `double  perimeter = width + actualDepth*(sqrt(1 + leftSideZ) + sqrt(1 + rightSideZ));` (he_hydrun.cpp, `HHU_GetTrapWetPerimeter`): should be `1 + Z*Z`.
* `double  wetArea     = (2.0/3.0)*(actualDepth*width);` (`HHU_GetParabWetArea`): uses the reference width, not the width at the depth.

**5.12 Report-only: infiltration subtracted twice.** `Fv_PipeData` already did `hm_pUnpeakFlow[i] -= hm_pInfiltrateFlow[i];` (he_staticrun.cpp:1255), and the report prints `(hm_pUnpeakFlow[i]- hm_pInfiltrateFlow[i])*hm_flowUnitFactor` (1896). In the "infil" variant the `.rpt` UnPeak column for CTS-020:CTS-010 reads 17.99 while the `.out` holds 18.83 (= 17.99 + 0.84). The `.out` and IEDB values are right; only the text report is wrong.

**5.13 Thresholds in input units.** `if (fabs(baseLoad) < HMAP_ZeroTolerance)  continue;` (he_input.cpp:4252) drops loads below 0.0001 in the *input* flow unit. That is 0.0001 cfs (about 65 gpd) in a CFS model and 0.0001 gpm in a GPM model. Extra loads are not filtered.

**5.14 Adverse pipes.** Full flow is -1 per barrel (he_hydrun.cpp:3618-3619), and `HHU_GetLinkFlow` returns that for any depth (3378). As a result, any adverse pipe with flow is "full" (depth D) and q/Q = -Q/barrels: for 125_21A, `CTS-452:CTS-451` reports q/Q = -11.470, rescap 0 and critical depth 6.7e-08 ft. `Fv_CheckSplitAdversePipe` (he_input.cpp:11831) rejects an adverse pipe as a *hydraulic* split leg.

**5.15 P = 0 peaking.** `pow(0, 0) = 1`, so with Peaking Factor p = 0 every link from a source manhole carries K cfs even with no load (he_hydrun.cpp:1438). This was found by reading, not run.

**5.16 Build-copy note.** `H_StaticRunManager`, `~He_StaticRun` and `Hx_EngineInputRun` contain `MessageBox(...)` debug calls ("Static End", he_staticrun.cpp:87-92 and 172-177). They are stubbed in the Linux build and have no effect on results.

---------------------------------------------------------------------------

## 6. Verification

### 6.1 Reference implementation

`/home/claude/ise/spec/steady_ref.py` (pure Python 3, no dependencies) parses the InfoSewer input text, applies sections 0-4, and writes the same 28 arrays as the engine's `.out`. Run it as `python3 -I steady_ref.py <in.txt> <engine.out>`; add `-v` to list every differing value.

### 6.2 Inputs

`/home/claude/ise/spec/tests/` holds `base.txt` (a copy of `runs/verify2/v.txt`, 125_21A: 87 pipes, 86 loading manholes, 1 outlet, CFS). `gen_variants.py base.txt .` regenerates the 23 variants below. Each was run through the unmodified engine `/home/claude/ise/port/build/ise_run <v>.txt <v>.rpt <v>.out steady`, and all returned status 0.

| variant | what it exercises |
|---|---|
| base | 125_21A as delivered: unpeakable loads, flat outlet pipe, adverse pipe, 2-way hydraulic split, backwater, part-full flags |
| peak | every load Peakable (K·B^P) |
| cover | every load Coverage (Babbitt a-e formula) |
| mixed | types 0/1/2 mixed, multiplier 1.7, non-default K, P, a-e, extra loads incl. one below 1e-4 |
| infil | infiltration types 1-4 on every pipe |
| pct1 / pct2 / pct0 / pct100 | percent + hydraulic, two percents (sum < 100, rescaled), 0 % (stays hydraulic), 100 % (hydraulic leg gets 0) |
| minor | minor-loss K = 0.25-0.55 at every manhole |
| parallel | 1-3 barrels, including both split legs |
| locked | Manhole Head Type 1, one node Unlocked, outlet tailwater 58 ft |
| submerged | outlet tailwater 57.5 ft, Unlocked caps |
| mgd / gpm / lps | unit conversion (US and SI) |
| pcurve | flow peaking curve + population curve extrapolated to negative flow (5.9) |
| pcurve2 | in-range curves, including the "through origin" low end |
| tol | Flow Tolerance 0.05 |
| heavy | loads x6 + 0.4: 84 of 87 pipes surcharged, 81 manholes at the rim |
| zero | no load from CTS-460 up: zero-flow pipes, ByInvert heads, the 5.1 path |
| topo | second outlet with submerged head, relief split to it, two-level split |
| topo2a / topo2b | as topo with a 104 ft tailwater, both outlet orders (5.10) |

### 6.3 Result

**All 24 inputs match the engine bit for bit.** For every array (20 pipe, 6 loading-manhole, 2 outlet) the maximum absolute and relative difference is 0.0. That includes velocity, depth, d/D, q/Q, critical depth, Froude, adjusted depth and velocity, heads, and the flow-type, partial, backwater, surcharge and jump flags.

The engine itself matches InfoSewer's stored IEDB results for 125_21A (`compare_m125.py`): flow relative error ≤ 4.8e-16, velocity/depth/d/D/head absolute error ≤ 5e-10. That is the DBF storage precision. So the reference matches InfoSewer to well within 1e-6.

### 6.4 Open points (not resolved, out of the tested scope)

* Curve-based and in/out-curve flow splits (`HU_GetFlowSplitPercentFromInFlow`): not traced.
* Non-circular channels: summarised in 2.10, not tested.
* Junction chambers, wet wells, pumps, force mains, siphons, rational-method storm loads, `hm_pressureNet = 1`: excluded by the brief. `steady_ref.py` refuses such inputs.
* A loading manhole with a single *incoming* pipe and nothing downstream: the topology checks (he_topology.cpp:416-506) presumably reject it, but this was not tested.
* The original Windows build (MSVC, possibly x87) could differ from this Linux/SSE2 build in the last bits. The IEDB comparison shows no difference above 5e-10.

---------------------------------------------------------------------------

## 7. Notes for the C implementation in SWMM 5.2.4

* Hold the InfoSewer pipe and manhole data alongside SWMM's objects: rim, minor-loss K, load list (base, type, cover units), barrels, per-pipe infiltration, split type and percent, head type, outlet submerged head, and the input Flow Unit for the 1e-4 load filter. SWMM's xsect/Manning routines must **not** be used. Use the 2.1/2.2 functions verbatim: 1.49, 0.6667, 3.1416, slope 1e-5 for flat, -1 full flow for adverse.
* Keep every quantity in cfs and ft internally, as InfoSewer does, and convert at the edges.
* Flow pass: any topological order; accumulate (UP, PK, CV, CC, INF) per link, peak at each node, split with 1.6.
* Hydraulics pass: full flag, critical depth, normal depth, velocity, filter, in that order (the head pass needs the *filtered* depth).
* Head pass: reproduce `He_RunHead` literally: outlet order, `U[k]` in conn order, DFS with split suspension, level loop and per-call node adjustment (5.10).
* Validation: run `steady_ref.py` and the C code on the 24 inputs in `tests/` and compare the 28 arrays. Expect 0 difference, or at most 1-ulp noise if the C code reorders a sum.
