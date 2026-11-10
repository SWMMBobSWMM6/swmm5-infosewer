# InfoSewer EPS gravity routing: a specification for re-implementation

**Scope.** This document covers how the InfoSewer engine (Innovyze/MWH Soft HeEngPkg, C++) runs an extended period simulation ("dynamic run", `runOption 2`) on a **gravity-only sanitary network**: loading manholes, circular gravity pipes and outlets. Wet wells, pumps, force mains, siphons, junction chambers, open channels, storm runoff, water quality, sediment and source tracing are out of scope. It is written so the same routing can be rebuilt in C as a flow-routing option inside EPA SWMM 5.2.4, with results that match InfoSewer bit for bit.

**Status.** A standalone Python version of this spec, `eps_ref.py`, reproduces the engine **exactly, with zero difference in IEEE doubles**. It matches every gravity-related array in the engine's binary `.out` (all pipes, nodes and outlets at all report times), the `[Pipe Results]` and `[Maximum Gravity Pipe Flows]` tables in the `.rpt`, the `[Mass Balance Calculation]` block, and all 506,160 routing sub-steps of the 125_21A model. It was tested on the real model and on 9 variants that exercise patterns, infiltration, Muskingum–Cunge, time offsets, parallel pipes, surcharge and GPM units (section 10).

**Source references.** `file:line` refers to the build copy in `/home/claude/ise/port/build/src/`. Its line numbers are identical to the original `/home/claude/ise/src/Code/HeEngPkg/He_DynamicRun.Cpp`; the build patches only change include case and `CString` casts.

The task brief mentioned an unused `He_DynamicRun_newengine.Cpp`. It **does not exist** in `/home/claude/ise/src` (searched the whole of `/home/claude`), so there is nothing to compare.

**Units.** Everything inside the engine is US customary: ft, s, cfs. Input diameters are in inches and are divided by 12 when read. Input flows are divided by the flow-unit factor; Manning uses the constant 1.49 whatever the user's units. All times are in seconds.

---

## 1. The big picture

InfoSewer EPS is **not a network time-stepper**. It is pipe-by-pipe **hydrograph routing**:

1. Nodes are processed in upstream-to-downstream (topological) order. For each node, every pipe leaving it is routed **for the whole simulation period** before the engine moves on.
2. The inflow hydrograph of a pipe is built from:
   - the local node load (base × pattern × multiplier);
   - infiltration assigned to that pipe;
   - the stored outflow hydrographs of the upstream pipes, which were sampled at the hydraulic step `H` (= pattern time step) and are re-interpolated linearly.
3. A pipe is cut into `n` segments and routed on its own fixed sub-step `dt`, close to the full-flow travel time of one segment:
   - default: pure translation (a lag);
   - optional ("Flow Attenuation = 1"): Muskingum–Cunge.
4. The pipe's outflow is sampled at multiples of `H` into `DF[0..nHyd-1]`, with `DF[0]` forced to 0. Downstream pipes read that array.
5. At report times, flows are linearly interpolated from the sub-steps. Depth is the Manning normal depth of the *reported* flow, or full when the flow is at or above full-pipe flow.
6. After routing, a separate backward HGL sweep runs from the outlets on each report-time snapshot (`He_RunHead`). It produces node heads, backwater flags, adjusted depths and adjusted velocities.

There is **no manhole storage, no backwater feedback into the flows, no surcharge storage and no iteration**. Downstream conditions never change upstream flows.

```
read input; compute Qfull for every pipe (x parallel count)
for node N in topological order (outlets skipped):
    if N has >= 2 downstream pipes:  S[i][h] = split of (load + upstream DF + infil) at every H step
    for each downstream pipe P of N:
        (n, dx, dt) = segmentation(P)
        I[j]  = inflow series at t_j = start + (j-1)dt, j = 1..K+1
        route n segments (translation or Muskingum-Cunge)          -> C_last[k]
        DF_P[h] = C_last at hydTime = h*H (DF_P[0] = 0)            -> stored for downstream
        report samples at r*R by linear interpolation             -> (U, Dn, V, d)
        max flow over sub-steps                                    -> range table
for each report step r: HGL sweep -> node heads, flags; write .out arrays
mass balance
```

---

## 2. Time base and the routing time loop

### 2.1 Global times (he_input.cpp:820-993)

| Quantity | Formula |
|---|---|
| `startTime`, `simDuration`, `patternTimeStep`, `reportTimeStep` | seconds. `Fv_GetDecimalTime` takes a number plus a unit prefix (Sec/Min/Hour/Day, default hours), or `hh:mm:ss` |
| pattern split | if `patternTimeStep > reportTimeStep + 1e-4`: `inc = max(1, (int)(patternTimeStep/reportTimeStep + 0.9995))`, `patternTimeStep /= inc`, and every pattern value is stored `inc` times when read (he_input.cpp:957-962, 1501ff) |
| `endTime` | `startTime + simDuration` |
| `reportStepCount` | `(int)(simDuration/reportTimeStep + 1.0 + 1e-4)` |
| **`H` = hm_hydraulicDeltaT** | `= patternTimeStep` after the split (he_input.cpp:985). The pump step plays no part for gravity |
| `nHyd` = mv_hydEventCount | `(int)(simDuration/H + 1.0 + 1e-4)` (he_dynamicrun.cpp:396) |

### 2.2 Segmentation and sub-step of each pipe (Fv_ManholePipeInit, he_dynamicrun.cpp:4946-5030)

```c
n  = (int)(L/minDeltaX - 0.0001) + 1;   if (n > maxSeg) n = maxSeg;
r  = 0.5*D;   Afull = pc*3.1415926*r*r;            /* pc = parallel count; pi = 3.1415926 (he_define.h:35) */
if (L > minDeltaX) dx = L/n;  else { n = 1; dx = L; }
if (fabs(S) <= 1e-5)  dt = H;                        /* HU_IsFlatGravityMain */
else                  dt = (double)(int)(Afull*dx/Qfull);   /* Qfull = hm_pFullFlow (x pc); truncation to whole seconds */
if (dt > H)        dt = H;
if (dt < minTravel) dt = minTravel;                  /* default 1 s */
multiple = (int)(H/dt) + 1;                          /* always +1, even when H/dt is a whole number */
dt = H/multiple;
if (dt < 1.0) dt = 1.0;
```

Consequences:
- `dt` is always `H/m` for a whole number `m >= 2`, so sub-steps land on multiples of `H`. For H = 1 h, `dt <= 1800 s`.
- `dt` is a little **less than** the full-flow travel time of one segment.
- When `n` is capped by `maxSeg`, segments are longer than `minDeltaX`.

Special cases:
- **Adverse pipes**: `Qfull = -1` (he_hydrun.cpp:3618), so `(int)` gives a negative number, the minimum-travel-time clamp makes it 1, `multiple = H+1`, and the final clamp sets `dt = 1 s`. That is 86,402 sub-steps per day.
- **Very steep or short pipes**: when `Afull*dx/Qfull < 1 s`, the same chain ends at `dt = 1 s`.

There is **no Courant-based adaptivity and no step alignment with pumps or reports**. The `[0.9, 1.25]` clip exists only inside the Muskingum–Cunge coefficients (§4.3), where it has no effect.

### 2.3 The sub-step loop (Fv_TimeLoopAttenuate, he_dynamicrun.cpp:5494-5793)

The same loop runs once per segment:

```c
simTime = startTime; hydTime = 0.0; k = 0;
while (simTime <= endTime + dt) {           /* 5539 */
    k++;  nextSimTime = simTime + dt;
    ...                                     /* step k: hydTime_k = (k-1)*dt, accumulated */
    simTime = nextSimTime;  hydTime += dt;  /* 5767-5768: accumulate, never k*dt */
}
```

- `K` is the number of passes, about `simDuration/dt + 2`.
- Array index `j` stands for time `(j-1)*dt`; the arrays are used over `j = 1..K+1`.
- `simTime` is absolute and drives the patterns. `hydTime` is time since start and drives the upstream hydrographs, reports and `DF`.
- For bit-exactness, accumulate `simTime` and `hydTime` with `+= dt` exactly as the C code does.

### 2.4 Initial conditions

The run starts **empty**:
- every segment flow starts at 0, and `C[0] = C[1] = 0` always;
- every upstream hydrograph has `DF[0] = 0` (he_dynamicrun.cpp:4957), so upstream flow ramps linearly from 0 to `DF[1]` over the first `H`;
- there is no steady-state start and no hot start.

At 00:00, only pipes that leave a loaded manhole **and** have one segment report a flow (see §6.1).

### 2.5 Order of computation for one pipe

1. Segment 0, steps `k = 1..K`: build the inflow (§3.5), then the segment outflow (§4).
2. Segments `s = 1..n-1`: apply the shift (§4.2), then the steps.
3. During the **last** segment's pass only, at every step `k`, in this order: `Fv_DownstreamFlow` (§4.4), `Fv_ReportPipeData` (§6.1, including the outlet-volume bookkeeping), then `Fv_OutputPipeFlowRange` (§6.5).

### 2.6 Report sampling

Report values are **instantaneous**, linearly interpolated between the two sub-steps around the report time. They are not averages (§6.1).

---

## 3. Loads

### 3.1 Base loads (he_input.cpp:3733-4371, he_data.cpp:2071-2174)

- A `[Loading Manholes]` record holds up to 10 load entries at tokens `18+4i .. 21+4i`: base load, load type, pattern ID, coverage count.
- An entry is skipped when it is `-` or `|base| < 1e-4` (he_input.cpp:4252).
- `[Extra Loads]` (`node base type pattern cover`) adds more entries.
- A pattern of `-` means the default load pattern. With no default it becomes the internal flat pattern `ZZZZZ<<>>FLATPATTERN`, for which `HU_GetPatternItemValue` returns 1.0 (he_input.cpp:11124-11127, 11259-11267; he_util.cpp:703).
- Base loads are divided by the flow-unit factor (HD_LoadData::F_ConvertUnit).

```
load(t) = SUM_e  base_e * (mult * pat_e(t))      F_GetDynamicTimeFlow  (used for pipe inflow)
load(t) = SUM_e  pat_e(t) * mult * base_e        F_GetDynamicTotalFlow (split nodes, outputs, report)
mult    = "Global Loading Multiplier" (hm_loadMultFactor, default 1.0)
```

Both forms sum the entries in input order. Keep the two multiplication orders if you want the last bit to match.

- **EPS ignores load type and coverage count.** Unpeakable, peakable and coverage loads are all `base × pattern × mult` (he_data.cpp:2126: "Dynamic runs assume no categorisation on loading types").
- **No peaking factor is applied in EPS.** The peaking k/p/a..e options are steady-state only.

### 3.2 Pattern evaluation (HU_GetPatternItemValue, he_util.cpp:635-868)

- The per-pattern **type in `[Patterns]` decides the interpolation**: 0 = stepwise, 1 = linear.
- The `Loading Pattern Type` option sets `hm_loadPatternType`, which is **never read** (he_input.cpp:10526).
- Time is the absolute `simTime`, so a non-zero start time shifts into the pattern. Patterns wrap with their own length.

```c
item = (int)(t/ps + 1e-4);  while (item >= size) item -= size;  if (item < 0) item = 0;   /* ps = patternTimeStep (after split) */
if (stepwise) return v[item];
/* linear */
i1 = (int)(t/(inc*ps))*inc;  while (i1 >= size) i1 -= size;
i2 = i1 + 1;  if (i2 >= size) i2 = 0;
if (i2*ps < t) { ref = min(i1,i2) + size;
                 while ((t - ref*ps) > ps*inc) ref += size;
                 if ((size - i1) <= 1) ref -= 1; }
else             ref = i1;
return v[i1] + (t - ref*ps)*(v[i2]-v[i1])/(ps*inc);
```

The linear form interpolates between the original (un-split) pattern points and wraps from the last value back to the first.

### 3.3 Infiltration (HD_PipeData::F_GetInfiltrateFlow, he_data.cpp:5420; conversion 5119-5140)

`[Infiltrate Rates]` takes 5 tokens: `pipe type rate count pattern`. The rate is divided by the flow-unit factor; types 2 and 3 are also multiplied by 12 (diameter unit).

| type | flow (cfs) |
|---|---|
| 0 | 0 |
| 1 by length | rate × L |
| 2 by length × diameter | rate × D_ft × L, with rate ×12, so rate × D_in × L |
| 3 by surface | rate × 3.1415926 × D_ft × L, with rate ×12 |
| 4 by count | rate × count |
| 5 by pattern | rate × pat(t) |

Every row is then multiplied by the parallel count. The infiltration of a pipe is added to **its own** inflow at the upstream node.

### 3.4 Storm

When `Storm Runoff Analysis = 1`, a storm hydrograph (hydrograph pattern × local area / hydrograph area) is added at nodes. This is out of scope and not verified.

### 3.5 How loads become first-segment inflow

**Single downstream pipe** (he_dynamicrun.cpp:5559-5622, 5801-5857). For `j = 1..K+1`, with `t_j = simTime` of step `j` and `τ_j = t_j - startTime`:

```
I[j] = load(t_j) + infil_P(t_j) + SUM_up interp(DF_up, τ_j)
```

Summation order: load, then infiltration, then upstream pipes in input-file order. This was checked: reversing the upstream order changes up to 391 values in v1 by about 1e-14, while input order matches exactly. At step `k`, `j = k+1` uses `nextSimTime`; `I[1]` is formed only at `k = 1`.

**Split node** (two or more downstream pipes; any loading manhole with more than one outgoing pipe becomes `EN_KernelType_FlowSplit`, he_input.cpp:5680-5710): `I[j] = interp(S_i, τ_j)`, using the split series of §4.5.

**Upstream-hydrograph interpolation** (5829-5851):

```c
i = (int)(τ/H + 0.001);  if (i > nHyd-2) i = nHyd-2;
δ = τ - i*H;  if (fabs(δ) <= 0.001) δ = 0;
value = DF[i] + δ*(DF[i+1]-DF[i])/H;
```

Beyond `(nHyd-1)H` this extrapolates with the slope of the last interval. Because the 0.001 is added to `τ/H`, a time within 0.001·H *before* a boundary uses the *next* interval's slope, extrapolated backward.

---

## 4. Routing

### 4.1 Notation

- `F[j]`: inflow to the current segment.
- `C[j]`: outflow of the current segment.
- `flowTol` = `Flow Tolerance` divided by the flow-unit factor (default 1e-4 cfs).

Segment 0 starts with `F = I`.

### 4.2 Default: pure translation, `Flow Attenuation = 0` (Fv_SegmentOutflow, 5920-5990)

At each step `k = 1..K`:

```c
C[k+1] = F[k];   if (C[k+1] <= flowTol) C[k+1] = 0;        /* 5960-5965 */
```

`C[0]` and `C[1]` are never written and stay 0. After every segment except the last, the next segment's inflow is built with an **extra one-step shift** (5776-5787):

```c
for (k = 1..K)  F[k] = C[k-1];       /* F[K+1] keeps its old value; unused */
```

So, apart from the tolerance cut-off, the pipe outflow is `O[k] = C_last[k] = I[k - (2n-1)]`, which is 0 when the index is below 1. The **lag is `(2n-1)·dt`**:

- one-segment pipes: about 0.88–1.0 × the full-flow travel time;
- multi-segment pipes: about 1.7–1.8 × (125_21A median 1.71).

There is no attenuation. Flows at or below `flowTol` are set to 0 at every segment.

### 4.3 `Flow Attenuation = 1`: Muskingum–Cunge (Fv_SegmentFlowAttenuate, 5998-6145)

**Pipes with `S <= 1e-5`** (flat and adverse) use the same translation as §4.2.

**Pipes with `S > 1e-5`.** On segment 0, at every step `k <= 10`, the coefficients are recomputed from `F[k]`. They are then **frozen** for the rest of segment 0 and for all later segments of the pipe.

```c
q = (F[k] < 0.001) ? 0.001 : F[k];
q = 3.5*q;                 /* 6040 */
q = 0.5*(q + 2.5*q);       /* 6042  -> q_ref = 6.125 * F[k] */
Qf = HHU_GetFullFlow(P);   /* single pipe, not x pc */   if (Qf < flowTol) Qf = q;
if (q >= Qf || slope not standard) {
    q = 0.25*Qf;  y = ynormal(q);
    if (y >= D) { y = 0.25*D;  q = Manning(y, D/2, S, n); }
} else y = ynormal(q);
A = A(y);  T = Twidth(y);              /* circular formulas, §4.6 */
c   = 1.25*(q/A);  qw = q/T;
dT  = dx/(q/A);                        /* 6072: a local "deltaT", NOT the routing step dt */
Cr  = c*(dT/dx);                       /* 6073: = 1.25 within rounding */
if (Cr > 1.25) { Cr = 1.25; c = Cr/(dT/dx); }  else if (Cr < 0.9) { Cr = 0.9; c = Cr/(dT/dx); }
Dn  = (qw/S)/(c*dx);
c1 = (1+Cr-Dn)/(1+Cr+Dn);  c2 = (-1+Cr+Dn)/(1+Cr+Dn);  c3 = (1-Cr+Dn)/(1+Cr+Dn);
```

At each step, with `C[1] = 0`:

```c
C[k+1] = c1*F[k] + c2*F[k+1] + c3*C[k];   if (C[k+1] <= flowTol) C[k+1] = 0;
```

This is Ponce's form: `c1` multiplies the upstream value at the old time, `c2` the upstream value at the new time, and `c3` the downstream value at the old time.

- The between-segment shift of §4.2 still applies, so the extra `dt` lag per interior segment remains.
- The Courant number is effectively a fixed **1.25**, because `dT` is the travel time at the reference velocity rather than the routing `dt`. The clip never matters.
- The reference flow is `6.125 × F[k]` taken in the first 10 sub-steps. On a cold start that is the local load only.
- Manning here is `HHU_GetFlow` with the pipe's own slope.

### 4.4 Outflow hydrograph passed downstream (Fv_DownstreamFlow, 10249-10331)

This runs only on the last segment. `h` starts at 1 and `DF[0] = 0`.

```c
if (h <= nHyd) {
  if (hydTime > h*H + 1e-4)            { DF[h] = max(0, DF[h-1] + H*(O[k]-DF[h-1])/(hydTime-(h-1)*H)); h++; }
  else if (fabs(hydTime - h*H) < 1e-4) { DF[h] = O[k]; h++; }
}
```

Here `O[k] = C_last[k]`, the outflow at `hydTime_k`. Because `dt = H/m`, the exact branch always fires, so `DF[h] = O(hH)`. Only `DF[0..nHyd-1]` is stored.

**This decimation to `H` (1 h for hourly patterns) followed by linear re-interpolation is a large part of InfoSewer's EPS behaviour.**

### 4.5 Flow-split manholes (Fv_SplitFlow 2010-2231; Fv_MultiHydSplitFlow 2343-2520)

The split is only computed when the node is `F_IsKernelFlowSplit()` (1649). For `h = 0,1,…` with `t = startTime + h·H` (accumulated `t += H`, loop `do … while (t <= endTime)`):

```
Qtot = load_total(t) + SUM_up DF_up[h] + SUM_{downstream pipes} infil_i(t)
```

Note that `DF_up[h]` is used directly, without interpolation.

1. **Percentage / curve pipes** (from `[Flow Splits]`; not used by 125_21A and **not verified**):
   - `Q_i = 0.01·Qtot·percent`, or the curve output for in/out curves;
   - if every pipe is percentage-based, they are rescaled to sum to `Qtot` when off by more than 1e-5.
2. **Hydraulic pipes** (the default) share `Qh = Qtot - SUM percentage`:
   - one pipe: `Q = Qh`;
   - `Qh < flowTol`: all 0;
   - `Qh >= SUM Qfull_i`: `Q_i = Qh·Qfull_i/SUM Qfull`;
   - otherwise bisection on a common water level `WL`, starting at `0.5·(min upstream invert + max crown)`:
     - `Q_i = pc_i·Manning_i(min(WL - inv_i, D_i))` (0 when `WL <= inv_i`);
     - halve `[minLevel, maxLevel]` while `|ΣQ - Qh|/ΣQ > iterTol`, where `iterTol = min(1e-5, flowTol)`;
     - when every pipe has depth `>= D - 1e-4`, set `Q_i = Qfull_i` and stop.

The exact first-guess and update order are in `eps_ref.py:multi_hyd_split`. Each downstream pipe then reads its own `S_i` through the interpolation of §3.5.

### 4.6 Hydraulic helpers (he_hydrun.cpp), circular pipes

| Function | Exact formula |
|---|---|
| wet area (2995) | `a2 = 2·3.1416 - 2·acos((y-r)/r)`;  `A = 0.125·(a2 - sin a2)·(2r)²` |
| wet perimeter (3944) | `P = 0.5·a2·2r` (same `a2`) |
| top width | `T = 2·sqrt(y·(2r - y))` |
| Manning (3824) | `Q = (1.49/n)·A·pow(A/P, 0.6667)·sqrt(S)`, with `S = 1e-5` when `S == 0` |
| full flow (3558) | `Manning(D)`, `-1` if `S < 0`; `hm_pFullFlow = full·pc` (hydrun.cpp:1563) |
| link flow (3312) | adverse pipe: returns `hm_pFullFlow` (−1); otherwise `Manning(y)` |
| normal depth (2164) | `|Q| < 1e-4`: 0. If `Manning(D) <= Q`: `D`. Otherwise bisection: `lo = 0`, `hi = D`, `y = D/2`, `e = Q - Manning(y)`; loop `{ if e·(Q-Q(hi)) <= 0: lo = y; else if e·(Q-Q(lo)) <= 0: hi = y; y = (lo+hi)/2; e = Q - Manning(y) } while |e|/Q >= 1e-4`. Non-standard slope: `y = min(D, y + 0.5·L·(-S))`, which is unreachable for circular pipes (see §9) |
| velocity (3079) | `A(y) < 1e-5` or `Q < 1e-4`: 0; otherwise `Q/A(y)` |
| filter (3262) | `y/D < 0.01`: `y = V = 0` |
| head loss, full pipe (he_data.cpp:5506) | `4.66·n²·L·Q²/D^5.33` (`Q` is the total link flow) |

Slope is `(fromInv - toInv)/L`. Slope type: `|S| < 1e-7` is flat, `S > 0` is standard, otherwise adverse (he_data.cpp:5611). Adverse flag: `fromInv < toInv`.

### 4.7 Special cases

| Case | Behaviour |
|---|---|
| zero flow | 0 is propagated; anything `<= flowTol` is cut to 0 at every segment |
| inflow above full capacity | **no limit**: pipes carry any flow at the same lag. Depth = `D` at report time |
| L <= minDeltaX | one segment, `dx = L`, lag about 1·dt |
| flat pipe (\|S\| <= 1e-5) | `dt = H/2`, lag `(2n-1)·H/2`; Manning with `S = 1e-5` |
| adverse pipe | `Qfull = -1`, `dt = 1 s`, depth always `D`; outputs overridden (§6.1) |
| negative MC coefficients | not guarded; only the `<= flowTol` cut-off applies |

---

## 5. Manholes

- **There is no manhole storage in EPS.** `Fv_SplitManholeStorage`, `Fv_SingleManholeFlow`, `Fv_SplitManholeFlow`, `Fv_ReadManholeResultFromFile`, `Fv_CheckLoadNodeSurcharge` and `Fv_LoadNodeVolume` (15264-15634) are "Redundant Functions from V.4". They are never called; the call sites at 2132, 5631 and 11254-11259 are commented out.
- `Fv_SetManholeResult` is called with `interpolateFlag = false` (1723) and only fills zeros.
- In the HGL sweep, the surcharge volume is hard-wired to 0 (11383).
- Manhole inflow becomes pipe inflow instantly (§3.5). The diameter, rim and loss coefficient of a manhole only affect heads.

---

## 6. Results at report times

### 6.1 Pipe values (Fv_ReportPipeData 10339; Fv_CalPipeReportData 10419)

These are recorded during the last segment's pass. Define:
- `Inr[k] = I[k]` for **one-segment** pipes;
- `Inr[k] = C_last[k]` (the outflow) for **multi-segment** pipes. Segments 1 and up overwrite `mv_pInflow` at 5968.

Recording:
- At `k = 1`: record report 0 as `U = Inr[1]`, `Dn = O[1] = 0`.
- At `k >= 2`: if `hydTime >= r·R` and `r <= reportStepCount`:
  - `δ = r·R + dt - hydTime` (0 when `|δ| <= 0.001`);
  - `U = Inr[k-1] + δ(Inr[k]-Inr[k-1])/dt`;
  - `Dn = O[k-1] + δ(O[k]-O[k-1])/dt`;
  - then `r++`.
- Velocity and depth come from **U**:
  - `|U| < flowTol`: both 0;
  - otherwise `q = U/pc`; `y = (U < Qfull) ? ynormal(q) : D`; `V = vel(y, q)`; then the filter.

**What lands in the .out** for each report step (Fv_DynamicData 11228-11344; gravity-only layout in `eps_ref.py:parse_out`):

- `flow = max(0, U)`.
- **Adverse pipes**: `U` is replaced by the sum of the `U` values of the pipes entering the adverse pipe's **from** node (12464-12478). This excludes that node's own load and the adverse pipe's own routed flow. The replacement is done in place, in input order.
- `vel` = `V` (adverse pipes: `vel(min(D,y), U_adverse)`).
- `depth = min(D, y)`; `dD = depth/D`.
- `qQ = max(0, flow/Manning(D))`, set to `max(qQ, 1)` when `dD >= 1`. For an adverse pipe `Manning(D)` = −1, so 0 before that rule.
- Froude, computed from the **raw** `V` and `y`:
  - `T = 2·sqrt(y(D-y))`, `A = A(y)`;
  - `T <= 0` and `y > 0`: `V/sqrt(32.2·y)`;
  - `A > 0` and `T > 0`: `V/sqrt(32.2·A/T)`;
  - otherwise 0.
- `overflow` = 0 (open channels only).
- `backwater`, `adjdepth`, `adjvel`: see §6.4.
- `flowtype`, `partial`: see §6.3.

Node arrays:
- `load = load_total(start + r·R)`; storm arrays = 0;
- `head` = §6.3;
- `status = (head >= rim)`;
- `jump = 1` if some upstream pipe has Froude > 1 and some downstream pipe has Froude < 1.

Outlet arrays: `outflow = SUM Dn` over the outlet pipes; `outhead` = §6.3.

**.rpt `[Pipe Results]`.** The "D/S Flow" column prints `U`, the *upstream* field, which is the true downstream value only for multi-segment pipes. It is the raw `U` (adverse pipes not replaced), times the flow-unit factor. The "d/D" column prints an unset Froude field (§9).

### 6.2 Depth

The depth is the normal depth of the reported flow, with no backwater. For multi-segment pipes it is the normal depth of the **outflow**, even though the engine labels it "upstream depth". Backwater appears only in the HGL and in the "adjusted" arrays.

### 6.3 Node heads: the HGL sweep (He_RunHead::F_Run, he_runhead.cpp:554-629)

Inputs: for every pipe, `y` (the reported depth) and `Q = U`, both raw.

**Defaults.** `hm_optHeadZeroLinkFlow = ByInvert` and `hm_optHeadUpstreamNode = Adjust`. Neither can be changed from the input file. The head type comes from `Manhole Head Type` (0 unlocked, the default; 1 locked) or from `[Load Node Head Type]`.

**Start at each outlet.**
- `head = max(submergeHead, max over incoming pipes of y + toInv)`. `submergeHead` comes from `[Outlet Heads]`, default −1e37.
- Then a depth-first search upstream (`Fv_DepthSearch`). Split nodes are suspended and later resumed level by level.

**For each pipe P, from upstream node U to downstream node N** (with `hN` = head of `N`), `Fv_GetPipeHead` (1404):

- **Standard slope, `Q < Qfull`:**
  - `hU = (Q <= 1e-4) ? max(hN, bottomRef(U)) : y + fromInv`. `bottomRef(U)` is the lowest invert among all pipes connected to U.
  - Then `Fv_AdjustPipeHead` (1591), which applies only when `0 < Q < Qfull - 1e-4`:
    1. if `hN > hU` and `hN < toInv + D`: `hU = max(hU, hN + hl(Q))`;
    2. if `hN > toInv + D + 1e-7`: `refL = (hN - toInv - y)/sin(atan S)`; `hU = max(hU, hN + (refL > L ? hl(Q) : 0))`.
- **Standard slope, `Q >= Qfull`:** `hU = max(hN, D + toInv) + hl(Q)`.
- **Flat or adverse, `Q > 1e-4`:**
  - `dd = max(hN, y + toInv)`;
  - `hU = max(fromInv + y, dd + (y >= D ? hl(Q) : 0))`;
  - then `Fv_AdjustPipeHead(P, dd, hU)`. This returns immediately for adverse pipes because `Qfull = -1`.
- **Flat or adverse, zero flow:** `hU = fromInv`.

**Node adjustment** (`Fv_AdjustNodeHead`, 1910):
1. For each downstream pipe of `U` carrying flow, `hU = max(hU, y + fromInv)`.
2. Add the minor loss `|K|·v²/(2·32.174)`, with `v = Q/A(D)`.
3. A split node takes the max over its downstream branches.
4. Output head = `min(h, rim)` when unlocked. The stored head is that capped value, so the cap propagates upstream.

**Flags:**
- `flowtype = 1` if `head(U) >= fromInv + D + 1e-7`, otherwise −1;
- `partial` = set inside `AdjustPipeHead`, or for zero-flow flat/adverse pipes when `max(hN, toInv) > toInv`; written as 1 or −1.

### 6.4 Backwater adjustment (Fv_PipeBackwaterAdjustment, 12760)

This uses the node heads and the **adverse-modified** flows and velocities.

- If `head(to) > y + toInv`:
  - `adjdepth = max(y, min(D, 0.5·((head(to) - toInv) + (head(from) - fromInv))))`;
  - `adjvel = Q/A(D)`, a full-pipe area even when not full;
  - `backwater = 1`.
- Otherwise `adjdepth = y`, `adjvel = vel`, `backwater = 0`.

### 6.5 Maximum table (Fv_OutputPipeFlowRange 10544; Fv_FlowRange 10879)

- Max flow is taken over **every sub-step** of `C_last[k]`. It starts with the value at `hydTime = 0` and is updated only when it exceeds the stored max by more than `flowTol`. The time recorded is the absolute `simTime`.
- Depth = `ynormal(Qmax/pc)`, which returns `D` at or above full flow. Velocity = `vel(depth, Qmax/pc)`. Then the filter; `d/D = depth/D`.
- The table time is printed by `HU_GetMilTime`: `t = h + 0.0002777`, hours and minutes truncated.
- The critical-depth column was not re-implemented.

---

## 7. Mass balance (Fv_DynamicReport, 14854-14866)

```
in   = SUM over single-outlet nodes [ (load+infil)(t_1)·dt + SUM_{k=2..K} (load+infil)(t_{k+1})·dt ]   (5583/5595/5607/5613)
     + SUM over split nodes SUM_{h=0..} [ load_total(t_h)+SUM infil(t_h) ]·H                        (2112, 2128)
out  = SUM over outlet pipes [ O[1]·dt + SUM_{k=2..K} O[k+1]·dt ]                                    (10350-10356)
stor = SUM pipes A(y_last)·pc·L,  y_last = report depth at the last report step (multi-segment: depth of the outflow)
       accumulated in Fv_ReadPipeMainResultFromFile, which runs once for the .out and again for the .rpt link table,
       then printed as mv_finalPipeStorage/2.0
flood (open conduits) = internIn - internOut = 0 for closed pipes
check = |100·(in + initWW - (out + flood + finalWW + stor + wwOverflow))/(in + initWW)|
```

- The integration runs over `K` sub-steps from `start` to `end + dt`. It skips `t_2` and includes `t_{K+1}`.
- Split nodes count `nHyd·H` seconds of load (90,000 s for a 24 h run).
- In 125_21A, `in` = 1,664,870.666 ft³ against 1,626,951.571 ft³ actually loaded in 24 h (+2.33 %). Of the 37,919 ft³ excess, 37,651 is split node CTS-460 counted one hour too long.
- The printed 4.894 % is therefore mostly a bookkeeping artefact.

**The uninitialised-variable defect (open channels only), he_dynamicrun.cpp:5676-5720:**

```c
        if (HU_CheckChannelOpenChannel(pipeIndex))
         {
			 double nodeoutflow;      //Feng Shang 05/19/2011
			 double nodeinflow;       //Feng Shang 05/19/2011
          if (mv_subTimeIndex <=1) 
           {
            if (mv_pInflow[mv_subTimeIndex] > hm_pFullFlow[pipeIndex])
             { ...  nodeinflow = mv_pInflow[mv_subTimeIndex]; ... nodeoutflow = mv_pInflow[mv_subTimeIndex]; }
           }
          else if (mv_subTimeIndex > 1)
           {
			 nodeinflow = mv_pInflow[mv_subTimeIndex]; //Feng Shang 05/19/2011
            ...
			nodeoutflow = mv_pInflow[mv_subTimeIndex+1]; //Feng Shang 05/19/2011
		  }
		  mv_internOutflow += nodeoutflow*mv_deltaT;  //flow leaving a node not the system //Feng Shang 05/19/2011
		  mv_internInflow += nodeinflow*mv_deltaT;  //Feng Shang 05/19/2011
		}
```

- At `k = 1` with no overflow, both variables are read uninitialised.
- For `k > 1`, `nodeinflow` is taken at step `k` while `nodeoutflow` is taken at `k+1` after capping, so "flooding" mixes two time levels.
- It never runs for circular pipes (`HU_CheckChannelOpenChannel` is false). For a gravity-only SWMM port, flooding = 0.

---

## 8. Inputs that affect EPS gravity routing

| Input (text file) | Variable | Default | Effect |
|---|---|---|---|
| `[Options] Flow Unit` | hm_flowUnitFactor | CFS = 1 | Loads, infiltration and `flowTol` are divided by it. GPM 448.831169, MGD 0.64632, IMGD 0.5382, AFD 1.9837 (he_define.h:40-48). SI units also convert lengths |
| `Global Loading Multiplier` | hm_loadMultFactor | 1.0 | scales every load |
| `Flow Tolerance` | hm_flowTolerance | 1e-4 cfs (must be > 0) | zero cut-off per segment, empty-pipe test, range hysteresis; `iterTol = min(1e-5, flowTol)` |
| `Flow Attenuation` | hm_flowAttenuation | 0 (off) | `>= 1e-4` turns Muskingum–Cunge on |
| `Minimum Pipe Segment` | hm_minDeltaX | 100 ft | segment count |
| `Maximum Number of Segments` | hm_maxSegmentNumber | 100 (a present value < 1 gives 1000) | segment cap |
| `Minimum Travel Time` | hm_minTravelTime | 1 s (values < 1 become 1) | lower clamp on `dt` before the `multiple` step |
| `Manhole Head Type` | hm_optLoadNodeHeadType | 0 unlocked | rim cap on heads |
| `Loading Pattern Type` | hm_loadPatternType | — | **no effect** |
| `Storm Runoff Analysis` | hm_stormRunoffAnalysis | 0 | out of scope |
| peaking k/p/a–e, design options | — | — | no effect in EPS |
| `[Time]` start, duration, pattern step, report step | §2.1 | 0 h, 24 h, 1 h, 1 h | time base, `H`, report times |
| `[Patterns]` `ID type v…` | per-pattern type | — | 0 stepwise / 1 linear; multi-line patterns append |
| `[Loading Manholes]` rim, diameter, loss K, loads | §3.1 | — | loads; rim/K affect heads only |
| `[Extra Loads]` | more load entries | — | same as above |
| `[Outlets]` ID rim diameter | — | — | |
| `[Outlet Heads]` | submergeHead | −1e37 | outlet HGL |
| `[Pipes]` ID from to invU invD L D(in) n parallel | §4.6 | — | parallel = extra pipes, so `pc = token + 1` |
| `[Infiltrate Rates]` | §3.3 | none | added to the pipe's inflow |
| `[Flow Splits]` | split type/percent/curve | hydraulic | §4.5 (percentage path not verified) |
| `[Report] Link : No` | hm_linkReport | yes | **halves the reported final storage** (§9 item 8) |

---

## 9. Oddities and apparent bugs (reproduce them; do not fix)

1. **The "upstream" flow, velocity and depth of a multi-segment pipe are at its downstream end.** For one-segment pipes they are at the upstream end. In he_dynamicrun.cpp:5968 `mv_pInflow[timeCounter] = mv_pCurrentSegmentFlow[timeCounter]; // MKM on 02/04/2008` overwrites the inflow on every segment pass, and only segment 0 restores it (5976-5980). This is why, at 00:00, `CTS-800:CTS-790` (1 segment) reports 0.17 while `CTS-820:CTS-810` (5 segments) reports 0.00.

2. **Extra one-step lag per interior segment.** At 5784, `mv_pFormerSegmentFlow[mv_subTimeIndex] = mv_pCurrentSegmentFlow[mv_subTimeIndex-1];` shifts again after 5960 `mv_pCurrentSegmentFlow[timeCounter+1] = mv_pFormerSegmentFlow[timeCounter];`. Total lag is `(2n-1)·dt`, about 1.7× the full-flow travel time.

3. **Upstream flow at t = 0 is dropped.** At 4956-4957, `mv_hydroIndex = 1;` and `mv_pDownstreamFlow[mv_hydroIndex-1] = 0.0;`. Every inter-pipe hydrograph starts at 0 and ramps over the first `H`.

4. **Outflow hydrographs are passed downstream only at the pattern step `H`** (10281-10285) and re-interpolated linearly (5831-5851). Peaks shorter than `H` are clipped. `ZERO_TIME_VALUE = 0.001` is added to `τ/H` (relative), but is used as an absolute 0.001 s elsewhere (10384).

5. **`dt` is truncated to whole seconds and then reduced by one extra multiple:** 4989 `mv_deltaT = (double) ((int) ((linkArea*mv_deltaX)/hm_pFullFlow[pipeIndex]));` and 5002 `double  multiple = int(hm_hydraulicDeltaT/mv_deltaT) + 1;`. The flat-pipe test (`|S| <= 1e-5`, he_util.cpp:269) differs from the flat slope type (`|S| < 1e-7`).

6. **Adverse pipes.** he_hydrun.cpp:3618-3619 `if(slope  < 0)    //If else added Feng Shang 10/26/2010` / `flow = -1.0;` makes `Qfull = -1`. As a result:
   - `dt = 1 s`;
   - depth is always `D` (10461 `if (mv_pLinkResult[reportIndex] < hm_pFullFlow[pipeIndex])` is never true);
   - the .out flow is replaced by the sum of the incoming pipes' flows at the from node (12470-12477, `flow += mv_pPipeMainResult[linkIndex].F_GetUpstreamFlow();`), which leaves out the node's own load;
   - the `y + 0.5·L·(-S)` correction in `HHU_GetFlowDepth` (2299-2305) is dead code for circular pipes.

7. **Muskingum–Cunge:**
   - 6040 `flow =  3.5*flow;`, 6042 `flow = 0.5*(flow + 2.5 * flow);` give a reference flow of 6.125 × the inflow;
   - 6072 `deltaT        = mv_deltaX/(flow/flowArea);` makes 6073 `C             = celerity * (deltaT/mv_deltaX);` equal to 1.25 always, so the `[0.9, 1.25]` clip (6076-6085) is inert;
   - the coefficients are computed only while `mv_segmentNumber <= 0` and `timeCounter <= 10` (6030-6032), i.e. from the first 10 sub-steps of a cold start.

8. **Final pipe storage is summed twice and printed halved.** 14741 `mv_finalPipeStorage += area*hm_pLink[pipeMainIndex].F_GetParallelCount()*pipelength;` runs in `Fv_ReadPipeMainResultFromFile`, which is called by both `Fv_DynamicData` (11255) and `Fv_LinkReport` (15014). 14862 then prints `mv_finalPipeStorage/2.0`.
   - Verified: with `[Report] Link : No`, storage drops from 71,094.408 to 35,547.204 ft³ and the check moves from 4.894 % to 7.029 %.
   - Also, inflow integrates `K` sub-steps (to `end + dt`), and split nodes integrate `nHyd·H` (§7).

9. **The .rpt pipe table prints the wrong fields.** 15037-15046 has 7 conversions for 8 arguments: `"%-20s  %7.4f  %6.2f  %9.2f  %9.2f  %9.2f  %9.2f\n"` with `... F_GetUpstreamFlow()*hm_flowUnitFactor, ... F_GetUpstreamVelocity()*..., ... F_GetUpstreamDepth()*..., mv_pPipeMainResult[i].F_GetUpstreamFroude(), mv_pPipeMainResult[i].F_GetUpstreamDoverD());`. The "d/D" column shows `F_GetUpstreamFroude()`, which is never set (always 0.00), and "D/S Flow" is the upstream field (item 1).

10. **Minor loss uses the loop counter as the pipe index.** he_runhead.cpp:2315 `double  velocity             = HHU_GetPipeVelocity(i,  diameter,  radius,  linkFlow);`. It is harmless for circular pipes because the area depends only on the diameter passed in.

11. **Backwater-adjusted velocity always uses the full-pipe area.** 12792 `mv_pLinkAdjustedVelocity[i] = HHU_GetPipeVelocity(i,  diameter,  radius,  linkFlow,  false);`

12. **Every pipe gets the global running max Froude.** 12671 `mv_pPipeRange[i].F_SetMaxFroude(mv_maxPipeFroude);`

13. **Constants:**
    - `pi` is 3.1415926 (he_define.h:35) for `dt` and infiltration, but 3.1416 in the wet-area, perimeter and width formulas (he_hydrun.cpp:2996);
    - Manning is `(1.49/friction)*wetArea*pow(hydRadius, 0.6667)*sqrt(slope)` (3824), not 1.486 and 2/3;
    - head loss is `4.66·n²·L·Q²/D^5.33`.

14. **Loads:** EPS ignores load type, coverage count and peaking, and the `Loading Pattern Type` option is dead (he_input.cpp:10526-10527).

15. **Hydraulic split:** the iteration uses `hm_pPipe[pipeIndex].F_GetFromInvert()` (2478) instead of `pUpstreamInvert[i]` (the same thing for gravity pipes). The convergence test `fabs(outflow-hydCalcFlow)/outflow` (2465) divides by a possibly zero outflow.

16. **Unreachable report index:** the report guard `(mv_reportIndex <= hm_reportStepCount)` (10380) allows one index past the last report; it is never reached.

17. **Unused interpolation branch:** `Fv_DownstreamFlow` would interpolate from the *stored* `DF[h-1]` rather than the previous sub-step (10272-10274). The branch is never used because `dt` divides `H`.

18. **HGL quirks:**
    - an empty pipe carries the downstream head upstream (`max(hN, bottomRef)`);
    - the rim cap propagates upstream;
    - heads use the reported (outflow-end) depth added to the *upstream* invert for multi-segment pipes.

---

## 10. Verification

**Method.**
- `eps_ref.py` reads the InfoSewer text input and runs §2–§7. It compares against the engine's binary `.out` (parsed by locating the time block, then the gravity-only layout of §6.1) and against the `.rpt` tables.
- The stock engine (`/home/claude/ise/port/build/ise_run`) produced all golden outputs.
- A private copy instrumented with trace prints, `/home/claude/ise/specwork/src` built into `/home/claude/ise/specwork/ise_run_tr` (its `.rpt` is byte-identical to the stock one), was used to compare every routing sub-step.

**Results.** In every row below the maximum absolute and relative difference is exactly 0.0 for all compared values. "All .out arrays" means pipe flow, velocity, depth, d/D, q/Q, Froude, backwater, adjusted depth and velocity, flow type, partial flag, overflow; node head, load, status, jump; outlet flow and head; and full flow.

| Run | What it exercises | Values compared | Result |
|---|---|---|---|
| `runs/eps2/e.txt` (125_21A) | constant loads, cold start, hydraulic split at CTS-460, flat outlet pipe, adverse pipe, `dt = 1 s` pipes | 19,662 .out values; 506,160 sub-steps (time, inflow, outflow); 2,175 rpt flows; 87-row max table; mass balance (4.894 %) | **all identical** |
| `variants/v1` | 24-h linear diurnal pattern, multiplier 1.25, infiltration types 1–5, 48 h | 38,454 values + rpt + max table + mass balance | identical |
| `v2` | v1 + Flow Attenuation (Muskingum–Cunge) | 38,454 | identical |
| `v3` | stepwise pattern, start 6 h, 30 h run, 30-min report (pattern split `inc = 2`, `H = 1800`), multiplier 0.8, 2 parallel pipes, min segment 50 ft, max 3 segments, min travel 30 s | 47,850 | identical |
| `v4` | v3 + Muskingum–Cunge | 47,850 | identical |
| `v5` | multiplier 3 (surcharged trunk, split above capacity), 2-h report | 10,266 | identical |
| `v6` | v2 with 40-min pattern step, start 1 h (`H = 2400`, report times between sub-steps) | 38,454 | identical |
| `v7` | v1 with 50-min pattern step, 2-h report | 19,662 | identical |
| `v8` | 125_21A with `[Report] Link : No` | 17,487 | identical (storage halving reproduced) |
| `v9` | v1 in GPM units | 38,454 | identical |

**Reproduce:**

```
cd /home/claude/ise/spec
python3 -I eps_ref.py /home/claude/ise/runs/eps2/e.txt --out /home/claude/ise/runs/eps2/e.out --rpt /home/claude/ise/runs/eps2/e.rpt
python3 -I eps_ref.py variants/v2.txt --out variants/v2.out --rpt variants/v2.rpt
python3 -I eps_ref.py <in.txt> --dump ref.csv          # per-report pipe results
python3 -I eps_ref.py <in.txt> --trace ref_trace.csv   # per-sub-step inflow/outflow (same rows as the engine trace)
```

Each run takes about 1 s for 125_21A.

**`eps_trace.csv`** (from the instrumented real engine, 125_21A). Columns:

`time_s, pipe_id, role, substep, dt_s, inflow_cfs, outflow_cfs, reported_flow_cfs, depth_ft, velocity_fps`

- `time_s` = `hydTime` (seconds since start).
- `inflow_cfs` = first-segment inflow `I[k]`.
- `outflow_cfs` = last-segment outflow `C_last[k]`.
- `reported_flow_cfs` = `Inr[k]`, the value the engine would report at that instant.
- `depth_ft` and `velocity_fps` are the engine's report-time formulas (§6.1) applied to `reported_flow`.

Three pipes, 5,885 rows:
- `CTS-820:CTS-810`, "upstream": the most upstream pipe that carries flow. `CTS-830:CTS-820` is above it but CTS-830 has no load, so it is zero throughout. 5 segments, dt 21.95 s.
- `CTS-460:CTS-453`, "mid": just below the hydraulic split. 4 segments, dt 45.57 s.
- `13`, "outlet_pipe": `CTS-010` to outlet `10`. Flat, 1 segment, dt 1800 s.

---

## 11. Notes for the SWMM 5.2.4 implementation

- **Structure.** The algorithm needs each upstream hydrograph at the *end* of the current `H` interval (`DF[h+1]`), which is up to `H` ahead of a downstream pipe's current time. A plain time-marching `flowrout_execute` therefore cannot reproduce it exactly. The straightforward exact route:
  1. at routing start, evaluate every node's external inflow on the InfoSewer time grid. For sanitary DWF these are known in advance; use InfoSewer's pattern evaluator, §3.2, not SWMM's stepwise DWF patterns;
  2. run §2–§7 for the whole horizon, in topological order;
  3. at each SWMM step, interpolate link flows from the stored sub-step series of each conduit (`Inr` or `O`, depending on which end you report), and node depths from the report-time HGL;
  4. if SWMM-generated inflows (RDII, runoff) must be included, use a lookahead of one `H`, or accept that results will differ.
- **Bit-exactness checklist:**
  - double precision everywhere;
  - accumulate times with `+= dt`;
  - reproduce the `(int)` truncations;
  - use 3.1415926 vs 3.1416, 1.49 and 0.6667;
  - use `flowTol` (user value) for routing cut-offs, the empty test and range hysteresis, but the compile-time 1e-4 (`HMAP_FlowTolerance`) inside normal depth, velocity and the HGL;
  - keep the bisection update order of normal depth exactly;
  - sum upstream pipes in input order.

  The reference matched the Linux gcc build (SSE2 doubles) bit for bit. The original Windows DLL may differ in the last bits if it was built for x87.
- InfoSewer pipes, loading manholes and outlets map one-to-one onto SWMM conduits, junctions and outfalls. For depth, use `y` (§6.1) or `head - invert` (§6.3) as required.

---

## 12. Open questions / not verified

- Percentage, curve and in/out `[Flow Splits]`; junction chambers; open-channel shapes, including their overflow logic and the uninitialised-variable defect; storm runoff; SI unit sets. All are described from the code only.
- Whether the shipped Windows InfoSewer (MSVC) gives identical last-bit results. A genuine InfoSewer EPS `.out` from Windows for 125_21A would settle it.
- `He_DynamicRun_newengine.Cpp` is not in the source tree provided.
- The critical-depth column of the max table (`HHU_GetCriticalDepth`) was not re-implemented.

---

*Disclosure: this specification, `eps_ref.py`, the variant inputs, the instrumentation of the private engine copy and `eps_trace.csv` were produced by Claude (AI) for SWMM Enablement on 2026-11-09, from the HeEngPkg C++ source. All numbers were checked by running the engine; nothing was edited under `/home/claude/ise/port` or `/home/claude/ise/src`.*
