# InfoSewer pressure networks: wet wells, pumps, force mains, junction chambers

How the InfoSewer engine (HeEngPkg, MWH Soft / Innovyze 1996-2011) runs pump stations and force mains in steady state and EPS, written from its C++ source as the basis for `src/solver/isewer.c` and `src/solver/isewer_eps.c`. File:line references are to the engine source. Rules marked QUIRK are reproduced on purpose.

Scope: the engine's default "siphon" mode (`hm_pressureNet = 0`). The option `Advanced Forcemain Network Support : YES` (he_input.cpp:9655) switches to a gradient solver in he_pressurenet.cpp; ise2inp refuses it. Pump controls (`[Controls]`), initial pump status (`[Status]`) and flow- or time-controlled pumps are also refused.

Companion documents: `spec_steady.md` (gravity steady state) and `spec_eps.md` (gravity EPS). Everything there still applies to the gravity part of a pressure network.

## 1. Input and internal order

### 1.1 Elements

| InfoSewer section | Tokens | Notes |
|---|---|---|
| `[Junction Chambers]` (Fv_ReadChamber 4379) | ID rim diameter lossK | no loads |
| `[Wet Wells]` (Fv_ReadWetWell 4639) | ID type bottom minLevel maxLevel initLevel diameter&#124;- curve&#124;- lossK | type 0 cylinder, 1 variable area |
| `[Force Mains]` (Fv_ReadForceMain 5722) | ID from to fromInv toInv length diameter C count | count ignored (always 1) |
| `[Pumps]` (Fv_ReadPump 6718) | ID wetwell toNode type extra cap h0 h1 q1 h2 q2 speed pattern&#124;- staticFlow baseFlow | type 0 capacity, 1 design point, 2 three point |

Checks that stop the engine, all reproduced:

- Force main C must exceed 10 (5899). So the Manning branch of the force main head loss (friction < 1) is never reached in analysis.
- Wet well: 0 <= min level, 0 <= init level, 0 <= max level, min < max, min <= init <= max (4782-4825, Error 1222 for init outside the range).
- Two force mains may not leave one node, and a force main may not leave a split manhole (Fv_CheckSplitForceMain 6196).
- A station must discharge into a junction chamber, and every curve pump's shut-off head must exceed the lift from the wet well's initial water level to the siphon's tail invert, `h0 <= toInvert(tail) - (initLevel + bottom)` being Error 1441 (Fv_CheckStationData 7800-7860).
- Capacity and curve pumps may not share a station (Fv_CheckParallelPumpV2 7174).
- Pump curve coefficients outside a >= 0, 0 < b <= 20 (Fv_SetPumpCoeffData 12515).

### 1.2 Order

Nodes: loading manholes, junction chambers, outlets, wet wells, each in input order (he_input.cpp:740). Links: gravity pipes, force mains (together the "pipe mains", `hm_pipeMainCount`), then one link element per pump at `pipeMain + p`. A pump station is a topology link at `pipeMain + s` (Fv_SetStationLinkData 7577), with ID `Station_n`.

QUIRK: station links and pump elements share the index range after the pipe mains, and steady state keeps both in the same flow arrays. Pump `p`'s element and station `p`'s link are the same slot. Whichever is written last wins (2.3).

### 1.3 Stations and siphons

A station is the set of pumps with the same wet well and to-node, in pump input order (Fv_SetWorkStationData 7415).

A siphon is a chain of force mains (Fv_GetSiphonSource 6283, Fv_SetWorkSiphonData 6350, Fv_SetSiphonElementData 6480):

1. Leaders: a force main leaving a loading manhole starts a standard siphon. A force main leaving a chamber that no force main enters starts a station siphon.
2. Each leader's chain grows by any unassigned force main whose from-node is the current end, rescanning from the first force main after each addition (`j = -1`). It stops right after adding a force main that leaves a loading manhole.
3. The station siphon of a station is the siphon whose head node is the station's chamber (Fv_SetStationSiphonData 7687).

swmm5-infosewer also refuses a force main left outside every siphon (for example one leaving a wet well) and a chamber-led siphon with no station. The engine has no flow rule for either.

### 1.4 Unit conversion

Pump capacity and flows are divided by the flow factor, heads by the length factor (HD_PumpData::F_ConvertUnit he_data.cpp:7538). One-point pumps get `h0 = 1.33334*h1`, `h2 = 0`, `q2 = 2*q1` before conversion. Coefficients come after conversion:

```
b = ln((h0 - h2)/(h0 - h1)) / ln(q2/q1)
a = (h0 - h1) / (n*q1)^b             n = extra + 1 parallel units
```

So a curve pump's flow is the flow of all its parallel units, while a capacity pump's capacity is one unit's.

Force main diameter is in inches (mm in SI), divided by 12 (304.8). Wet well levels, bottom and diameter are in length units.

Wet well volumes (needed in EPS):

- Cylinder: computed at read time in input units, `area = 0.25*pi*d*d`, `V = area*level`, then divided by the volume factor (Fv_ReadWetWell 4836, F_ConvertUnit he_data.cpp:3382). pi is the macro 3.1415926. QUIRK: in SI the volume factor is `M3perFT3 = 0.028317`, not 0.3048³, so volumes and levels disagree by 5e-6 relative.
- Variable area: the curve is level (x) against volume (y). x is scaled by `1/lengthFactor`, y by `1/volumeFactor` (he_input.cpp:12218), then min, max and initial volumes are interpolated from it (Fv_SetWetWellVolumeData 12355).

The bottom reference level of a loading manhole is the lowest invert of every gravity pipe and force main touching it (F_SetBottomRefLevel keeps the minimum, he_data.cpp:1417; he_input.cpp:5260, 5963).

## 2. Steady state

### 2.1 Flow pass (Fv_StaticRun he_staticrun.cpp:344)

Nodes in topological order, as for gravity networks, with these changes:

- The single-source rule (one connecting link) applies to loading manholes only. If that link is a force main, its standard siphon is solved at once with the manhole's load (HHU_StaticTopSiphonLink he_hydrun.cpp:863).
- A node's own load counts only at loading manholes. Preceding flow components are summed over every visited connecting link (pipes, force mains, stations).
- The zero-inflow shortcut is skipped at wet wells and chambers (he_staticrun.cpp:606).
- A single downstream link is copied (pipe), solved as a siphon (force main, Fc_StaticGetSiphonFlow he_hydrun.cpp:4127) or solved as a station (Fc_StaticGetStationFlowV2 4365).

### 2.2 Force mains and siphons

Head loss of one force main (HD_ForceMainData::F_GetHeadloss he_data.cpp:6339):

```
station siphon, or C >= 1:  hL = 4.727 L Q^1.852 / (C^1.852 D^4.871)
standard siphon and C < 1:  hL = 4.66 n² L Q² / D^5.33      (unreachable, 1.1)
```

Velocity `Q/(0.25*pi*D²)`, pi = 3.1415926, and both are zero unless length, diameter and friction are positive and the force main has a siphon (F_GetHydraulicData 6443).

Siphon helpers (HD_SiphonData he_data.cpp:8082-8180): tail head = tail force main's to-invert + diameter; siphon head loss = the sum over its force mains at one flow; Newton factor = `Σ 4.727·1.852·L·Q^0.852/(C^1.852·D^4.871)`.

Standard siphon (Fc_StaticGetStandardSiphonFlow 4273): every force main carries the inflow and its components.

Station siphon (Fc_StaticGetStationSiphonFlow 4313): `stationFlow = Σ Q[pump element]`. QUIRK: that is one unit per capacity pump, because a capacity pump's element flow is per unit. The force mains get Q = unpeakable = stationFlow, with velocity and loss at stationFlow. Then the "Pump Flow Continuity" option decides:

- 0, Upstream (default): Q is overwritten with the station's inflow and the components copied from it. Velocity and head loss stay at stationFlow.
- anything else, Station: unpeakable = Q, the other components zero. QUIRK he_input.cpp:9712 `else if (option = HM_StationFlow_Station)` assigns rather than compares, so every non-zero value means Station.

### 2.3 Stations (Fc_StaticGetStationFlowV2 he_hydrun.cpp:4365)

1. All pumps are off if the wet well head (initial level + bottom) is at or below the bottom and the inflow is within the flow tolerance. A pump with speed below 1e-10 is off. An off pump's element gets Q = all components = 0 and head 0.
2. Capacity pumps that are on: `Q[element] = capacity` (one unit); `totalCap += parallel·capacity`; head = `max(0, tailHead + siphonLoss(parallel·capacity) − wellHead)`.
3. If the station's LAST pump is a curve pump, the on pumps are solved together (2.4), starting from the FIRST pump.
4. Station link components: Upstream mode copies the inflow components; Station mode sets unpeakable = totalCap (or the solved flow) and the rest to zero. QUIRK (1.2): these writes land in pump elements' slots.

### 2.4 Curve pumps (Fc_StaticGetPointPumpFlowV2 4698)

Newton on the station flow Q with the siphon:

```
Q0 = Σ_on parallel(first pump) · q1(pump)
limit = Σ_i speed_i² h0_i − (tailHead − wellHead)          (Σ over all pumps, QUIRK)
loop, at most 150 iterations, until |ΔQ/Q| <= 1e-4:
   while siphonLoss(Q) > limit: Q = Q/2
   for each on pump: f = wellHead − tailHead − loss + speed² h0
                     if f < 0: stop the pump loop
                     A' = a / speed^(b−2);  q_i = (f/A')^(1/b)
                     derivative term (1/b)(f/A')^(1/b − 1) · NewtonFactor(Q)/A'
   Newton step on Q − Σ q_i, with the engine's halving fall-backs
   new Q < 0.001: all curve pumps zero
bisection (Fc_PumpFlowBisection 5039) from [flow tolerance, q2 of the first pump] if Newton
   did not converge; 100 or more bisection steps stop the run (Error 3900)
```

Outputs: `q_i = (f/A')^(1/b)` with the loss of the last Newton evaluation, head `max(0, speed² h0 − A'·|q_i|^b)`.

### 2.5 Heads (he_runhead.cpp)

Wet well heads are fixed at initial level + bottom (he_staticrun.cpp:937). The reverse sweep from each outlet adds:

- Force main reached (Fv_GetForceMainHead 1722): once per siphon (per-run flags). For a station siphon, `down = wetWellHead + pumpHead(first pump) − siphonLoss(Q of the force main given)`; if that exceeds the head at the tail node, the tail node is re-adjusted to it and it becomes the reference. Then from the tail force main upwards: `h = ref + hL(Q)`, raised to `fromInvert + D` when flowing, or set to `fromInvert` at zero flow; the from-node is adjusted and its head becomes the next reference.
- Station reached (Fv_GetStationHead 1871): wet well head = downstream head − pump head of the station's first pump, then the wet well adjustment overrides it with the fixed level.
- Node adjustment (Fv_AdjustNodeHead 1910): loading manholes as before, but the downstream adjustment also counts a flowing force main as `fromInvert + D` (2152), and minor loss uses `Q/(3.14 r²)` for force mains (2258). Chambers take the minor loss and no rim cap. Wet wells keep their fixed head.
- An outlet fed by anything but gravity pipes stops the run (1296).

QUIRK: with zero flow, the "by invert" head reads `hm_pLoadNode[upstream node index]`, a node index used as a load-node index. It is right for loading manholes because they come first, and reads past the array for a pipe leaving a chamber or wet well.

### 2.6 Outputs (Fv_StaticOut he_staticrun.cpp:987)

After the gravity pipe and loading manhole arrays: chamber head; outlet flow and head; wet well head; force main flow, unpeakable − infiltration, peakable, coverage, infiltration, storm, velocity, head loss; pump flow, head, power `max(0, Q·H·62.4/550)`, usage (parallel count) and speed. Power, usage and speed are zero when |Q| <= 1e-7. The hydraulic jump flag looks at gravity pipes only.

## 3. EPS

### 3.1 Driver (Fv_FlowRouting he_dynamicrun.cpp:1551)

The node loop of the gravity EPS routes each downstream link: a station (Fv_CheckStationLink 2528), a force main (Fv_CheckForceMainLink 4503; a siphon is routed once, when its first force main is reached) or a gravity pipe. Each routed link writes one outflow record, at the hydraulic step H, to a scratch file.

Preceding flows are read back by each upstream link's topology order (Fv_GetPrecedeFlow 1934, Fv_LinkOrder he_topology.cpp:1525). QUIRK: a siphon writes the records of all its force mains together, so for a chain through a chamber or manhole the records of links ordered after a later force main are shifted. swmm5-infosewer keeps the write order and reads by topology order, as the engine does.

The pump time step is `min(Simulation Time Step for Pump, H)`, 60 s by default (he_input.cpp:984).

### 3.2 Wet well and station (Fv_ComputeStationFlow 2594)

For the whole run, on the pump time step:

```
V = max(0, initial volume); status_i = usage (1, or 0 if speed < 1e-10); off if speed <= 0
   status_i = 1 if V > Vmax − 1e-4; 0 if V < Vmin + 1e-4          (Fv_UpdatePumpInitStatus 2863)
loop while simTime <= end + Δt:
   inflow = linear interpolation of the preceding records at hydTime
            (QUIRK: hydro index 0 uses 0 for the previous value)
   station flow (3.3) at V
   Δt = Fv_GetPumpTimeIncrement (3.4)
   V' = V + (inflow − station flow)·Δt
   station outflow record at H steps, interpolated when a step overshoots
   at report times: pump flow, head, usage (parallel count if on and flowing), speed;
      wet well volume (capped at Vmax, the excess reported as overflow), %, level
   V = V'
```

Initial wet well storage adds the initial volumes; final storage adds the last reported volumes.

### 3.3 Station flow (Fv_GetStationFlowV2 3020)

Differs from steady state:

- A capacity pump's flow is `parallel·capacity`, and its head uses the loss at `parallel·(parallel·capacity)` (QUIRK).
- Curve pumps are solved from the station's LAST pump (Fv_GetPointPumpFlowV2 3434): start `Σ parallel(last)·q1`; `limit = min_i speed_i² h0_i − (tail − wellHead)` (minimum, not the steady sum); Newton up to 50 iterations with tolerance 1e-4; the derivative term is `(1/b)·f^(1/b−1)/A'^(1/b) · NewtonFactor`; a zero derivative steps by the flow tolerance; new Q < 1e-4 means zero; bisection from [1e-5, 5·q2(last pump)] with tolerance 1e-3 and at most 100 steps (100 stops the run, Error 3900); outputs use the loss at the final flow.
- A pump with zero speed but on gets speed 1. A zero shut-off head is replaced by `1.33·h1`.
- Wet well head = `level(min(V, Vmax)) + bottom` (HU_GetWellHeadFromVolume he_util.cpp:1589); level from volume by cylinder area `pi·r·r` or by the curve (1539, 1708).

### 3.4 Pump time increment (Fv_GetPumpTimeIncrement 4232, no controls)

```
Δt = pump step; rate = inflow − station flow; V' = max(0, V + rate·Δt)
if |rate| < 1e-4 and the wet well has inflow links: return Δt        (usage not updated)
usage = status
if (H·hydroIndex − hydTime) > 0: Δt = min(Δt, that)
if V' > Vmax − 1e-4: t = (Vmax − V)/rate; if t < Δt − 1e-4: Δt = t + 1e-4, turn pumps on
   (all pumps, or with more than one pump only those not off by speed; Fv_UpdateWellFullPumpStatus 4429)
else if V' < Vmin + 1e-4: t = −(V − Vmin)/rate; if t < Δt − 1e-4: Δt = t + 1e-4, all pumps off
if Δt < 1e-4: return the pump step
```

V' always uses the full pump step, even after Δt is shortened.

### 3.5 Siphons (Fv_GetStandardSiphonFlow 4588, Fv_GetStationSiphonFlow 4820)

Standard: the record is the sum of the head node's preceding records plus its load `pattern·mult·base` at each H. Report j uses hydro index `j·(report step/H)`. Every force main in the chain gets the same flow, its own velocity and head loss, and the same record.

Station: the record is the station's outflow record. QUIRK: report j uses the station's hydro-step value j, so it is a report-time value only when H equals the report step. Gravity pipes entering the station's chamber are ignored by the siphon.

### 3.6 Outputs (Fv_DynamicData 11228)

Per report time: loading manhole arrays; chamber head; outlet flow and head; wet well level, head, volume, % full, overflow; the 10 gravity pipe arrays; force main flow, velocity, head loss; pump flow, head, power, usage, speed (speed zero when usage is zero); then the two HGL flags over gravity pipes.

The HGL (Fv_DynamicGetHead 11352) uses the force main report flows, a station's link flow = the sum of its pump flows, the pump report heads, and wet well head = report level + bottom. The jump flag looks at gravity pipes only (14762).

Mass balance (he_dynamicrun.cpp:14855): wet well initial and final storage and wet well flooding (last report's overflow) join the gravity terms, `check = |100·(in + init − (out + final + pipe storage + flooding))/(in + init)|`. Loads at manholes that feed force mains are not counted as inflow, which is part of why Thornton's check is 14.6 %.

## 4. SWMM representation

| InfoSewer | SWMM |
|---|---|
| junction chamber | `[JUNCTIONS]` + `[INFOSEWER_CHAMBERS] Node Diameter`; its K into `[LOSSES]` of the links leaving it |
| wet well | `[STORAGE]`: invert = bottom, max depth = max level, initial depth = initial level, `CYLINDRICAL d d` or `TABULAR` area curve `WWA_<id>`; `[INFOSEWER_WETWELLS] Node MinLevel Curve&#124;* [Diameter]` |
| force main | `[CONDUITS]` + `[XSECTIONS] FORCE_MAIN D C`, `FORCE_MAIN_EQUATION H-W` |
| pump | `[PUMPS]` with curve `PC_<id>` (PUMP3 head-flow, PUMP4 for capacity), startup/shutoff = max/min level; `[INFOSEWER_PUMPS] Link Type Extra Capacity Head0 Head1 Flow1 Head2 Flow2 Speed Pattern StaticFlow BaseFlow` |
| Pump Flow Continuity | `INFOSEWER_PUMP_CONTINUITY` |
| Simulation Time Step for Pump | `INFOSEWER_PUMP_STEP` |

Pipe and force main inverts are taken as written in the `.inp`, before SWMM drops an offset below a node invert (WARNING 03), because InfoSewer keeps a pipe invert below a wet well bottom.

In the SWMM state a force main runs full, a pump carries its element flow (EPS: the station contribution, capacity pumps times parallel units), and a storage node's volume is SWMM's own volume at InfoSewer's level. In steady state each node's lateral inflow is InfoSewer's increment there, so continuity is 0.000 %. In EPS SWMM's continuity error reflects InfoSewer's own: a station siphon passes on the station flow sampled at H, not the volume pumped between samples (Thornton: InfoSewer check 14.6 %, SWMM 4.5 %).
