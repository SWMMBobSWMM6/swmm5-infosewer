# swmm5-infosewer

InfoSewer's gravity-sewer engine, rebuilt as a flow routing option inside EPA SWMM 5.2.4. Set `FLOW_ROUTING INFOSEWER` in a normal `.inp` and SWMM runs InfoSewer's steady-state or extended period (EPS) calculations: loads by type, peaking, Manning normal depth, InfoSewer's hydraulic flow splits, its pipe-by-pipe hydrograph routing and its HGL sweep.

The results match the InfoSewer engine (HeEngPkg, MWH Soft / Innovyze, 1996-2011) bit for bit. On 47 test models, 1,276,720 of 1,277,022 output values are identical in every IEEE double. The other 302 are all in the one model that uses LPS, and differ by at most 7.9e-12 (see Verification).

InfoSewer left the market with ArcMap. Its models did not have to. Stock SWMM 5.2.4 refuses the 125_21A test network outright under KINWAVE (`ERROR 133` at the split manhole, `ERROR 115` at the adverse pipe). swmm5-infosewer runs it and gives InfoSewer's numbers.

## Quick start

```bash
bash build.sh                       # bin/swmm5_is, plus bin/swmm5_epa (unmodified EPA 5.2.4) for the null test
python3 -I tools/iedb2ise.py MODEL.IEDB model.txt          # InfoSewer database -> engine text input
python3 -I tools/ise2inp.py model.txt model.inp --mode eps # or --mode steady
bin/swmm5_is model.inp model.rpt model.out
bash tests/run_all.sh               # build and run every verification (about 30 s)
```

`ise2inp.py --stock DYNWAVE` writes the same network for stock SWMM 5.2.4, with the loads as `[DWF]` and `[INFLOWS]`, so the two engines can be compared on one network.

## The .inp

### Options

| Option | Default | InfoSewer option |
|---|---|---|
| `FLOW_ROUTING INFOSEWER` | | turns it on |
| `INFOSEWER_ANALYSIS STEADY` or `EPS` | STEADY | run type |
| `INFOSEWER_LOAD_MULT` | 1.0 | Global Loading Multiplier |
| `INFOSEWER_PEAK_K`, `INFOSEWER_PEAK_P` | 2.4, 0.89 | peakable loads, Q = K·Qbase^P (cfs) |
| `INFOSEWER_PEAK_A` ... `_E` | 5, 0, 0.2, 0, 1 | coverage-based peaking |
| `INFOSEWER_PEAK_CURVE`, `INFOSEWER_COVER_CURVE` | none | peaking factor curve, coverage based peaking curve |
| `INFOSEWER_FLOW_TOLERANCE` | 0.0001 (flow units) | Flow Tolerance |
| `INFOSEWER_PATTERN_STEP` | 1:00:00 | Pattern Time Step |
| `INFOSEWER_MIN_SEGMENT`, `INFOSEWER_MAX_SEGMENTS` | 100 ft, 100 | EPS pipe segmentation |
| `INFOSEWER_MIN_TRAVEL_TIME` | 1 s | EPS minimum sub-step |
| `INFOSEWER_FLOW_ATTENUATION` | 0 | 1 = Muskingum-Cunge instead of translation (EPS) |
| `INFOSEWER_RESULTS_FILE` | none | CSV of every InfoSewer output array at full precision |

### InfoSewer-only data

```
[INFOSEWER_LOADS]          Node  Base  UNPEAK|PEAK|COVER  Pattern|*  CoverUnits
[INFOSEWER_PATTERNS]       Name  STEPWISE|LINEAR  f1 f2 ...
[INFOSEWER_CURVES]         Name  X  Y
[INFOSEWER_SPLITS]         Link  Percent
[INFOSEWER_INFILTRATION]   Link  LENGTH|DIAMETER|SURFACE|COUNT|PATTERN  Rate  [Count|Pattern]
```

### The network is plain SWMM

| InfoSewer | SWMM |
|---|---|
| loading manhole | `[JUNCTIONS]`; rim = invert + max depth |
| head type Locked | junction surcharge depth > 0 (head capped at rim + surcharge depth; ise2inp writes 1000) |
| manhole minor loss K | `[LOSSES]` inlet loss of each conduit leaving the manhole |
| outlet, outlet head | `[OUTFALLS]` FREE, or FIXED at that head |
| gravity pipe, parallel count | `[CONDUITS]` + `[XSECTIONS] CIRCULAR`, barrels |

ise2inp writes `LINK_OFFSETS ELEVATION` and every number with Python's shortest round-trip form, so SWMM reads back the doubles InfoSewer uses. Junction and conduit order follow the InfoSewer input, which keeps InfoSewer's topological order and therefore its results.

## What SWMM reports

**Steady state.** The InfoSewer solution is computed once and held for the whole run. Each manhole gets a lateral inflow equal to InfoSewer's flow increment there (outflow minus inflow). That increment is negative at a confluence where peaking lowers the combined flow, which is how InfoSewer's peaked flows are not conserved. With it, SWMM's routing continuity error is 0.000 %. The `.rpt` gets two InfoSewer tables: pipe flows by load type with velocity, d/D, q/Q and surcharge, and manhole loads and HGL.

**EPS.** InfoSewer EPS is not a network time-stepper. It routes each pipe for the whole run, upstream first, on its own sub-step. So swmm5-infosewer routes the whole run at start-up, then each SWMM step reads the stored results: exactly at report times, linearly interpolated between them. A conduit's SWMM volume is the water in transit (inlet minus outlet flow, integrated), because that is what InfoSewer's translation routing holds. On 125_21A SWMM's continuity error is -0.033 %; on the diurnal variant v1 it is -0.008 %. The `.rpt` adds InfoSewer's maximum pipe flow table and its own mass balance as InfoSewer computes it.

In both modes, link flow, depth and velocity and node head in SWMM's `.out` equal InfoSewer's to float32 precision (6e-8).

## Verification

`tests/run_all.sh` runs three checks.

| Check | Models | Values | Result |
|---|---|---|---|
| steady vs InfoSewer engine `.out` (`tests/steady/`) | 32 | 94,752 | 31 identical; LPS model within 7.9e-12 |
| EPS vs InfoSewer engine `.out`, every report time (`tests/eps/`) | 15 | 1,182,270 | all identical; InfoSewer's mass balance block identical |
| null test vs EPA SWMM 5.2.4 (`tests/null/`) | 14 | | `.out` byte-identical, `.rpt` identical except clock lines |

The steady set covers peakable, coverage and mixed loads, curves, multiplier, extra loads, infiltration types 1-5, percent and hydraulic splits, nested splits, minor losses, parallel pipes, locked heads, tailwater, heavy surcharge, zero flow, CFS, GPM, MGD and LPS, and two gravity catchments cut from the 2009 Thornton model, upstream of outlets 1444 and 4723 (98 and 938 pipes). The EPS set covers 125_21A, linear and stepwise patterns, pattern infiltration, Muskingum-Cunge, start-time offsets, 30-minute and 2-hour reporting, 40- and 50-minute pattern steps, parallel pipes, segment and travel-time limits, a surcharged split, percent splits and the two Thornton catchments.

The reference `.out` files come from the InfoSewer engine compiled from its own C++ source on Linux. That build matches the results InfoSewer stored in the 125_21A `.IEDB` to 5e-10 (the DBF storage precision). SI models agree to 1e-11 rather than exactly, because SWMM stores a diameter in metres and InfoSewer converts from millimetres.

The full InfoSewer behaviour is written up in `docs/spec_steady.md` and `docs/spec_eps.md`, with a file:line reference to the engine source for every rule, and reproduced in Python in `docs/reference/`.

## InfoSewer quirks kept on purpose

Bit-for-bit means the engine's oddities come along. Each one is marked `QUIRK` in the C. The main ones:

- Manning uses 1.49 and exponent 0.6667, areas use pi = 3.1416, but full-pipe head loss uses 4.66 (that is 1.486), and infiltration by surface uses 3.1415926.
- The critical-depth bisection has a sign slip (`0.5*(d1 - d0)`), accurate to about 4e-7 anyway.
- Head loss, minor-loss velocity and q/Q use the total flow of parallel barrels where per-barrel flow is meant.
- Backwater below the crown adds full-length, full-pipe friction loss to a part-full pipe.
- Nested flow splits can leave a stale HGL upstream of a split.
- An adverse pipe has full flow -1, so it always reports full, and in EPS its flow is replaced by the sum of the pipes entering its upstream manhole.
- EPS reports a multi-segment pipe's "upstream" flow, depth and velocity at its downstream end, adds an extra one-step lag between segments, and passes hydrographs downstream decimated to the pattern step.
- In EPS, Muskingum-Cunge's Courant number is always 1.25, and load type and peaking are ignored.
- InfoSewer's own EPS mass balance counts a split manhole's load for one extra pattern step and integrates to the end time plus one sub-step.

## Not supported yet

Wet wells, pumps and force mains, junction chambers, open channels, storm (rational method) loads, curve-based flow splits, design mode and water quality. `ise2inp.py` and the engine stop with a message on any of these.

## Files

| Path | What |
|---|---|
| `src/solver/isewer.c` | input, network, steady state, HGL sweep, SWMM hooks |
| `src/solver/isewer_eps.c` | EPS routing and report-time outputs |
| `patches/swmm5-infosewer-on-epa-5.2.4.patch` | everything against EPA SWMM 5.2.4: 76 added and 3 changed lines in 13 EPA files, plus the 3 new files |
| `epa524/` | unmodified EPA SWMM 5.2.4 solver, for the null test |
| `tools/ise2inp.py`, `tools/iedb2ise.py` | converters |
| `tests/` | the three checks and their reference files |
| `docs/` | the two specifications and the Python references |

## Licence and provenance

EPA SWMM 5.2.4 is in the public domain. The swmm5-infosewer code was written from a reading of the InfoSewer engine source (copyright 1996-2011 MWH Soft and Innovyze, now Autodesk). It re-implements the engine's algorithms; no InfoSewer source is in this repository, and the specifications quote only short excerpts with line references. The repository is private until Autodesk agrees to an open release.

The test models include 125_21A and catchments of the Thornton 2009 model, both real utility networks. Remove `tests/` model data, or get the owners' agreement, before this repository goes public.

## AI use

swmm5-infosewer was built with Claude, Anthropic's AI model, in Claude Code sessions in November 2026, under the direction of Robert (Bob) Dickinson. He chose the project and its scope (gravity networks, steady state and EPS), supplied the InfoSewer engine source and the 125_21A and Thornton models from his files, and is responsible for the work.

Claude ported the InfoSewer engine to Linux to produce the reference results, wrote the two specifications and the Python references from the engine source, wrote the C in `src/solver/isewer*.c` and the EPA SWMM hooks, wrote the converters, the test variants, the check scripts and this README.

The checks are mechanical: every value is compared with the InfoSewer engine's own binary output, InfoSewer's stored IEDB results anchor the engine build, and the null test compares byte for byte with EPA SWMM 5.2.4. No person has yet reviewed the C line by line, run swmm5-infosewer on a model outside the test set, or checked the models with wet wells, pumps or force mains, which it does not run.
