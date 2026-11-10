#!/usr/bin/env python3
"""iedb2ise.py - write the InfoSewer engine (HeEngPkg) text input file from an
InfoSewer .IEDB database folder, the job InfoSewer's HeDataPkg did with
CHe_InputRun (Code/HeDataPkg/He_InputRun.Cpp).  Each section below is a
line-for-line port of the matching Fv_Make* function: same field numbers,
same order, same printf formats (%-lf = %f, %-ld = %d).

Usage:
  python3 -I iedb2ise.py <IEDB folder> <out.txt> [--mode steady|design|eps]
                         [--scenario ID] [--full] [--keep-report-path]
                         [--not-pro] [--epa-design-flow]

What the C++ got from its caller, and where this script gets it:
  report / time / option set IDs   SCENARIO.DBF row of the scenario
                                   (TIME_SET, REPORT_SET, OPTION_SET), else
                                   Project.ini [Scenario], else "BASE"
  data folder                      the IEDB root (the active data sets as
                                   InfoSewer last saved them)
  use-active flag                  ANODE.DBF non-empty and Project.ini
                                   [Facility] Full Network = FALSE (--full
                                   forces the whole network)
  IsPro / IsAppLicSuite            assumed true (InfoSewer Pro); --not-pro
  m_bUseEPALogic (design only)     --epa-design-flow, default 0

Relations: the C++ walks a CodeBase relation (master + related tables joined
on ID, in ID-index order).  Here every related table must hold the master's
ID (inner join), which is how the ANODE/ALINK tables filter the network.

Written by Claude (AI) for SWMM Enablement, 2026-11-09, from the HeDataPkg
source; IEDB contents are read as data only.
"""
import os, sys, glob, argparse, configparser

try:
    import dbfread
except ImportError:
    sys.exit("needs dbfread:  pip install dbfread")

# ---- enums (He_InputRun.h, He_Common.h) -----------------------------------
RUN_STATIC, RUN_DESIGN, RUN_DYNAMIC = 0, 1, 2
PAT_STEP, PAT_LINEAR, PAT_HYDRO = 0, 1, 4
CUR_WETWELL, CUR_CRITERIA, CUR_COST, CUR_SPLIT, CUR_PEAKFLOW, CUR_INTENSITY, \
    CUR_HYDGRAPH, CUR_PEAKPOP, CUR_MINSLOPE = range(9)
LOAD_COVERAGE, LOAD_CEILING = 2, 3
JCT_LOADING, JCT_CHAMBER, JCT_OUTLET = 0, 1, 2
PIPE_GRAVITY, PIPE_FORCE = 0, 1
CTRL_VOLUME, CTRL_LEVEL, CTRL_DISCHARGE, CTRL_INFLOW = 0, 1, 2, 3
SPLIT_HYDRAULIC, SPLIT_CEILING = 0, 4
Q_H2S, Q_CORROSION = 6, 7
DESIGN_NEW = 0
K_SECOND, K_HOUR, K_CLOCK, K_CEILING = 0, 2, 4, 7
TIME_UNITS = ["Seconds", "Minutes", "Hours", "Days", "", "AM", "PM"]
FLOW_UNITS = ["CFS", "GPM", "MGD", "IMGD", "AFD", "LPS", "LPM", "MLD", "CMH", "CMD"]


class WriterError(Exception):
    pass


# ---- CHnDBField semantics --------------------------------------------------
def blank(v):
    return v is None or (isinstance(v, str) and v.strip() == "")

def S(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "T" if v else "F"
    if isinstance(v, float):
        return ("%f" % v).rstrip("0").rstrip(".")
    return str(v).strip()

def D(v):
    if blank(v):
        return 0.0
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0

def L(v):
    if blank(v):
        return 0
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0

def B(v):
    if isinstance(v, bool):
        return v
    if blank(v):
        return False
    return str(v).strip().upper() in ("T", "Y", "1", "TRUE")

def f(x):   # %-lf
    return "%f" % x

def d(x):   # %-ld / %i
    return "%d" % x


# ---- DBF access ------------------------------------------------------------
class Table:
    """One DBF: records in physical order, field names in order, ID index."""
    def __init__(self, path):
        self.path = path
        self.fields, self.recs = [], []
        if path and os.path.exists(path):
            t = dbfread.DBF(path, encoding="latin-1", ignore_missing_memofile=True,
                            char_decode_errors="replace")
            self.fields = [x.name.upper() for x in t.fields]
            for r in t:
                self.recs.append([r[n] for n in t.field_names])
        self.ix = {}
        if self.fields and self.fields[0] == "ID":
            for r in self.recs:
                k = S(r[0]).upper()
                if k and k not in self.ix:
                    self.ix[k] = r

    def exists(self):
        return bool(self.fields)

    def has(self, name):
        return name.upper() in self.fields

    def col(self, name):
        """1-based field number like FindField; 0 if missing."""
        try:
            return self.fields.index(name.upper()) + 1
        except ValueError:
            return 0

    def get(self, rec, num):
        """FetchField(field, num): None when the field number does not exist."""
        if rec is None or num <= 0 or num > len(self.fields):
            return None
        return rec[num - 1]

    def byname(self, rec, name):
        n = self.col(name)
        return self.get(rec, n) if n else None

    def seek(self, key):
        return self.ix.get(S(key).upper())

    def sorted_recs(self):
        return sorted((r for r in self.recs if not blank(r[0])), key=lambda r: S(r[0]))


class Relate:
    """Master table plus related tables joined on ID (inner join)."""
    def __init__(self, tables):
        self.t = tables

    def rows(self):
        m = self.t[0]
        for r in m.sorted_recs():
            k = S(r[0]).upper()
            row = [r]
            ok = True
            for t in self.t[1:]:
                rr = t.ix.get(k)
                if rr is None:
                    ok = False
                    break
                row.append(rr)
            if ok:
                yield row

    def get(self, row, num, table=0):
        return self.t[table].get(row[table], num)

    def find(self, name):
        """relate.FindField(name): (table, number) of the first table that has it."""
        for i, t in enumerate(self.t):
            n = t.col(name)
            if n:
                return i, n
        return None, 0

    def byname(self, row, name):
        i, n = self.find(name)
        return None if i is None else self.t[i].get(row[i], n)


# ---- the writer ------------------------------------------------------------
class InputRun:
    TABLES = {
        "Junction": "MANHOLE", "JctHydData": "MHHYD", "Tank": "WWELL",
        "TankHydData": "WWELLHYD", "Pipe": "PIPE", "PipeHydData": "PIPEHYD",
        "Pump": "PUMP", "PumpHydData": "PUMPHYD", "Node": "NODE", "Link": "LINK",
        "ActiveNode": "ANODE", "ActiveLink": "ALINK", "Time": "TIME",
        "Report": "REPORT", "SimulationOption": "RUNOPT", "Status": "STATUS",
        "Control": "CONTROL", "Xload": "XLOAD", "PipeSplit": "PIPESPLT",
        "PipeDesign": "PIPEDGN", "PipeConstraint": "PIPECNST", "PipeII": "PIPEII",
        "Pattern": "PATTERN", "PatternDetail": "PATNDATA", "Curve": "CURVE",
        "CurveDetail": "CUVDATA", "ManCover": "MANCOVER", "DischargeHGL": "DISCGHGL",
        "NodeSource": "NDSOURCE", "MHDesign": "MHOFFSET", "Hydrograph": "HYDROGRH",
        "PipeQuality": "PIPEWQ", "PipeCrossSection": "NCPIXSEC",
        "PipeDesignNC": "NCPIDGN", "Scenario": "SCENARIO",
    }

    def __init__(self, folder, mode, scenario=None, full=False, keep_report_path=False,
                 pro=True, epa_design_flow=False):
        self.folder = folder
        self.mode = mode
        self.pro = pro
        self.suite = pro
        self.keep_report_path = keep_report_path
        self.epa_design_flow = epa_design_flow
        self.warnings = []
        self._files = {os.path.basename(p).upper(): p for p in glob.glob(os.path.join(folder, "*"))}
        self._cache = {}
        self.ini = configparser.ConfigParser(strict=False, interpolation=None)
        ini = self._files.get("PROJECT.INI")
        if ini:
            try:
                self.ini.read(ini, encoding="latin-1")
            except configparser.Error as e:
                self.warn("Project.ini not readable: %s" % e)
        self.scenario = scenario or self._ini("Scenario", "Active Scenario ID") or "BASE"
        sc = self.T("Scenario").seek(self.scenario)
        def pick(col, inikey):
            v = self.T("Scenario").byname(sc, col) if sc is not None else None
            if not blank(v):
                return S(v)
            v = self._ini("Scenario", inikey)
            return v if v else "BASE"
        self.timeSet = pick("TIME_SET", "Time Set")
        self.reportSet = pick("REPORT_SET", "Report Set")
        self.optionSet = pick("OPTION_SET", "Simulation Option Set")
        fullnet = (self._ini("Facility", "Full Network") or "FALSE").upper() == "TRUE"
        self.useActive = not (full or fullnet) and len(self.T("ActiveNode").recs) >= 1
        self.qualityType = 0
        self.out = []

    def _ini(self, sec, key):
        try:
            return self.ini.get(sec, key).strip()
        except (configparser.Error, KeyError):
            return ""

    def warn(self, m):
        self.warnings.append(m)

    def T(self, name):
        if name not in self._cache:
            fn = self.TABLES[name] + ".DBF"
            p = self._files.get(fn)
            self._cache[name] = Table(p)
        return self._cache[name]

    def w(self, s):
        self.out.append(s)

    def opt_id(self, v):          # Fv_SetOptionalIdField
        s = S(v)
        self.w("\t-" if s == "" else "\t%s" % s.upper())

    def opt_data(self, v, string=False):   # Fv_SetOptionalDataField
        if blank(v):
            self.w("\t-")
        elif string:
            self.w("\t%s" % S(v))
        else:
            self.w("\t" + f(D(v)))

    def need(self, rel, *names):
        for n in names:
            if rel.find(n)[0] is None:
                raise WriterError("%s: field %s missing" % (rel.t[0].path, n))

    # -- setup (Fv_SetupPatternType / Fv_SetupCurveType) --
    def setup(self):
        pat, opt = self.T("Pattern"), self.T("SimulationOption")
        optrec = opt.seek(self.optionSet)
        self.patList, self.patType = [], []
        if optrec is None:
            self.warn("option set %s not found in RUNOPT" % self.optionSet)
        for r in pat.recs:
            self.patList.append(S(r[0]).upper())
            self.patType.append(PAT_STEP)
        if optrec is not None and L(opt.byname(optrec, "LD_PATN_TP")) != PAT_STEP:
            self.patType = [PAT_LINEAR if t == PAT_STEP else t for t in self.patType]
        n = len(pat.recs)
        hyd = self.T("Hydrograph")
        for r in hyd.recs:
            pid = S(hyd.byname(r, "PATTERN")).upper()
            if pid not in self.patList and len(self.patList) < n:
                self.patList.append(pid)
                self.patType.append(PAT_HYDRO)

        cur = self.T("Curve")
        self.curList = [S(r[0]).upper() for r in cur.recs]
        self.curType = [-1] * len(self.curList)

        def mark(cid, ctype):
            cid = S(cid).upper()
            if not cid:
                return
            if cid in self.curList:
                i = self.curList.index(cid)
                if self.curType[i] >= 0 and self.curType[i] != ctype:
                    raise WriterError("Curve %s has already been assigned with a different application type" % cid)
                self.curType[i] = ctype

        ww = self.T("TankHydData")
        for r in ww.recs:
            mark(ww.byname(r, "CURVE"), CUR_WETWELL)
        dg = self.T("PipeDesign")
        for r in dg.recs:              # C++ 'continue' skips CURVE_ANA when CURVE_DGN is blank
            if blank(dg.byname(r, "CURVE_DGN")):
                continue
            mark(dg.byname(r, "CURVE_DGN"), CUR_CRITERIA)
            mark(dg.byname(r, "CURVE_ANA"), CUR_CRITERIA)
        for r in dg.recs:
            if blank(dg.byname(r, "COST_REPL")):
                continue
            mark(dg.byname(r, "COST_REPL"), CUR_COST)
            mark(dg.byname(r, "COST_DUP"), CUR_COST)
        if optrec is not None:
            o = lambda n: opt.byname(optrec, n)
            mark(o("DGN_CURVE"), CUR_CRITERIA)
            mark(o("ANA_CURVE"), CUR_CRITERIA)
            mark(o("COST_REPL"), CUR_COST)
            mark(o("COST_DUP"), CUR_COST)
            mark(o("SLP_CURVE"), CUR_MINSLOPE)
            mark(o("PEAK_CURVE"), CUR_PEAKFLOW)
            mark(o("PF_CURVE"), CUR_PEAKPOP)
            mark(o("IDF_CURVE"), CUR_INTENSITY)
            mark(o("HYETOGRAPH"), CUR_HYDGRAPH)
        sp = self.T("PipeSplit")
        for r in sp.recs:
            mark(sp.byname(r, "CURVE"), CUR_SPLIT)
        mh = self.T("JctHydData")
        for r in mh.recs:
            mark(mh.byname(r, "HYETOGRAPH"), CUR_HYDGRAPH)
        for r in dg.recs:
            mark(dg.byname(r, "SLP_CURVE"), CUR_MINSLOPE)

    # -- sections --
    def title(self):
        t = self.T("Report")
        r = t.seek(self.reportSet)
        self.w("[Titles]\n")
        if r is None:
            return
        for n in ("TITLE1", "TITLE2", "TITLE3"):
            v = t.byname(r, n)
            if t.has(n) and not blank(v):
                self.w("%s\n" % S(v))
        self.w("\n")

    def options(self):
        t = self.T("SimulationOption")
        r = t.seek(self.optionSet)
        if r is None:
            return
        def g(n):
            return t.byname(r, n) if t.has(n) else None
        def present(n):
            return t.has(n) and not blank(g(n))
        w = self.w
        w("[Options]\n")
        if present("LOAD_PATN"):
            w("Default Manhole Diurnal Pattern ID : %s\n" % S(g("LOAD_PATN")).upper())
        if present("DGN_CURVE"):
            w("Default Design Criteria Curve ID : %s\n" % S(g("DGN_CURVE")).upper())
        if present("ANA_CURVE"):
            w("Default Analysis Criteria Curve ID : %s\n" % S(g("ANA_CURVE")).upper())
        if present("COST_REPL"):
            w("Default Pipe Replacement Cost ID : %s\n" % S(g("COST_REPL")).upper())
        if present("COST_DUP"):
            w("Default Pipe Duplication Cost ID : %s\n" % S(g("COST_DUP")).upper())
        if present("LOAD_SCALE"):
            w("Global Loading Multiplier : %s\n" % f(D(g("LOAD_SCALE"))))
        if present("FLOW_UNIT"):
            u = L(g("FLOW_UNIT"))
            if 0 <= u < len(FLOW_UNITS):
                w("Flow Unit : %s\n" % FLOW_UNITS[u])
        if present("FM_NETWORK"):
            w("Advanced Forcemain Network Support : %s\n" % ("YES" if B(g("FM_NETWORK")) else "NO"))
        for n, lab in (("PF_K", "Peaking Factor k"), ("PF_P", "Peaking Factor p"),
                       ("PF_A", "Peaking Parameter a"), ("PF_B", "Peaking Parameter b"),
                       ("PF_C", "Peaking Parameter c"), ("PF_D", "Peaking Parameter d")):
            if present(n):
                w("%s : %s\n" % (lab, f(D(g(n)))))
        if present("PEAK_CURVE"):
            w("Peaking Factor Curve : %s\n" % S(g("PEAK_CURVE")))
        if t.has("PF_E"):
            w("Peaking Parameter e : %s\n" % f(1.0 if blank(g("PF_E")) else D(g("PF_E"))))
        if present("PF_CURVE"):
            w("Coverage Based Peaking Curve : %s\n" % S(g("PF_CURVE")))
        designMode = -1
        if present("DGN_MODE"):
            designMode = L(g("DGN_MODE"))
            w("Design Mode : %d\n" % designMode)
        if self.mode == RUN_DESIGN:
            w("Design Flow Type : %d\n" % (1 if self.epa_design_flow else 0))
        for n, lab in (("MIN_VEL", "Minimum Design Velocity"), ("MAX_VEL", "Maximum Design Velocity"),
                       ("MANNING", "Design Manning Coefficient")):
            if present(n):
                w("%s : %s\n" % (lab, f(D(g(n)))))
        if present("DGNEXIST"):
            w("Duplication Design Option : %d\n" % (1 if B(g("DGNEXIST")) else 0))
        if present("MAX_DROP"):
            w("Maximum Crown Drop : %s\n" % f(D(g("MAX_DROP"))))
        if present("MIN_SLOPE"):
            w("Minimum Design Slope : %s\n" % f(D(g("MIN_SLOPE"))))
        if present("MAX_SLOPE"):
            w("Maximum Design Slope : %s\n" % f(D(g("MAX_SLOPE"))))
        if present("SLP_CURVE"):
            w("Default Minimum Slope Curve : %s\n" % S(g("SLP_CURVE")))
        if present("MIN_COVER"):
            w("Minimum Design Cover : %s\n" % f(D(g("MIN_COVER"))))
        if present("MAX_COVER"):
            w("Maximum Design Cover : %s\n" % f(D(g("MAX_COVER"))))
        if designMode == DESIGN_NEW:
            if present("ADD_ERCOST") and B(g("ADD_ERCOST")):
                w("Enable Cost Equations : 1\n")
                for n, lab in (("EXCA_ANGLE", "Default Pipe Excavation Angle"),
                               ("CLR_CTF", "Clearance Cutoff Diameter"), ("CLR_SMALL", "Small Pipe Clearance"),
                               ("CLR_LARGE", "Big Pipe Clearance"), ("UVCF_CTF", "Cost Cutoff Depth"),
                               ("UVCF_SHLW", "Shallow Pipe Cost"), ("UVCF_DEEP", "Deep Pipe Cost")):
                    if present(n):
                        w("%s : %s\n" % (lab, f(D(g(n)))))
            else:
                w("Enable Cost Equations : 0\n")
            if present("FIXED_IVT") and B(g("FIXED_IVT")):
                w("Fixed Invert Elevation : 1\n")
            if present("DROP_CTF"):
                w("Drop Manhole Limiting Diameter : %s\n" % f(D(g("DROP_CTF"))))
        if present("LD_PATN_TP"):
            w("Loading Pattern Type : %d\n" % L(g("LD_PATN_TP")))
        if present("MIN_TRAVEL"):
            w("Minimum Pipe Segment : %s\n" % f(D(g("MIN_TRAVEL"))))
        if present("ACCURACY"):
            w("Flow Tolerance : %s\n" % f(D(g("ACCURACY"))))
        if present("LOCKMCVR"):
            w("Manhole Head Type : %d\n" % (1 if B(g("LOCKMCVR")) else 0))
        if present("FLOWCONT"):
            w("Pump Flow Continuity : %d\n" % (1 if B(g("FLOWCONT")) else 0))
        if present("FLOW_ATN"):
            w("Flow Attenuation : %d\n" % (1 if B(g("FLOW_ATN")) else 0))
        if present("MAX_SEG"):
            w("Maximum Number of Segments : %d\n" % L(g("MAX_SEG")))
        if present("MIN_TLSP"):
            w("Minimum Travel Time : %d\n" % L(g("MIN_TLSP")))
        for n, lab in (("INV_OFFSET", "Invert Level Offset"), ("BODRRC", "Reaction Rate Constant"),
                       ("BODTEMP", "Temperature"), ("PTSETVEL", "Settling Velocity"),
                       ("LMTFLWVEL", "Scouring Velocity")):
            if present(n):
                w("%s : %s\n" % (lab, f(D(g(n)))))
        if present("QUALITY"):
            q = L(g("QUALITY"))
            if q in (Q_H2S, Q_CORROSION) and not self.suite:
                q = 0
            self.qualityType = q
            w("Quality Type : %d\n" % q)
            if q in (Q_H2S, Q_CORROSION):
                for n, lab in (("H2S_FLUX", "Sulfide Flux Coefficient"), ("H2S_LOSS", "Sulfide Loss Coefficient"),
                               ("H2S_PH", "pH"), ("H2S_ION", "Ionization Constant"),
                               ("H2S_SLB", "Percent of Soluble Sulfide")):
                    if present(n):
                        w("%s : %s\n" % (lab, f(D(g(n)))))
            if q == Q_CORROSION:
                for n, lab in (("PIPE_MAT", "Default Pipe Material"), ("ALKALINITY", "Default Pipe Alkalinity"),
                               ("ACID_REACT", "Acid Reaction Efficiency"), ("H2S_SAT", "Relative Saturation of H2S")):
                    if present(n):
                        w("%s : %s\n" % (lab, f(D(g(n)))))
        if present("TRACE_ND"):
            w("Quality Source Node : %s\n" % S(g("TRACE_ND")))
        if self.pro:
            if present("RUN_STORM"):
                w("Storm Runoff Analysis : %s\n" % ("1" if B(g("RUN_STORM")) else "0"))
            if present("MIN_TOC"):
                w("Minimum Concentration Time : %s\n" % f(D(g("MIN_TOC"))))
            if present("IDF_CURVE"):
                w("Intensity Duration Curve : %s\n" % S(g("IDF_CURVE")))
            if present("HYETOGRAPH"):
                w("Default Rainfall Curve : %s\n" % S(g("HYETOGRAPH")))
            if present("IDF_IUNIT"):
                w("Intensity Unit : %d\n" % L(g("IDF_IUNIT")))
            if present("IDF_DUNIT"):
                w("Duration Unit : %d\n" % L(g("IDF_DUNIT")))
            if present("HYDROGRAPH"):
                w("Default Unit Hydrograph : %s\n" % S(g("HYDROGRAPH")))
            if present("HYDRO_DUR"):
                w("Unit Hydrograph Duration : %s\n" % f(D(g("HYDRO_DUR"))))
            if present("UH_VOLADJ"):
                v = D(g("UH_VOLADJ"))
                w("Unit Hydrograph Volume Threshold : %s\n" % ("0.0" if v < 0.0 else f(v)))
            else:
                w("Unit Hydrograph Volume Threshold : 5.0\n")
            if present("SD_METHOD") and B(g("SD_METHOD")):
                w("San Diego Rational : 1\n")
            if present("PREC_6HR") and D(g("PREC_6HR")) > 0.0:
                w("Six Hour Rainfall : %s\n" % f(D(g("PREC_6HR"))))
        w("\n")

    def _grouped(self, table, sortcol):
        groups, order = {}, []
        for r in table.sorted_recs():
            k = S(r[0]).upper()
            if k not in groups:
                groups[k] = (S(r[0]), [])
                order.append(k)
            groups[k][1].append(r)
        for k in order:
            rid, recs = groups[k]
            yield k, rid, sorted(recs, key=lambda r: D(r[sortcol]))

    def patterns(self):
        t = self.T("PatternDetail")
        self.w("[Patterns]\n")
        for k, rid, recs in self._grouped(t, 1):     # OrderBy("SEQ"), value At(2)
            if k not in self.patList:
                raise WriterError("pattern %s in PATNDATA but not in PATTERN" % rid)
            ptype = self.patType[self.patList.index(k)]
            if ptype < 0:
                continue
            n = len(recs)
            for i, r in enumerate(recs):
                if i % 6 == 0:
                    self.w("%s %d" % (rid, ptype))
                self.w("\t" + f(D(r[2])))
                if i == n - 1 or i % 6 == 5:
                    self.w("\n")
        self.w("\n")

    def curves(self):
        t = self.T("CurveDetail")
        self.w("[Curves]\n")
        for k, rid, recs in self._grouped(t, 1):     # OrderBy("X")
            if k not in self.curList:
                raise WriterError("curve %s in CUVDATA but not in CURVE" % rid)
            ctype = self.curType[self.curList.index(k)]
            if ctype < 0:
                continue
            for i, r in enumerate(recs):
                if i == 0:
                    self.w("%s %d\t%s\t%s\n" % (rid, ctype, f(D(r[1])), f(D(r[2]))))
                else:
                    self.w("\t\t%s\t%s\n" % (f(D(r[1])), f(D(r[2]))))
        self.w("\n")

    def design_tables(self):
        if not self.pro:
            return
        t = self.T("PipeDesignNC")
        cols = {}
        for n in ("SHAPE", "GEOM1", "GEOM2", "GEOM3", "GEOM4", "ANA_DD", "DGN_DD", "RPL_COST", "DUP_COST"):
            c = t.col(n) - 1
            if c <= 0:
                raise WriterError("NCPIDGN: field %s missing" % n)
            cols[n] = c
        self.w("[Design Tables]\n")
        recs = sorted((r for r in t.recs if not blank(r[0])), key=lambda r: L(r[cols["SHAPE"]]))
        def flush(oldID, rs):
            for i, r in enumerate(rs):
                if i == 0:
                    self.w("%s \t%d" % (oldID, L(r[cols["SHAPE"]])))
                else:
                    self.w("\t")
                for n in ("GEOM1", "GEOM2", "GEOM3", "GEOM4", "ANA_DD", "DGN_DD", "RPL_COST", "DUP_COST"):
                    v = r[cols[n]]
                    self.w("\t" + f(0.0 if blank(v) else D(v)))
                self.w("\n")
        oldID, oldShape, rs = "", -1, []
        for r in recs:
            shape = L(r[cols["SHAPE"]])
            if oldShape != shape:
                flush(oldID, rs)
                rs = []
            oldID, oldShape = S(r[0]), shape
            rs.append(r)
        if rs:
            flush(oldID, rs)
        self.w("\n")

    def hydrograph(self):
        if not self.pro:
            return
        t = self.T("Hydrograph")
        names = ["TYPE", "AREA", "PATTERN", "RFVOL_EFF", "RFVOL1", "RFVOL2", "PEAKTM1", "PEAKTM2",
                 "PEAKTM3", "RECESS1", "RECESS2", "RECESS3", "NRCS_LAGTM", "NRCS_HYDLN", "NRCS_CUVNO",
                 "NRCS_SLOPE", "WS_LEN", "WS_CENTRO", "WS_SLOPE", "WS_PERIMP"]
        for n in names:
            if not t.has(n):
                raise WriterError("HYDROGRH: field %s missing" % n)
        self.w("[Unit Hydrograph]\n")
        defaults = [-1, 0.0, "-", 30.0, 20.0, 40.0, 1.0, 4.0, 12.0, 2.0, 2.0, 2.0,
                    0.0, 0.0, 75.0, 0.0, 5000.0, 2500.0, 0.01, 0.75]
        for r in t.recs:
            if blank(r[0]):
                continue
            vals = []
            for n, dv in zip(names, defaults):
                v = t.byname(r, n)
                if blank(v):
                    vals.append(dv)
                elif n == "TYPE":
                    vals.append(L(v))
                elif n == "PATTERN":
                    vals.append(S(v))
                else:
                    vals.append(D(v))
            line = S(r[0]) + "\t%d" % vals[0] + "\t" + f(vals[1]) + "\t%s" % vals[2]
            line += "".join("\t" + f(x) for x in vals[3:])
            self.w(line + "\n")
        self.w("\n")

    def load_nodes(self):
        files = [self.T("JctHydData"), self.T("Junction")]
        if self.useActive:
            files.append(self.T("ActiveNode"))
        rel = Relate(files)
        uh = ["SWS_AREA", "SWS_RUNOFF", "SWS_SLOPE", "SWS_LEN", "SWS_TOC", "HYETOGRAPH", "HYDROGRAPH",
              "SWS_PERIMP", "SWS_PERV", "SWS_IMPERV", "SWS_INFIL1", "SWS_INFIL2", "SWS_DECAY", "SWS_REGEN"]
        self.need(rel, *uh)
        mh = self.T("JctHydData")
        self.w("[Loading Manholes]\n")
        for row in rel.rows():
            m, e = row[0], row[1]
            if blank(rel.get(row, 1, 1)) or blank(rel.get(row, 1, 0)):
                continue
            if L(rel.get(row, 3, 1)) != JCT_LOADING:
                continue
            s = S(rel.get(row, 1, 1))
            s += "\t" + f(D(mh.get(m, 3))) + "\t" + f(D(mh.get(m, 2))) + "\t" + f(D(mh.get(m, 4)))
            for n in ("SWS_AREA", "SWS_RUNOFF", "SWS_SLOPE", "SWS_LEN", "SWS_TOC"):
                s += "\t" + f(D(mh.byname(m, n)))
            for n in ("HYETOGRAPH", "HYDROGRAPH"):
                v = S(mh.byname(m, n))
                s += "\t-" if v == "" else "\t" + v.upper()
            for n, dv in (("SWS_PERIMP", 0.75), ("SWS_PERV", 0.0), ("SWS_IMPERV", 0.0),
                          ("SWS_INFIL1", 0.0), ("SWS_INFIL2", 0.0), ("SWS_DECAY", 0.0), ("SWS_REGEN", 0.0)):
                v = mh.byname(m, n)
                s += "\t" + f(dv if blank(v) else D(v))
            for i in range(10):
                ix = i * 4 + 5
                vals = [mh.get(m, ix + j) for j in range(4)]
                if ix + 3 > len(mh.fields):
                    continue
                ltype = L(vals[1])
                if ltype < 0 or ltype >= LOAD_CEILING:
                    ltype = 0
                cov = 0.0 if (ltype == LOAD_COVERAGE and blank(vals[3])) else D(vals[3])
                p = S(vals[2])
                s += "\t" + f(D(vals[0])) + "\t%d" % ltype + ("\t-" if p == "" else "\t" + p.upper()) + "\t" + f(cov)
            self.w(s + "\n")
        self.w("\n")

    def _node_rel(self):
        files = []
        if self.useActive:
            files.append(self.T("ActiveNode"))
        files += [self.T("Junction"), self.T("JctHydData")]
        return Relate(files), len(files) - 2, len(files) - 1

    def chambers(self):
        rel, et, mt = self._node_rel()
        self.w("[Junction Chambers]\n")
        for row in rel.rows():
            if blank(rel.get(row, 1, et)) or blank(rel.get(row, 1, mt)):
                continue
            if L(rel.get(row, 3, et)) != JCT_CHAMBER:
                continue
            self.w(S(rel.get(row, 1, et)) + "\t" + f(D(rel.get(row, 3, mt))) + "\t" +
                   f(D(rel.get(row, 2, mt))) + "\t" + f(D(rel.get(row, 4, mt))) + "\n")
        self.w("\n")

    def outlets(self):
        rel, et, mt = self._node_rel()
        self.w("[Outlets]\n")
        for row in rel.rows():
            if blank(rel.get(row, 1, et)) or blank(rel.get(row, 1, mt)):
                continue
            if L(rel.get(row, 3, et)) != JCT_OUTLET:
                continue
            self.w(S(rel.get(row, 1, et)) + "\t" + f(D(rel.get(row, 3, mt))) + "\t" +
                   f(D(rel.get(row, 2, mt))) + "\n")
        self.w("\n")

    def wet_wells(self):
        files = [self.T("TankHydData"), self.T("Node")]
        if self.useActive:
            files.append(self.T("ActiveNode"))
        rel = Relate(files)
        if not files[0].exists():
            raise WriterError("WWELLHYD.DBF missing")
        self.need(rel, "TYPE", "BTM_ELEV", "MIN_LEVEL", "MAX_LEVEL", "INIT_LEVEL", "DIAMETER", "CURVE", "HEADLOSS")
        t = files[0]
        self.w("[Wet Wells]\n")
        for row in rel.rows():
            r = row[0]
            g = lambda n: t.byname(r, n)
            c = S(g("CURVE"))
            self.w(S(r[0]) + "\t%d" % L(g("TYPE")) + "".join("\t" + f(D(g(n))) for n in
                   ("BTM_ELEV", "MIN_LEVEL", "MAX_LEVEL", "INIT_LEVEL", "DIAMETER")) +
                   ("\t" + c.upper() if c else "\t-") + "\t" + f(D(g("HEADLOSS"))) + "\n")
        self.w("\n")

    def _pipes(self, header, ptype):
        files = [self.T("PipeHydData"), self.T("Pipe"), self.T("Link")]
        if self.useActive:
            files.append(self.T("ActiveLink"))
        rel = Relate(files)
        self.w(header + "\n")
        for row in rel.rows():
            if L(rel.get(row, 3, 1)) != ptype:
                continue
            g = lambda n, tb=0: rel.get(row, n, tb)
            self.w(S(g(1)) + "\t" + S(g(2, 2)) + "\t" + S(g(3, 2)) + "".join("\t" + f(D(g(n))) for n in
                   (2, 3, 4, 5, 6)) + "\t%d" % L(g(7)) + "\n")
        self.w("\n")

    def pipes(self):
        self._pipes("[Pipes]", PIPE_GRAVITY)

    def channel_shape(self):
        if not self.pro:
            return
        t = self.T("PipeCrossSection")
        for n in ("SHAPE", "GEOM1", "GEOM2", "GEOM3", "GEOM4"):
            if not t.has(n):
                raise WriterError("NCPIXSEC: field %s missing" % n)
        self.w("[Channel Shape]\n")
        for r in t.recs:
            if blank(r[0]):
                continue
            g = lambda n: t.byname(r, n)
            shape = L(g("SHAPE"))      # C++ leaves 'shape' uninitialised when blank
            if blank(g("SHAPE")):
                self.warn("channel shape %s has no SHAPE; written as 0" % S(r[0]))
            self.w(S(r[0]) + "\t%d" % shape + "".join("\t" + f(0.0 if blank(g(n)) else D(g(n)))
                   for n in ("GEOM1", "GEOM2", "GEOM3", "GEOM4")) + "\n")
        self.w("\n")

    def force_mains(self):
        self._pipes("[Force Mains]", PIPE_FORCE)

    def pumps(self):
        files = [self.T("PumpHydData"), self.T("Link")]
        if self.useActive:
            files.append(self.T("ActiveLink"))
        rel = Relate(files)
        st, ct = self.T("Status"), self.T("Control")
        self.w("[Pumps]\n")
        for row in rel.rows():
            g = lambda n, tb=0: rel.get(row, n, tb)
            pid = S(g(1))
            speed = 1.0
            sr = st.seek(pid)
            if sr is not None and st.has("SETTING") and not blank(st.byname(sr, "SETTING")):
                speed = D(st.byname(sr, "SETTING"))
            desired, base, patID = 0.0, 0.0, ""
            cr = ct.seek(pid)
            if cr is not None:
                if not blank(ct.byname(cr, "STA_FLOW")):
                    desired = D(ct.byname(cr, "STA_FLOW"))
                if not blank(ct.byname(cr, "BASE_FLOW")):
                    base = D(ct.byname(cr, "BASE_FLOW"))
                if not blank(ct.byname(cr, "PATTERN")):
                    patID = S(ct.byname(cr, "PATTERN"))
            outPat = "-"
            if pid and patID and cr is not None:
                if L(ct.byname(cr, "METHOD")) not in (CTRL_VOLUME, CTRL_LEVEL, CTRL_INFLOW):
                    outPat = patID
            self.w(pid + "\t" + S(g(2, 1)) + "\t" + S(g(3, 1)) + "\t%d\t%d" % (L(g(2)), L(g(3))) +
                   "".join("\t" + f(D(g(n))) for n in (4, 5, 6, 7, 8, 9)) + "\t" + f(speed) +
                   "\t" + outPat + "\t" + f(desired) + "\t" + f(base) + "\n")
        self.w("\n")

    def _link_rel(self, first, second):
        files = [self.T(first), self.T(second)]
        if self.useActive:
            files.append(self.T("ActiveLink"))
        return Relate(files)

    def _node_extra(self, first):
        files = [self.T(first), self.T("Junction")]
        if self.useActive:
            files.append(self.T("ActiveNode"))
        return Relate(files)

    def status(self):
        rel = self._link_rel("Status", "Pump")
        self.w("[Status]\n")
        for row in rel.rows():
            db = L(rel.get(row, 2))
            kernel = 0 if db == 0 else 1
            self.w(S(rel.get(row, 1)) + "\t%d\n" % kernel)
        self.w("\n")

    def controls(self):
        rel = self._link_rel("Control", "Pump")
        self.need(rel, "METHOD", "ON", "OFF", "MIN_SPEED", "MAX_SPEED")
        self.w("[Controls]\n")
        for row in rel.rows():
            b = lambda n: rel.byname(row, n)
            typ = L(b("METHOD"))
            s = S(rel.get(row, 1)) + "\t%d" % typ + "\t" + f(D(b("ON"))) + "\t" + f(D(b("OFF")))
            if typ in (CTRL_DISCHARGE, CTRL_INFLOW):
                for n in ("MIN_SPEED", "MAX_SPEED"):
                    v = b(n)
                    s += "\t-" if blank(v) else "\t" + f(D(v))
            self.w(s + "\n")
        self.w("\n")

    def extra_loads(self):
        rel = self._node_extra("Xload")
        self.w("[Extra Loads]\n")
        for row in rel.rows():
            g = lambda n: rel.get(row, n)
            if B(g(2)):
                continue
            lt = L(g(4))
            if lt < 0 or lt >= LOAD_CEILING:
                continue
            if lt == LOAD_COVERAGE and blank(g(5)):
                continue
            p = S(g(5))
            self.w(S(g(1)) + "\t" + f(D(g(3))) + "\t%d" % lt + ("\t-" if p == "" else "\t" + p.upper()) +
                   "\t" + f(D(g(6))) + "\n")
        self.w("\n")

    def load_node_head(self):
        rel = self._node_extra("ManCover")
        self.w("[Load Node Head Type]\n")
        for row in rel.rows():
            self.w("%s %d\n" % (S(rel.get(row, 1)), 1 if B(rel.get(row, 2)) else 0))
        self.w("\n")

    def outlet_heads(self):
        rel = self._node_extra("DischargeHGL")
        self.w("[Outlet Heads]\n")
        for row in rel.rows():
            self.w(S(rel.get(row, 1)) + "\t" + f(D(rel.get(row, 2))) + "\n")
        self.w("\n")

    def quality_source(self):
        files = [self.T("NodeSource")]
        if self.useActive:
            files.append(self.T("ActiveNode"))
        rel = Relate(files)
        t = files[0]
        qcol = t.col("QUALITY") or t.col("SETTING")
        if not qcol:
            raise WriterError("NDSOURCE: QUALITY/SETTING missing")
        if not t.col("UBOD"):
            raise WriterError("NDSOURCE: UBOD missing")
        self.w("[Initial Concentration]\n")
        for row in rel.rows():
            s = S(rel.get(row, 1)) + "\t" + f(D(rel.get(row, qcol)))
            if self.qualityType in (Q_H2S, Q_CORROSION):
                s += "\t" + f(D(t.byname(row[0], "UBOD")))
            self.w(s + "\n")
        self.w("\n")

    def pipe_corrosion(self):
        if self.qualityType != Q_CORROSION:
            return
        files = [self.T("PipeQuality")]
        if self.useActive:
            files.append(self.T("ActiveLink"))
        rel = Relate(files)
        self.need(rel, "PIPE_MAT", "ALKALINITY")
        self.w("[Pipe Corrosion]\n")
        for row in rel.rows():
            self.w(S(rel.get(row, 1)) + "\t" + f(D(rel.byname(row, "PIPE_MAT"))) + "\t" +
                   f(D(rel.byname(row, "ALKALINITY"))) + "\n")
        self.w("\n")

    def flow_splits(self):
        rel = self._link_rel("PipeSplit", "Pipe")
        self.w("[Flow Splits]\n")
        for row in rel.rows():
            g = lambda n: rel.get(row, n)
            if blank(g(2)) and blank(g(3)):
                continue
            typ = L(g(2))
            if typ <= SPLIT_HYDRAULIC or typ >= SPLIT_CEILING:
                continue
            self.w(S(g(1)) + "\t%d" % typ)
            self.opt_data(g(3))
            self.opt_data(g(4), True)
            self.w("\n")
        self.w("\n")

    def infiltration(self):
        rel = self._link_rel("PipeII", "Pipe")
        self.w("[Infiltrate Rates]\n")
        for row in rel.rows():
            g = lambda n: rel.get(row, n)
            if blank(g(2)) or blank(g(3)):
                continue
            self.w(S(g(1)) + "\t%d" % L(g(2)) + "\t" + f(D(g(3))))
            self.opt_data(g(4))
            self.opt_id(g(5))
            self.w("\n")
        self.w("\n")

    def design_curves(self):
        rel = self._link_rel("PipeDesign", "Pipe")
        names = ["CURVE_DGN", "CURVE_ANA", "COST_REPL", "COST_DUP", "ML_DROP", "SLP_CURVE", "MAX_DIA", "EXCA_ANGLE"]
        t = rel.t[0]
        for n in names:
            if not t.has(n):
                raise WriterError("PIPEDGN: field %s missing" % n)
        self.w("[Design Curves]\n")
        for row in rel.rows():
            r = row[0]
            self.w(S(r[0]))
            for n in names:
                if n in ("ML_DROP", "MAX_DIA", "EXCA_ANGLE"):
                    self.opt_data(t.byname(r, n))
                else:
                    self.opt_id(t.byname(r, n))
            self.w("\n")
        self.w("\n")

    def design_offsets(self):
        rel = self._node_extra("MHDesign")
        self.need(rel, "OFFSET")
        self.w("[Design Offsets]\n")
        for row in rel.rows():
            self.w(S(rel.get(row, 1)) + "\t" + f(D(rel.byname(row, "OFFSET"))) + "\n")
        self.w("\n")

    def design_constraints(self):
        rel = self._link_rel("PipeConstraint", "Pipe")
        names = ["V_MIN", "V_MAX", "MIN_SLOPE", "MAX_SLOPE", "MIN_COVER", "MAX_COVER"]
        self.need(rel, *names)
        t = rel.t[0]
        self.w("[Design Constraints]\n")
        for row in rel.rows():
            self.w(S(row[0][0]))
            for n in names:
                self.opt_data(t.byname(row[0], n))
            self.w("\n")
        self.w("\n")

    @staticmethod
    def check_time(s, unit):
        """Fv_CheckTimeData: (ok, decimal time)."""
        if unit < K_SECOND or unit >= K_CEILING:
            return False, 0.0
        parts = [p for p in s.split(":") if p != ""][:4]
        vals = []
        for p in parts:
            try:
                vals.append(float(p.strip()))
            except ValueError:
                return False, 0.0
        if not vals:
            return False, 0.0
        while len(vals) < 3:
            vals.append(0.0)
        if unit < K_CLOCK:
            if len(parts) != 1 or vals[0] < 0.0:
                return False, 0.0
            return True, vals[0]
        if vals[0] < 0.0 or not (0.0 <= vals[1] <= 60.0) or not (0.0 <= vals[2] <= 60.0):
            return False, 0.0
        if vals[0] >= 13.0:    # C++ tests 'unitFlag = kClockAM' (assignment): always the AM check
            return False, 0.0
        return True, vals[0] + vals[1] / 60.0 + vals[2] / 3600.0

    def time(self):
        t = self.T("Time")
        r = t.seek(self.timeSet)
        if r is None:
            return
        dyn = self.mode == RUN_DYNAMIC
        self.w("[Time]\nInput Time Data\n")
        def entry(col, ucol, label, dflt, needpos, only_if_present=False):
            if not t.has(col):
                return
            v = t.byname(r, col)
            if blank(v):
                if dyn and dflt and not only_if_present:
                    self.w("%s : %s\n" % (label, dflt))
                return
            s = S(v)
            u = t.byname(r, ucol)
            unit = K_HOUR if blank(u) else L(u)
            ok, dec = self.check_time(s, unit)
            if not ok:
                raise WriterError("bad time %r (unit %d) in TIME.%s" % (s, unit, col))
            if needpos and dec <= 0.0:
                return
            self.w("%s : %s %s\n" % (label, s, TIME_UNITS[unit]))
        entry("RUNST", "RUNST_UNIT", "Simulation Start Time", "0 Hours", False)
        entry("DUR", "DUR_UNIT", "Simulation Duration", "24 Hours", False)
        entry("PNSP", "PNSP_UNIT", "Pattern Time Step", "1 Hours", True)
        entry("UPNSP", "UPNSP_UNIT", "Pump Pattern Time Step", "1 Hours", True)
        entry("RPSP", "RPSP_UNIT", "Report Time Step", "1 Hours", True)
        entry("PTSP", "PTSP_UNIT", "Simulation Time Step for Pipe", None, True, True)
        entry("UTSP", "UTSP_UNIT", "Simulation Time Step for Pump", None, True, True)
        self.w("\n")

    def report(self):
        t = self.T("Report")
        r = t.seek(self.reportSet)
        if r is None:
            return
        self.w("[Report]\n")
        g = lambda n: t.byname(r, n) if t.has(n) else None
        if not blank(g("FOLDER")):
            if self.keep_report_path:
                self.w("Output Report Folder : %s\n" % S(g("FOLDER")))
            else:
                self.warn("REPORT.FOLDER %r left out (Windows path); --keep-report-path keeps it" % S(g("FOLDER")))
        if not blank(g("FILENAME")):
            if self.keep_report_path:
                self.w("Output Report File : %s\n" % S(g("FILENAME")))
            else:
                self.warn("REPORT.FILENAME %r left out; --keep-report-path keeps it" % S(g("FILENAME")))
        if not blank(g("NODE")):
            self.w("Node Report : %s\n" % ("No" if L(g("NODE")) == 0 else "Yes"))
        if not blank(g("LINK")):
            self.w("Link Report : %s\n" % ("No" if L(g("LINK")) == 0 else "Yes"))
        self.w("\n")

    def make(self):
        self.setup()
        self.title()
        self.options()
        self.patterns()
        self.curves()
        if self.mode == RUN_DESIGN:
            self.design_tables()
        self.hydrograph()
        self.load_nodes()
        self.chambers()
        self.outlets()
        self.wet_wells()
        self.pipes()
        self.channel_shape()
        self.force_mains()
        self.pumps()
        self.status()
        self.controls()
        self.extra_loads()
        self.load_node_head()
        self.outlet_heads()
        self.quality_source()
        self.pipe_corrosion()
        self.flow_splits()
        self.infiltration()
        if self.mode == RUN_DESIGN:
            self.design_curves()
            self.design_offsets()
            self.design_constraints()
        self.time()
        self.report()
        self.w("\n[End]\n")
        return "".join(self.out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("iedb")
    ap.add_argument("out")
    ap.add_argument("--mode", default="steady", choices=["steady", "design", "eps"])
    ap.add_argument("--scenario")
    ap.add_argument("--full", action="store_true", help="whole network, ignore ANODE/ALINK")
    ap.add_argument("--keep-report-path", action="store_true")
    ap.add_argument("--not-pro", action="store_true")
    ap.add_argument("--epa-design-flow", action="store_true")
    a = ap.parse_args()
    mode = {"steady": RUN_STATIC, "design": RUN_DESIGN, "eps": RUN_DYNAMIC}[a.mode]
    if not os.path.isdir(a.iedb):
        sys.exit("not a folder: %s" % a.iedb)
    w = InputRun(a.iedb, mode, a.scenario, a.full, a.keep_report_path, not a.not_pro, a.epa_design_flow)
    try:
        text = w.make()
    except WriterError as e:
        sys.exit("input file not written: %s" % e)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="latin-1", newline="\n") as fh:
        fh.write(text)
    sys.stderr.write("scenario %s  time/report/option sets %s/%s/%s  active network %s\n" %
                     (w.scenario, w.timeSet, w.reportSet, w.optionSet, "yes" if w.useActive else "no"))
    for m in w.warnings:
        sys.stderr.write("warning: %s\n" % m)


if __name__ == "__main__":
    main()
