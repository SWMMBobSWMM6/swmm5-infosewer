//-----------------------------------------------------------------------------
//   isewer.c
//
//   Project:  swmm5-infosewer (EPA SWMM 5.2.4 + InfoSewer routing)
//   Purpose:  InfoSewer's sewer calculations (gravity pipes, wet wells, pump
//             stations, force mains, junction chambers) as a SWMM5 flow
//             routing option:  FLOW_ROUTING INFOSEWER
//               INFOSEWER_ANALYSIS STEADY  - InfoSewer steady state (loads by
//                                            type, peaking, normal depth, HGL)
//               INFOSEWER_ANALYSIS EPS     - InfoSewer extended period routing
//                                            (isewer_eps.c)
//
//   This file: input sections and options, the InfoSewer network built from
//   SWMM's conduits, pumps, junctions, storage nodes and outfalls, siphons and
//   pump stations, the steady-state solver, the HGL sweep shared with EPS,
//   and the SWMM state and report hooks.
//
//   The calculations follow the InfoSewer engine (HeEngPkg, MWH Soft/Innovyze,
//   1996-2011) as written up in docs/spec_steady.md, docs/spec_eps.md and
//   docs/spec_pressure.md, which cite the
//   C++ by file:line. Function comments name the engine function each one
//   reproduces. Engine quirks are reproduced on purpose so results match
//   InfoSewer bit for bit; each one is marked QUIRK.
//
//   Written by Claude (AI) for SWMM Enablement, November 2026.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"
#include "isewer.h"

#ifndef ISEWER_BUILD
#define ISEWER_BUILD "dev"
#endif

//-----------------------------------------------------------------------------
//  Options and shared data
//-----------------------------------------------------------------------------
int     IsAnalysis;            // IS_STEADY or IS_EPS
double  IsLoadMult;            // global loading multiplier
static double IsPeakK, IsPeakP;            // peakable: Q = K * Qbase^P (cfs)
static double IsPeakA, IsPeakB, IsPeakC, IsPeakD, IsPeakE;   // coverage
static char   IsPeakCurveId[MAXMSG+1];     // peaking factor curve
static char   IsCoverCurveId[MAXMSG+1];    // coverage based peaking curve
static int    IsPeakCurve, IsCoverCurve;   // curve indexes, -1 = none
double  IsMinSegment;          // minimum pipe segment (ft)
int     IsMaxSegments;         // maximum number of segments
double  IsMinTravelTime;       // minimum travel time (sec)
int     IsFlowAttenuation;     // 1 = Muskingum-Cunge
static double IsFlowTolIn;     // flow tolerance as read (flow units)
double  IsFlowTol, IsIterTol;  // flow tolerances (cfs)
double  IsFlowFactor;          // InfoSewer flow unit factor
double  IsPatternStep;         // loading pattern time step (sec)
double  IsPumpStep;            // pump station time step (sec), HMAP_StationTimeStep
int     IsLinkReport;          // InfoSewer [Report] Link : Yes (1) or No (0)
char    IsResultsFile[MAXFNAME+1];  // optional full-precision CSV

TIsLoad*    IsLoads;   int IsNloads;   static int MaxLoads;
TIsPattern* IsPats;    int IsNpats;    static int MaxPats;
static TIsCurve* IsCurves; static int IsNcurves, MaxCurves;

TIsNode* IsNode;
TIsPipe* IsPipe;
TIsPump* IsPump;
TIsWetWell* IsWell;
TIsStation* IsStation;
TIsSiphon* IsSiphon;
int      IsNN, IsNP, IsNF, IsNPM, IsNPump, IsNWell, IsNStation, IsNSiphon, IsNT, IsNLE;
int*     IsNodeOfSwmm;
int*     IsLinkOfSwmm;
int*     IsOrder;
static int  IsNorder;
static double* IsRim;        // junction rim saved before SWMM raises full depth
static double* IsFullDepth0; // node full depth as read (wet well max level)
static double* IsLinkInv;    // conduit end inverts as read (2 per link, ft)
static int  IsPumpContinuity;    // "Pump Flow Continuity": 0 upstream, else station

// pressure-network input ([INFOSEWER_CHAMBERS], _WETWELLS, _PUMPS)
typedef struct { int node; double diam; } TIsChamberIn;
typedef struct { int node; double minLevel, diam; char curve[MAXMSG+1]; } TIsWellIn;
typedef struct { int link; int type; int extra; double v[7]; char pat[MAXMSG+1];
                 double qs, qd; } TIsPumpIn;   // v: cap h0 h1 q1 h2 q2 speed
static TIsChamberIn* ChamberIn; static int NchamberIn, MaxChamberIn;
static TIsWellIn*    WellIn;    static int NwellIn, MaxWellIn;
static TIsPumpIn*    PumpIn;    static int NpumpIn, MaxPumpIn;
static double *SFmLoss;       // force main headloss (ft), per force main
static double *SPumpHead;     // pump head increase (ft), per pump
static double *SWellHead;     // wet well head (ft), per wet well

// steady-state results (internal units: cfs, ft)
static double *SQ, *SUP, *SPK, *SCV, *SCC, *SINF;  // link flow and components
static double *SDepth, *SVel, *SCrit;
static int    *SFull;
static double *SStatic;                             // node static load
static double *SHead;                               // node head
static int    *SFtype, *SPartial, *SAdjFull;
static double *SInc;                                // node load increment
static int    *IsAllocList;   // scratch

// per-link split data read from [INFOSEWER_SPLITS] / [INFOSEWER_INFILTRATION]
typedef struct { int link; double pct; } TIsSplitIn;
typedef struct { int link; int type; double rate; double count;
                 char pat[MAXMSG+1]; } TIsInfilIn;
static TIsSplitIn* SplitIn;  static int NsplitIn, MaxSplitIn;
static TIsInfilIn* InfilIn;  static int NinfilIn, MaxInfilIn;

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int    findPattern(const char* id);
static int    findCurve(const char* id);
static int    growArray(void** a, int* nmax, int n, size_t size, int first);
static double flowFactor(void);
static int    buildNetwork(void);
static int    classifyNodes(void);
static int    classifyLinks(void);
static int    buildSiphons(void);
static int    buildTopology(void);
static char*  linkId(int l);
static int    netError(int code, char* id);
static void   steadySolve(void);
static double staticGetLoading(int k, double* up, double* pk, double* cv,
                               double* cc, double* inf);
static double staticTotal(int k);
static double peakFlow(double base);
static double coverFlow(double cflow, double cunit);
static double curveY(int c, double x);
static void   splitLinks(int k, double inflow, double* Q, int* visit);
static void   multiHydSplit(int* hyd, int nhyd, double qh, double* Q,
                            int* visit);
static void   setSteadyState(void);
static void   writeSteadyReport(void);
static void   writeSteadyCsv(void);
static void   freeAll(void);

//=============================================================================
//  Options and input
//=============================================================================

char* isewer_methodName(void)
{
    static char name[80];
    sprintf(name, "INFOSEWER %s (swmm5-infosewer %s)",
            IsAnalysis == IS_EPS ? "EPS" : "STEADY", ISEWER_BUILD);
    return name;
}

void isewer_setDefaults(void)
//  Defaults are the InfoSewer engine's (he_input.cpp:457-535).
{
    IsAnalysis = IS_STEADY;
    IsLoadMult = 1.0;
    IsPeakK = 2.4;  IsPeakP = 0.89;
    IsPeakA = 5.0;  IsPeakB = 0.0;  IsPeakC = 0.2;  IsPeakD = 0.0;  IsPeakE = 1.0;
    IsPeakCurveId[0] = '\0';  IsCoverCurveId[0] = '\0';
    IsPeakCurve = IsCoverCurve = -1;
    IsMinSegment = 100.0;
    IsMaxSegments = 100;
    IsMinTravelTime = 1.0;
    IsFlowAttenuation = 0;
    IsFlowTolIn = IS_FLOW_TOL;
    IsPatternStep = 3600.0;
    IsPumpStep = 60.0;
    IsLinkReport = 1;
    IsResultsFile[0] = '\0';
    IsLoads = NULL;  IsNloads = MaxLoads = 0;
    IsPats = NULL;   IsNpats = MaxPats = 0;
    IsCurves = NULL; IsNcurves = MaxCurves = 0;
    SplitIn = NULL;  NsplitIn = MaxSplitIn = 0;
    InfilIn = NULL;  NinfilIn = MaxInfilIn = 0;
    ChamberIn = NULL;  NchamberIn = MaxChamberIn = 0;
    WellIn = NULL;     NwellIn = MaxWellIn = 0;
    PumpIn = NULL;     NpumpIn = MaxPumpIn = 0;
    IsPumpContinuity = 0;
    IsNode = NULL;  IsPipe = NULL;  IsOrder = NULL;  IsRim = NULL;  IsFullDepth0 = NULL;
    IsPump = NULL;  IsWell = NULL;  IsStation = NULL;  IsSiphon = NULL;
    IsNodeOfSwmm = NULL;  IsLinkOfSwmm = NULL;
    IsNN = IsNP = IsNF = IsNPM = IsNPump = IsNWell = IsNStation = IsNSiphon = IsNT = IsNLE = 0;
    SFmLoss = SPumpHead = SWellHead = NULL;
    IsAllocList = NULL;
    SQ = SUP = SPK = SCV = SCC = SINF = SDepth = SVel = SCrit = NULL;
    SStatic = SHead = SInc = NULL;
    SFull = SFtype = SPartial = SAdjFull = NULL;
}

static int growArray(void** a, int* nmax, int n, size_t size, int first)
{
    void* p;
    if ( n < *nmax ) return 0;
    *nmax = *nmax ? 2 * (*nmax) : first;
    p = realloc(*a, (*nmax) * size);
    if ( !p ) return ERR_MEMORY;
    *a = p;
    return 0;
}

int isewer_readOption(char* s1, char* s2)
//  INFOSEWER_<name> <value> lines from [OPTIONS]
{
    double x;
    char*  k = s1 + 10;
    if ( match(k, "ANALYSIS") )
    {
        if ( match(s2, "STEADY") ) IsAnalysis = IS_STEADY;
        else if ( match(s2, "EPS") || match(s2, "DYNAMIC") ) IsAnalysis = IS_EPS;
        else return error_setInpError(ERR_KEYWORD, s2);
        return 0;
    }
    if ( match(k, "PEAK_CURVE") )
    {
        sstrncpy(IsPeakCurveId, s2, MAXMSG);
        return 0;
    }
    if ( match(k, "COVER_CURVE") )
    {
        sstrncpy(IsCoverCurveId, s2, MAXMSG);
        return 0;
    }
    if ( match(k, "LINK_REPORT") )
    {
        if ( match(s2, "NO") ) IsLinkReport = 0;
        else if ( match(s2, "YES") ) IsLinkReport = 1;
        else return error_setInpError(ERR_KEYWORD, s2);
        return 0;
    }
    if ( match(k, "RESULTS_FILE") )
    {
        sstrncpy(IsResultsFile, s2, MAXFNAME);
        return 0;
    }
    if ( match(k, "PATTERN_STEP") || match(k, "PUMP_STEP") )
    {
        DateTime t;
        double* step = match(k, "PATTERN_STEP") ? &IsPatternStep : &IsPumpStep;
        if ( getDouble(s2, &x) ) { *step = x; return 0; }
        if ( !datetime_strToTime(s2, &t) ) return error_setInpError(ERR_DATETIME, s2);
        *step = floor(t * SECperDAY + 0.5);
        return 0;
    }
    if ( !getDouble(s2, &x) ) return error_setInpError(ERR_NUMBER, s2);
    if      ( match(k, "LOAD_MULT") )          IsLoadMult = x;
    else if ( match(k, "PEAK_K") )             IsPeakK = x;
    else if ( match(k, "PEAK_P") )             IsPeakP = x;
    else if ( match(k, "PEAK_A") )             IsPeakA = x;
    else if ( match(k, "PEAK_B") )             IsPeakB = x;
    else if ( match(k, "PEAK_C") )             IsPeakC = x;
    else if ( match(k, "PEAK_D") )             IsPeakD = x;
    else if ( match(k, "PEAK_E") )             IsPeakE = x;
    else if ( match(k, "MIN_SEGMENT") )        IsMinSegment = x;
    else if ( match(k, "MAX_SEGMENTS") )       IsMaxSegments = (int)x;
    else if ( match(k, "MIN_TRAVEL_TIME") )    IsMinTravelTime = x;
    else if ( match(k, "FLOW_ATTENUATION") )   IsFlowAttenuation = (x >= IS_ZERO_TOL);
    else if ( match(k, "FLOW_TOLERANCE") )     IsFlowTolIn = x;
    else if ( match(k, "PUMP_CONTINUITY") )    IsPumpContinuity = ((int)x != 0);  // QUIRK he_input.cpp:9712
    else return error_setInpError(ERR_KEYWORD, s1);
    return 0;
}

int isewer_readLoad(char* tok[], int ntoks)
//  [INFOSEWER_LOADS]  Node  Base  Type(UNPEAK|PEAK|COVER)  Pattern|*  [CoverUnits]
{
    int    j, type;
    double base, cover = 0.0;
    TIsLoad* ld;
    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");
    j = project_findObject(NODE, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);
    if ( !getDouble(tok[1], &base) ) return error_setInpError(ERR_NUMBER, tok[1]);
    if      ( match(tok[2], "UNPEAK") ) type = IS_UNPEAK;
    else if ( match(tok[2], "PEAK") )   type = IS_PEAK;
    else if ( match(tok[2], "COVER") )  type = IS_COVER;
    else return error_setInpError(ERR_KEYWORD, tok[2]);
    if ( ntoks > 4 && !getDouble(tok[4], &cover) )
        return error_setInpError(ERR_NUMBER, tok[4]);
    if ( growArray((void**)&IsLoads, &MaxLoads, IsNloads, sizeof(TIsLoad), 256) )
        return ERR_MEMORY;
    ld = &IsLoads[IsNloads++];
    ld->node = j;  ld->baseIn = base;  ld->base = 0.0;  ld->type = type;
    ld->cover = (type == IS_COVER) ? cover : 0.0;
    ld->pat = -1;
    ld->pattern[0] = '\0';
    if ( ntoks > 3 && strcmp(tok[3], "*") != 0 ) sstrncpy(ld->pattern, tok[3], MAXMSG);
    return 0;
}

int isewer_readPattern(char* tok[], int ntoks)
//  [INFOSEWER_PATTERNS]  Name  [STEPWISE|LINEAR]  f1 f2 ...
//  (lines with the same name continue the pattern; default STEPWISE)
{
    int i, p, i0 = 1;
    double x;
    TIsPattern* pt;
    if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
    p = findPattern(tok[0]);
    if ( p < 0 )
    {
        if ( growArray((void**)&IsPats, &MaxPats, IsNpats, sizeof(TIsPattern), 16) )
            return ERR_MEMORY;
        p = IsNpats++;
        pt = &IsPats[p];
        sstrncpy(pt->id, tok[0], MAXMSG);
        pt->f = NULL;  pt->n = pt->nmax = 0;  pt->linear = 0;
    }
    pt = &IsPats[p];
    if ( match(tok[1], "LINEAR") )        { pt->linear = 1; i0 = 2; }
    else if ( match(tok[1], "STEPWISE") ) { pt->linear = 0; i0 = 2; }
    for ( i = i0; i < ntoks; i++ )
    {
        if ( !getDouble(tok[i], &x) ) return error_setInpError(ERR_NUMBER, tok[i]);
        if ( growArray((void**)&pt->f, &pt->nmax, pt->n, sizeof(double), 24) )
            return ERR_MEMORY;
        pt->f[pt->n++] = x;
    }
    return 0;
}

int isewer_readCurve(char* tok[], int ntoks)
//  [INFOSEWER_CURVES]  Name  X  Y   (lines with the same name continue it)
//  A curve named by INFOSEWER_PEAK_CURVE gives the peaked flow for a base
//  flow (both in flow units); one named by INFOSEWER_COVER_CURVE gives the
//  peaking factor for a coverage count.
{
    int c, i;
    double x[2];
    TIsCurve* cv;
    if ( ntoks < 3 || (ntoks - 1) % 2 ) return error_setInpError(ERR_ITEMS, "");
    c = findCurve(tok[0]);
    if ( c < 0 )
    {
        if ( growArray((void**)&IsCurves, &MaxCurves, IsNcurves, sizeof(TIsCurve), 8) )
            return ERR_MEMORY;
        c = IsNcurves++;
        cv = &IsCurves[c];
        sstrncpy(cv->id, tok[0], MAXMSG);
        cv->x = cv->y = NULL;  cv->n = cv->nmax = 0;
    }
    cv = &IsCurves[c];
    for ( i = 1; i + 1 < ntoks; i += 2 )
    {
        if ( !getDouble(tok[i], &x[0]) )   return error_setInpError(ERR_NUMBER, tok[i]);
        if ( !getDouble(tok[i+1], &x[1]) ) return error_setInpError(ERR_NUMBER, tok[i+1]);
        if ( cv->n == cv->nmax )
        {
            int nmax = cv->nmax ? 2 * cv->nmax : 8;
            double* px = (double*) realloc(cv->x, nmax * sizeof(double));
            double* py;
            if ( !px ) return ERR_MEMORY;
            cv->x = px;
            py = (double*) realloc(cv->y, nmax * sizeof(double));
            if ( !py ) return ERR_MEMORY;
            cv->y = py;
            cv->nmax = nmax;
        }
        cv->x[cv->n] = x[0];  cv->y[cv->n] = x[1];  cv->n++;
    }
    return 0;
}

int isewer_readSplit(char* tok[], int ntoks)
//  [INFOSEWER_SPLITS]  Link  Percent
//  Percent share of the flow at a split manhole; a pipe without a percent
//  (or with 0) takes a hydraulic share (InfoSewer [Flow Splits] type 1).
{
    int j;
    double pct;
    if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
    j = project_findObject(LINK, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);
    if ( !getDouble(tok[1], &pct) || pct < 0.0 || pct > 100.0 )
        return error_setInpError(ERR_NUMBER, tok[1]);
    if ( growArray((void**)&SplitIn, &MaxSplitIn, NsplitIn, sizeof(TIsSplitIn), 16) )
        return ERR_MEMORY;
    SplitIn[NsplitIn].link = j;
    SplitIn[NsplitIn].pct = pct;
    NsplitIn++;
    return 0;
}

int isewer_readInfil(char* tok[], int ntoks)
//  [INFOSEWER_INFILTRATION]  Link  Type  Rate  [Count | Pattern]
//    Type LENGTH       rate x length              (flow units per ft or m)
//         DIAMETER     rate x diameter x length   (per in-ft or mm-m)
//         SURFACE      rate x pi x diameter x length
//         COUNT        rate x count
//         PATTERN      rate x pattern factor
//  The result is multiplied by the number of barrels.
{
    int j, t;
    double rate, count = 0.0;
    TIsInfilIn* in;
    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");
    j = project_findObject(LINK, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);
    if      ( match(tok[1], "NONE") )     t = IS_INF_NONE;
    else if ( match(tok[1], "LENGTH") )   t = IS_INF_LEN;
    else if ( match(tok[1], "DIAMETER") ) t = IS_INF_LENDIA;
    else if ( match(tok[1], "SURFACE") )  t = IS_INF_SURF;
    else if ( match(tok[1], "COUNT") )    t = IS_INF_COUNT;
    else if ( match(tok[1], "PATTERN") )  t = IS_INF_PATTERN;
    else return error_setInpError(ERR_KEYWORD, tok[1]);
    if ( !getDouble(tok[2], &rate) ) return error_setInpError(ERR_NUMBER, tok[2]);
    if ( t == IS_INF_COUNT )
    {
        if ( ntoks < 4 || !getDouble(tok[3], &count) )
            return error_setInpError(ERR_NUMBER, ntoks > 3 ? tok[3] : "");
    }
    if ( t == IS_INF_PATTERN && ntoks < 4 ) return error_setInpError(ERR_ITEMS, "");
    if ( growArray((void**)&InfilIn, &MaxInfilIn, NinfilIn, sizeof(TIsInfilIn), 16) )
        return ERR_MEMORY;
    in = &InfilIn[NinfilIn++];
    in->link = j;  in->type = t;  in->rate = rate;  in->count = count;
    in->pat[0] = '\0';
    if ( t == IS_INF_PATTERN ) sstrncpy(in->pat, tok[3], MAXMSG);
    return 0;
}

int isewer_readChamber(char* tok[], int ntoks)
//  [INFOSEWER_CHAMBERS]  Node  [Diameter]
//  Marks a junction as an InfoSewer junction chamber (force main node).
{
    int j;
    double d = 0.0;
    if ( ntoks < 1 ) return error_setInpError(ERR_ITEMS, "");
    j = project_findObject(NODE, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);
    if ( ntoks > 1 && !getDouble(tok[1], &d) ) return error_setInpError(ERR_NUMBER, tok[1]);
    if ( growArray((void**)&ChamberIn, &MaxChamberIn, NchamberIn, sizeof(TIsChamberIn), 16) )
        return ERR_MEMORY;
    ChamberIn[NchamberIn].node = j;
    ChamberIn[NchamberIn].diam = d;
    NchamberIn++;
    return 0;
}

int isewer_readWetWell(char* tok[], int ntoks)
//  [INFOSEWER_WETWELLS]  Node  MinLevel  [VolumeCurve|*  [Diameter]]
//  For a SWMM storage node: bottom = invert, maximum level = max depth,
//  initial level = initial depth. A volume curve (level vs volume,
//  [INFOSEWER_CURVES]) makes it a variable-area well; otherwise it is a
//  cylinder of the given diameter (default: from the CYLINDRICAL area).
{
    int j;
    double x;
    TIsWellIn* w;
    if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
    j = project_findObject(NODE, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);
    if ( !getDouble(tok[1], &x) ) return error_setInpError(ERR_NUMBER, tok[1]);
    if ( growArray((void**)&WellIn, &MaxWellIn, NwellIn, sizeof(TIsWellIn), 16) )
        return ERR_MEMORY;
    w = &WellIn[NwellIn++];
    w->node = j;  w->minLevel = x;  w->curve[0] = '\0';  w->diam = -1.0;
    if ( ntoks > 2 && strcmp(tok[2], "*") != 0 ) sstrncpy(w->curve, tok[2], MAXMSG);
    if ( ntoks > 3 )
    {
        if ( !getDouble(tok[3], &w->diam) || w->diam <= 0.0 )
            return error_setInpError(ERR_NUMBER, tok[3]);
    }
    return 0;
}

int isewer_readPump(char* tok[], int ntoks)
//  [INFOSEWER_PUMPS]  Link  CAPACITY|POINT|THREEPOINT  Extra  Capacity  Head0
//                     Head1  Flow1  Head2  Flow2  Speed  Pattern|-  StaticFlow  BaseFlow
//  InfoSewer's own pump definition ([Pumps] fields 3-14); flows in flow
//  units, heads in ft (m). The SWMM pump curve is for stock SWMM only.
{
    int j, i;
    TIsPumpIn* pp;
    if ( ntoks < 13 ) return error_setInpError(ERR_ITEMS, "");
    j = project_findObject(LINK, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);
    if ( growArray((void**)&PumpIn, &MaxPumpIn, NpumpIn, sizeof(TIsPumpIn), 16) )
        return ERR_MEMORY;
    pp = &PumpIn[NpumpIn];
    pp->link = j;
    if      ( match(tok[1], "CAPACITY") )   pp->type = IS_PUMP_CAPACITY;
    else if ( match(tok[1], "POINT") )      pp->type = IS_PUMP_POINT;
    else if ( match(tok[1], "THREEPOINT") ) pp->type = IS_PUMP_3POINT;
    else return error_setInpError(ERR_KEYWORD, tok[1]);
    pp->extra = atoi(tok[2]);
    for ( i = 0; i < 7; i++ )
        if ( !getDouble(tok[3+i], &pp->v[i]) ) return error_setInpError(ERR_NUMBER, tok[3+i]);
    pp->pat[0] = '\0';
    if ( strcmp(tok[10], "-") != 0 && strcmp(tok[10], "*") != 0 ) sstrncpy(pp->pat, tok[10], MAXMSG);
    if ( !getDouble(tok[11], &pp->qs) ) return error_setInpError(ERR_NUMBER, tok[11]);
    if ( !getDouble(tok[12], &pp->qd) ) return error_setInpError(ERR_NUMBER, tok[12]);
    NpumpIn++;
    return 0;
}

static int findPattern(const char* id)
{
    int p;
    for ( p = 0; p < IsNpats; p++ ) if ( strcomp((char*)id, IsPats[p].id) ) return p;
    return -1;
}

static int findCurve(const char* id)
{
    int c;
    for ( c = 0; c < IsNcurves; c++ ) if ( strcomp((char*)id, IsCurves[c].id) ) return c;
    return -1;
}

void isewer_saveRims(void)
//  Called by project_validate before link_validate raises a node's full
//  depth to the crown of its highest conduit: InfoSewer caps unlocked
//  manhole heads at the rim the user gave, and a wet well's maximum level is
//  the storage depth the user gave, so keep both. Also keeps each conduit's
//  end inverts as read: link_validate drops an offset below the node invert
//  (WARNING 03), but InfoSewer uses the pipe inverts as given, even below a
//  wet well bottom.
{
    int i, j;
    FREE(IsRim);
    FREE(IsFullDepth0);
    FREE(IsLinkInv);
    IsRim = (double*) calloc(Nobjects[NODE] + 1, sizeof(double));
    IsFullDepth0 = (double*) calloc(Nobjects[NODE] + 1, sizeof(double));
    IsLinkInv = (double*) calloc(2 * Nobjects[LINK] + 2, sizeof(double));
    if ( !IsRim || !IsFullDepth0 || !IsLinkInv ) { report_writeErrorMsg(ERR_MEMORY, ""); return; }
    for ( j = 0; j < Nobjects[NODE]; j++ )
    {
        IsRim[j] = Node[j].invertElev + Node[j].fullDepth;
        IsFullDepth0[j] = Node[j].fullDepth;
    }
    for ( i = 0; i < Nobjects[LINK]; i++ )
    {
        double e1 = Node[Link[i].node1].invertElev, e2 = Node[Link[i].node2].invertElev;
        double o1 = Link[i].offset1, o2 = Link[i].offset2;
        if ( LinkOffsets == ELEV_OFFSET )
        {
            IsLinkInv[2*i]   = (o1 <= MISSING) ? e1 : o1;
            IsLinkInv[2*i+1] = (o2 <= MISSING) ? e2 : o2;
        }
        else
        {
            IsLinkInv[2*i]   = e1 + o1;
            IsLinkInv[2*i+1] = e2 + o2;
        }
    }
}

static double flowFactor(void)
//  InfoSewer's flow unit factors (he_define.h); they differ from SWMM's Qcf
//  in the 6th digit, and loads must be divided by InfoSewer's to match it.
{
    switch ( FlowUnits )
    {
      case CFS: return 1.0;
      case GPM: return 448.831169;
      case MGD: return 0.64632;
      case LPS: return 28.316847;
      case MLD: return 2.446575;
      default:  return UCF(FLOW);     // CMS: not an InfoSewer unit
    }
}

//=============================================================================
//  Circular pipe hydraulics (he_hydrun.cpp, he_data.cpp)
//=============================================================================

double isewer_wetArea(double d, double r)
//  HHU_GetWetArea he_hydrun.cpp:2995 (3.1416, not pi)
{
    double angle = (d - r) / r;
    double angle2 = (2.0 * 3.1416) - (2.0 * acos(angle));
    double dp = 2.0 * r;
    return 0.125 * (angle2 - sin(angle2)) * (dp * dp);
}

double isewer_wetPerimeter(double d, double r)
//  HHU_GetWetPerimeter he_hydrun.cpp:3945
{
    double angle = (d - r) / r;
    double angle2 = (2.0 * 3.1416) - (2.0 * acos(angle));
    return (0.5 * angle2) * (2.0 * r);
}

double isewer_wetWidth(double d, double r)
//  HHU_GetWetWidth he_hydrun.cpp:2873
{
    return 2.0 * sqrt(d * (2 * r - d));
}

static double getFlow(double d, double r, double slope, double n)
//  HHU_GetFlow he_hydrun.cpp:3767-3824 (1.49 and 0.6667)
{
    double A, P, R;
    if ( slope == 0.0 ) slope = 0.00001;
    A = isewer_wetArea(d, r);
    P = isewer_wetPerimeter(d, r);
    R = A / P;
    return (1.49 / n) * A * pow(R, 0.6667) * sqrt(slope);
}

double isewer_fullFlowBarrel(int i)
//  HHU_GetFullFlow he_hydrun.cpp:3588-3621, per barrel.
//  QUIRK: an adverse pipe has full flow -1 (per barrel).
{
    TIsPipe* p = &IsPipe[i];
    double D = p->dia, r = 0.5 * p->dia, s = p->slope;
    if ( s == 0.0 ) s = 0.00001;
    if ( s < 0 ) return -1.0;
    return getFlow(D, r, s, p->n);
}

double isewer_linkFlow(int i, double d, double r, double n)
//  HHU_GetLinkFlow he_hydrun.cpp:3312 (circular). r < 0 and n <= 0 mean
//  "the pipe's own".
{
    TIsPipe* p = &IsPipe[i];
    double s;
    if ( p->adverse ) return p->fullFlow;
    s = p->slope;
    if ( s == 0.0 ) s = 0.00001;
    if ( r < 0.0 ) r = 0.5 * p->dia;
    if ( n <= 0.0 ) n = p->n;
    return getFlow(d, r, s, n);
}

double isewer_flowDepth(int i, double q)
//  HHU_GetFlowDepth he_hydrun.cpp:2164 (circular): bisection on [0, D]
//  starting at D/2, stopping on a 1e-4 relative flow error.
{
    TIsPipe* p = &IsPipe[i];
    double D, n, r, dL, qL, dU, qU, depth, d1, d2, av;
    long   iter = 0;
    if ( fabs(q) < IS_FLOW_TOL ) return 0.0;
    D = p->dia;  n = p->n;  r = 0.5 * D;
    dL = 0.0;  qL = 0.0;
    dU = D;
    qU = isewer_linkFlow(i, dU, r, n);
    if ( qU <= q ) return D;
    depth = 0.5 * (dU + dL);
    d1 = q - qL;
    d2 = q - qU;
    av = q - isewer_linkFlow(i, depth, r, n);
    for (;;)
    {
        if ( av * d2 <= 0.0 )
        {
            dL = depth;
            qL = isewer_linkFlow(i, dL, r, n);
            d1 = q - qL;
        }
        else if ( av * d1 <= 0.0 )
        {
            dU = depth;
            qU = isewer_linkFlow(i, dU, r, n);
            d2 = q - qU;
        }
        depth = 0.5 * (dL + dU);
        av = q - isewer_linkFlow(i, depth, r, n);
        if ( !((fabs(av) / q) >= 0.0001) ) break;
        if ( ++iter > 100000 ) break;           // safety only; never reached
    }
    if ( p->slopeType != IS_SLOPE_STANDARD )    // he_hydrun.cpp:2299-2305
    {
        depth = depth + 0.5 * p->length * (-1.0 * p->slope);
        depth = MIN(depth, D);
    }
    return depth;
}

double isewer_criticalDepth(int i, double q, double D)
//  HHU_GetCriticalDepth he_hydrun.cpp:2586 (circular 2642-2679).
//  QUIRK: the bisection's lower update is 0.5*(d1 - d0), not +; g = 9.81/0.3048.
{
    double r, d0, d1, depth, ref, A, T, Fr;
    (void)i;
    if ( fabs(q) < IsFlowTol ) return 0.0;
    r = 0.5 * D;
    d0 = 0.0;  d1 = D;
    depth = 0.5 * (d1 - d0);
    ref = 0.0;
    A = isewer_wetArea(depth, r);
    if ( A == 0.0 ) return 0.0;
    Fr = (q / A) / sqrt(A / isewer_wetWidth(depth, r) * 9.81 * 1 / 0.3048);
    while ( fabs(depth - ref) > IS_DEPTH_TOL )
    {
        ref = depth;
        if ( Fr < 1 )
        {
            d1 = depth;
            depth = 0.5 * (d1 - d0);
        }
        else
        {
            d0 = depth;
            depth = d0 + 0.5 * (d1 - d0);
        }
        A = isewer_wetArea(depth, r);
        if ( A == 0.0 ) return 0.0;
        T = isewer_wetWidth(depth, r);
        Fr = (q / A) / sqrt(A / T * 9.81 * 1 / 0.3048);
    }
    return depth;
}

double isewer_pipeVelocity(double d, double r, double q)
//  HHU_GetPipeVelocity he_hydrun.cpp:3138 (checkCutoff = false)
{
    double A = isewer_wetArea(d, r);
    if ( A < IS_WETAREA_TOL ) return 0.0;
    if ( q < IS_ZERO_TOL ) return 0.0;
    return q / A;
}

double isewer_headloss(int i, double q)
//  HD_PipeData::F_GetHeadloss he_data.cpp:5480 (circular).
//  QUIRK: 4.66 (i.e. 1.486) here vs 1.49 in Manning; uses the TOTAL flow of
//  parallel barrels.
{
    TIsPipe* p = &IsPipe[i];
    double varN, varD, fp, dp, qp;
    if ( p->n < IS_MAX_MANNING )
    {
        varN = 4.66 * p->n * p->n * p->length * q * q;
        varD = pow(p->dia, 5.33);
        return fabs(varD) > IS_SYS_ZERO ? varN / varD : 0.0;
    }
    fp = pow(p->n, 1.852);
    dp = pow(p->dia, 4.871);
    qp = pow(q, 1.852);
    varN = 4.727 * p->length * qp;
    varD = fp * dp;
    return fabs(varD) > IS_SYS_ZERO ? varN / varD : 0.0;
}

double isewer_infiltration(int i)
//  HD_PipeData::F_GetInfiltrateFlow he_data.cpp:5420 at time 0 (steady);
//  isewer_eps.c evaluates PATTERN infiltration over time itself.
{
    TIsPipe* p = &IsPipe[i];
    double f;
    switch ( p->infType )
    {
      case IS_INF_NONE:    return 0.0;
      case IS_INF_LEN:     f = p->infRate * p->length; break;
      case IS_INF_LENDIA:  f = p->infRate * p->dia * p->length; break;
      case IS_INF_SURF:    f = p->infRate * IS_PI_DEF * p->dia * p->length; break;
      case IS_INF_COUNT:   f = p->infRate * p->infCount; break;
      case IS_INF_PATTERN:
          f = p->infRate * ((p->infPat >= 0 && IsPats[p->infPat].n > 0) ?
                            IsPats[p->infPat].f[0] : 1.0);
          break;
      default:             f = 0.0;
    }
    return f * p->count;
}

static double curveY(int c, double x)
//  HD_CurveData::F_GetInterpolatedY he_data.cpp:10224.
//  QUIRK: extrapolation past the last point is unbounded.
{
    double *X, *Y, dx, dy;
    int n, i, last;
    if ( c < 0 ) return 0.0;
    X = IsCurves[c].x;  Y = IsCurves[c].y;  n = IsCurves[c].n;
    if ( n <= 0 || x < 0.0 ) return 0.0;
    last = n - 1;
    if ( x <= X[0] )
    {
        if ( X[0] == 0.0 ) return 0.0;
        return x * Y[0] / X[0];
    }
    else if ( x >= X[last] )
    {
        if ( n == 1 ) return Y[0];
        dx = X[last] - X[last-1];  dy = Y[last] - Y[last-1];
        if ( fabs(dx) < IS_FLOW_TOL ) return Y[last];
        return Y[last] + (x - X[last]) * dy / dx;
    }
    if ( n == 1 ) return Y[0];
    for ( i = 1; i < n; i++ )
    {
        if ( X[i] >= x )
        {
            dx = X[i] - X[i-1];  dy = Y[i] - Y[i-1];
            if ( fabs(dx) < IS_FLOW_TOL ) return Y[i];
            return Y[i] - (X[i] - x) * dy / dx;
        }
    }
    return 0.0;
}

//=============================================================================
//  Network
//=============================================================================

double isewer_wellVolumeFromLevel(int w, double level)
//  HU_GetVarAreaWellVolumeFromLevel he_util.cpp:1661 (curve wells only)
{
    TIsCurve* c;
    int i;
    if ( IsWell[w].curve < 0 ) return 0.0;
    c = &IsCurves[IsWell[w].curve];
    if ( c->n <= 0 ) return 0.0;
    if ( c->n == 1 ) return c->y[0];
    for ( i = 1; i < c->n; i++ )
        if ( c->x[i] > level )
        {
            double x1 = c->x[i-1], x2 = c->x[i], y1 = c->y[i-1], y2 = c->y[i];
            double deltaX = level - x1;
            return y1 + ((y2 - y1)/(x2 - x1))*(deltaX);
        }
    return c->y[c->n - 1];
}

double isewer_wellLevelFromVolume(int w, double volume)
//  HU_GetWellLevelFromVolume he_util.cpp:1539
{
    TIsWetWell* ww = &IsWell[w];
    if ( volume <= 0.0 ) return 0.0;
    if ( ww->curve < 0 )
    {
        double radius = 0.5*ww->diameter;
        double area = IS_PI_DEF*radius*radius;
        return volume/area;
    }
    else
    {
        // HU_GetVarAreaWellLevelFromVolume he_util.cpp:1708
        TIsCurve* c = &IsCurves[ww->curve];
        int i;
        if ( c->n <= 0 ) return 0.0;
        if ( c->n == 1 ) return c->x[0];
        for ( i = 1; i < c->n; i++ )
            if ( c->y[i] > volume )
            {
                double x1 = c->x[i-1], x2 = c->x[i], y1 = c->y[i-1], y2 = c->y[i];
                double deltaY = volume - y1;
                return x1 + ((x2 - x1)/(y2 - y1))*(deltaY);
            }
        return c->x[c->n - 1];
    }
}

void isewer_init(void)
//  Called by flowrout_init in place of validateTreeLayout. Builds the
//  InfoSewer network and, for steady state, solves it; for EPS it routes the
//  whole run. SWMM's initNodes/initLinks then compute volumes from the
//  depths set here.
{
    if ( ErrorCode ) return;
    IsFlowFactor = flowFactor();
    if ( buildNetwork() ) return;
    if ( buildTopology() ) return;
    if ( IsAnalysis == IS_STEADY )
    {
        steadySolve();
        if ( ErrorCode ) return;
        setSteadyState();
        writeSteadyReport();
        if ( IsResultsFile[0] ) writeSteadyCsv();
    }
    else
    {
        if ( isewer_epsInit() ) return;
        isewer_epsSetState(0.0);
        isewer_epsReport();
    }
}

static int netError(int code, char* id)
{
    report_writeErrorMsg(code, id);
    return ErrorCode ? ErrorCode : code;
}

static char* linkId(int l)
//  ID of a topology link: the SWMM conduit, or a station's first pump.
{
    if ( l < IsNPM ) return Link[IsPipe[l].swmm].ID;
    return Link[IsPump[IsStation[IsPipe[l].station].pump[0]].swmm].ID;
}

static int classifyNodes(void)
//  InfoSewer node order (he_input.cpp:740): loading manholes, junction
//  chambers, outlets, wet wells, each in SWMM input order.
{
    int j, k, m, N = Nobjects[NODE], pass, kind;
    int *isChamber;

    isChamber = (int*) calloc(N + 1, sizeof(int));
    if ( !isChamber ) return netError(ERR_MEMORY, "");
    for ( m = 0; m < NchamberIn; m++ ) isChamber[ChamberIn[m].node] = 1;
    for ( j = 0; j < N; j++ )
    {
        if ( Node[j].type == JUNCTION || Node[j].type == OUTFALL ) continue;
        if ( Node[j].type == STORAGE ) continue;
        FREE(isChamber);
        return netError(ERR_ISEWER_NODE, Node[j].ID);
    }
    IsNN = N;
    IsNWell = 0;
    k = 0;
    for ( pass = 0; pass < 4; pass++ )
    {
        for ( j = 0; j < N; j++ )
        {
            if ( Node[j].type == JUNCTION ) kind = isChamber[j] ? IS_NODE_CHAMBER : IS_NODE_LOAD;
            else if ( Node[j].type == OUTFALL ) kind = IS_NODE_OUTLET;
            else kind = IS_NODE_WETWELL;
            if ( kind != pass ) continue;
            IsNode[k].kind = kind;
            IsNode[k].swmm = j;
            IsNode[k].isOutlet = (kind == IS_NODE_OUTLET);
            IsNode[k].wetwell = -1;
            if ( kind == IS_NODE_WETWELL ) IsNode[k].wetwell = IsNWell++;
            IsNodeOfSwmm[j] = k;
            k++;
        }
    }
    FREE(isChamber);
    return 0;
}

static int classifyLinks(void)
//  InfoSewer link order (he_input.cpp:747): gravity pipes, force mains (the
//  "pipe mains"), then pump elements; station links share the pump element
//  index range (QUIRK: hm_pipeMainCount + station == hm_pipeMainCount + pump).
{
    int i, k, L = Nobjects[LINK], pass;
    IsNP = IsNF = IsNPump = 0;
    for ( i = 0; i < L; i++ )
    {
        if ( Link[i].type == CONDUIT )
        {
            if ( Link[i].xsect.type == FORCE_MAIN ) IsNF++;
            else if ( Link[i].xsect.type == CIRCULAR ) IsNP++;
            else return netError(ERR_ISEWER_LINK, Link[i].ID);
        }
        else if ( Link[i].type == PUMP ) IsNPump++;
        else return netError(ERR_ISEWER_LINK, Link[i].ID);
    }
    if ( IsNF > 0 && ForceMainEqn != H_W ) return netError(ERR_ISEWER_FM, "");
    IsNPM = IsNP + IsNF;
    k = 0;
    for ( pass = 0; pass < 3; pass++ )
    {
        for ( i = 0; i < L; i++ )
        {
            int kind;
            if ( Link[i].type == PUMP ) kind = 2;
            else kind = (Link[i].xsect.type == FORCE_MAIN) ? 1 : 0;
            if ( kind != pass ) continue;
            IsLinkOfSwmm[i] = k++;
        }
    }
    return 0;
}

static int buildNetwork(void)
//  InfoSewer nodes and links from SWMM junctions, outfalls, storage nodes,
//  conduits and pumps. Mirrors he_input.cpp: Fv_ReadPipe (slope, bottom
//  reference), Fv_ReadForceMain, Fv_ReadWetWell, Fv_ReadPump, the siphon and
//  station set-ups, Fv_UnitConversion, Fv_SetPumpCoeffData and the option
//  conversions of he_input.cpp:12143-12339.
{
    int    i, j, k, m, n1, n2, N = Nobjects[NODE], L = Nobjects[LINK];
    double lenF, diaF;

    if ( UnitSystem == US ) { lenF = 1.0;    diaF = 12.0; }
    else                    { lenF = 0.3048; diaF = 304.8; }

    IsNode = (TIsNode*) calloc(N + 1, sizeof(TIsNode));
    IsOrder = (int*) calloc(N + 1, sizeof(int));
    IsNodeOfSwmm = (int*) calloc(N + 1, sizeof(int));
    IsLinkOfSwmm = (int*) calloc(L + 1, sizeof(int));
    if ( !IsNode || !IsOrder || !IsNodeOfSwmm || !IsLinkOfSwmm )
        return netError(ERR_MEMORY, "");
    if ( classifyNodes() ) return ErrorCode;
    if ( classifyLinks() ) return ErrorCode;

    // --- stations: pumps grouped by (wet well, to-node) (Fv_SetWorkStationData)
    IsPump = (TIsPump*) calloc(IsNPump + 1, sizeof(TIsPump));
    IsWell = (TIsWetWell*) calloc(IsNWell + 1, sizeof(TIsWetWell));
    IsStation = (TIsStation*) calloc(IsNPump + 1, sizeof(TIsStation));
    if ( !IsPump || !IsWell || !IsStation ) return netError(ERR_MEMORY, "");
    for ( i = 0; i < L; i++ )
    {
        int p;
        if ( Link[i].type != PUMP ) continue;
        p = IsLinkOfSwmm[i] - IsNPM;
        IsPump[p].swmm = i;
        IsPump[p].frm = IsNodeOfSwmm[Link[i].node1];
        IsPump[p].to  = IsNodeOfSwmm[Link[i].node2];
        IsPump[p].station = -1;
        if ( IsNode[IsPump[p].frm].kind != IS_NODE_WETWELL ) return netError(ERR_ISEWER_PUMP, Link[i].ID);
        if ( IsNode[IsPump[p].to].kind == IS_NODE_WETWELL ||
             IsNode[IsPump[p].to].kind == IS_NODE_OUTLET ) return netError(ERR_ISEWER_PUMP, Link[i].ID);
        IsPump[p].wetwell = IsNode[IsPump[p].frm].wetwell;
    }
    IsNStation = 0;
    for ( i = 0; i < IsNPump; i++ )
    {
        int s, cnt = 0;
        if ( IsPump[i].station >= 0 ) continue;
        s = IsNStation++;
        IsPump[i].station = s;
        for ( j = 0; j < IsNPump; j++ )
        {
            if ( j == i || IsPump[j].station >= 0 ) continue;
            if ( IsPump[j].frm == IsPump[i].frm && IsPump[j].to == IsPump[i].to ) IsPump[j].station = s;
        }
        for ( j = 0; j < IsNPump; j++ ) if ( IsPump[j].station == s ) cnt++;
        IsStation[s].pump = (int*) calloc(cnt + 1, sizeof(int));
        if ( !IsStation[s].pump ) return netError(ERR_MEMORY, "");
        IsStation[s].npump = 0;
    }
    for ( i = 0; i < IsNPump; i++ )
    {
        TIsStation* st = &IsStation[IsPump[i].station];
        st->pump[st->npump++] = i;
    }
    IsNT  = IsNPM + IsNStation;
    IsNLE = IsNPM + MAX(IsNPump, IsNStation);
    IsPipe = (TIsPipe*) calloc(IsNT + 1, sizeof(TIsPipe));
    IsAllocList = (int*) calloc(8 * (IsNT + 1) + IsNloads + 1, sizeof(int));  // conn, down, up, uph: 2 per link end
    if ( !IsPipe || !IsAllocList ) return netError(ERR_MEMORY, "");

    // --- nodes
    for ( k = 0; k < IsNN; k++ )
    {
        TIsNode* nd = &IsNode[k];
        j = nd->swmm;
        nd->rim = IsRim ? IsRim[j] : Node[j].invertElev + Node[j].fullDepth;
        nd->cap = nd->rim + Node[j].surDepth;
        nd->bottomRef = IS_BIG;
        nd->submerge = -IS_BIG;
        if ( nd->kind == IS_NODE_OUTLET )
        {
            m = Node[j].subIndex;
            if ( Outfall[m].type == FIXED_OUTFALL ) nd->submerge = Outfall[m].fixedStage;
        }
        if ( nd->kind == IS_NODE_WETWELL )
        {
            // Fv_ReadWetWell he_input.cpp:4639, HD_WetWellData::F_ConvertUnit
            TIsWetWell* w = &IsWell[nd->wetwell];
            int found = -1;
            w->node = k;
            w->curve = -1;
            w->bottom = Node[j].invertElev;
            w->maxLevel = IsFullDepth0 ? IsFullDepth0[j] : Node[j].fullDepth;
            w->initLevel = Node[j].initDepth;
            for ( m = 0; m < NwellIn; m++ ) if ( WellIn[m].node == j ) found = m;
            if ( found < 0 ) return netError(ERR_ISEWER_WELL, Node[j].ID);
            w->minLevel = WellIn[found].minLevel / lenF;
            // Fv_ReadWetWell he_input.cpp:4782-4825 level checks
            if ( w->maxLevel < 0.0 || w->minLevel < 0.0 || w->initLevel < 0.0 ||
                 w->minLevel >= w->maxLevel || w->initLevel < w->minLevel ||
                 w->initLevel > w->maxLevel )
                return netError(ERR_ISEWER_WELL, Node[j].ID);
            if ( WellIn[found].curve[0] )
            {
                w->curve = findCurve(WellIn[found].curve);
                if ( w->curve < 0 ) return netError(ERR_NAME, WellIn[found].curve);
                w->diameter = 0.0;
            }
            else
            {
                // cylinder (Fv_ReadWetWell he_input.cpp:4836): volumes from the
                // input-unit diameter and levels, then / volume unit factor
                double d = WellIn[found].diam, area, volF;
                m = Node[j].subIndex;
                if ( Storage[m].shape != CYLINDRICAL ) return netError(ERR_ISEWER_WELL, Node[j].ID);
                if ( d <= 0.0 ) d = 2.0 * sqrt(Storage[m].a0 / PI) * UCF(LENGTH);
                volF = (UnitSystem == US) ? 1.0 : IS_M3_PER_FT3;
                area = 0.25*IS_PI_DEF*d*d;
                w->minVolume  = area*WellIn[found].minLevel / volF;
                w->maxVolume  = area*(w->maxLevel*lenF) / volF;
                w->initVolume = area*(w->initLevel*lenF) / volF;
                w->diameter = d / lenF;
            }
        }
    }

    // --- pipes and force mains
    for ( i = 0; i < L; i++ )
    {
        TIsPipe* p;
        if ( Link[i].type != CONDUIT ) continue;
        p = &IsPipe[IsLinkOfSwmm[i]];
        k = Link[i].subIndex;
        n1 = IsNodeOfSwmm[Link[i].node1];  n2 = IsNodeOfSwmm[Link[i].node2];
        p->swmm = i;
        p->kind = (Link[i].xsect.type == FORCE_MAIN) ? IS_LINK_FM : IS_LINK_PIPE;
        p->station = -1;
        p->siphon = -1;
        p->frm = n1;  p->to = n2;
        if ( IsLinkInv )
        {
            p->fromInv = IsLinkInv[2*i];
            p->toInv   = IsLinkInv[2*i+1];
        }
        else
        {
            p->fromInv = Node[Link[i].node1].invertElev + Link[i].offset1;
            p->toInv   = Node[Link[i].node2].invertElev + Link[i].offset2;
        }
        p->length  = Conduit[k].length;
        p->dia     = Link[i].xsect.yFull;
        p->n       = Conduit[k].roughness;
        p->count   = Conduit[k].barrels;
        p->kLoss   = fabs(Link[i].cLossInlet);
        p->split   = IS_SPLIT_HYD;
        p->splitPct = 0.0;
        p->infType = IS_INF_NONE;  p->infRate = 0.0;  p->infCount = 0.0;
        p->infPat  = -1;
        if ( p->kind == IS_LINK_FM )
        {
            // Fv_ReadForceMain he_input.cpp:5722: count is always 1
            p->count = 1;
            p->friction = Link[i].xsect.rBot;
            if ( IsNode[n2].kind == IS_NODE_OUTLET ) return netError(ERR_ISEWER_LINK, Link[i].ID);
            continue;
        }
        if ( IsNode[n1].kind == IS_NODE_CHAMBER ) return netError(ERR_ISEWER_LINK, Link[i].ID);

        // Fv_UpdateSlopeData he_data.cpp:5600
        p->slope = (p->fromInv - p->toInv) / p->length;
        if ( fabs(p->slope) < IS_SYS_ZERO_FLOAT ) p->slopeType = IS_SLOPE_FLAT;
        else if ( p->slope > 0.0 )                p->slopeType = IS_SLOPE_STANDARD;
        else                                      p->slopeType = IS_SLOPE_ADVERSE;
        p->adverse = p->fromInv < p->toInv;       // Fv_CheckPipeInvert

        // bottom reference levels (he_input.cpp:5256-5270)
        if ( IsNode[n1].kind != IS_NODE_OUTLET )
            IsNode[n1].bottomRef = MIN(IsNode[n1].bottomRef, p->fromInv);
        if ( IsNode[n2].kind != IS_NODE_OUTLET )
            IsNode[n2].bottomRef = MIN(IsNode[n2].bottomRef, p->toInv);
    }
    // a force main lowers the bottom reference of a loading manhole at either
    // end too (he_input.cpp:5960-5972; F_SetBottomRefLevel keeps the minimum)
    for ( i = IsNP; i < IsNPM; i++ )
    {
        TIsPipe* p = &IsPipe[i];
        if ( IsNode[p->frm].kind == IS_NODE_LOAD )
            IsNode[p->frm].bottomRef = MIN(IsNode[p->frm].bottomRef, p->fromInv);
        if ( IsNode[p->to].kind == IS_NODE_LOAD )
            IsNode[p->to].bottomRef = MIN(IsNode[p->to].bottomRef, p->toInv);
    }

    // --- station links (Fv_SetStationLinkData he_input.cpp:7577)
    for ( i = 0; i < IsNStation; i++ )
    {
        TIsPipe* p = &IsPipe[IsNPM + i];
        TIsPump* first = &IsPump[IsStation[i].pump[0]];
        p->kind = IS_LINK_STATION;
        p->swmm = -1;
        p->station = i;
        p->siphon = -1;
        p->frm = first->frm;
        p->to = first->to;
        p->count = 1;
        IsStation[i].link = IsNPM + i;
        IsStation[i].wetwell = first->wetwell;
        IsStation[i].frm = first->frm;
        IsStation[i].to = first->to;
        IsStation[i].siphon = -1;
    }

    // --- percent splits (Fv_ReadFlowSplit; F_SetSplitFlowPercent)
    for ( i = 0; i < NsplitIn; i++ )
    {
        TIsPipe* p = &IsPipe[IsLinkOfSwmm[SplitIn[i].link]];
        if ( Link[SplitIn[i].link].type != CONDUIT ) continue;
        p->splitPct = SplitIn[i].pct;
        if ( p->splitPct > 0.0 ) p->split = IS_SPLIT_PCT;
    }

    // --- infiltration (Fv_ReadInfiltrateRate; HD_PipeData::F_ConvertUnit)
    for ( i = 0; i < NinfilIn; i++ )
    {
        TIsPipe* p;
        if ( Link[InfilIn[i].link].type != CONDUIT ) continue;
        p = &IsPipe[IsLinkOfSwmm[InfilIn[i].link]];
        p->infType  = InfilIn[i].type;
        p->infRate  = InfilIn[i].rate;
        p->infCount = InfilIn[i].count;
        if ( p->infType == IS_INF_PATTERN )
        {
            p->infPat = findPattern(InfilIn[i].pat);
            if ( p->infPat < 0 ) return netError(ERR_NAME, InfilIn[i].pat);
        }
        p->infRate /= IsFlowFactor;
        if ( p->infType == IS_INF_LEN ) p->infRate *= lenF;
        else if ( p->infType == IS_INF_LENDIA || p->infType == IS_INF_SURF )
        {
            p->infRate *= lenF;
            p->infRate *= diaF;
        }
    }

    // --- loads (HD_LoadData::F_ConvertUnit) and their patterns
    for ( i = 0; i < IsNloads; i++ )
    {
        IsLoads[i].node = IsNodeOfSwmm[IsLoads[i].node];
        IsLoads[i].base = IsLoads[i].baseIn / IsFlowFactor;
        if ( IsLoads[i].pattern[0] )
        {
            IsLoads[i].pat = findPattern(IsLoads[i].pattern);
            if ( IsLoads[i].pat < 0 ) return netError(ERR_NAME, IsLoads[i].pattern);
        }
        if ( IsNode[IsLoads[i].node].kind != IS_NODE_LOAD )
            return netError(ERR_ISEWER_NODE, Node[IsNode[IsLoads[i].node].swmm].ID);
    }

    // --- pumps (Fv_ReadPump he_input.cpp:6718, HD_PumpData::F_ConvertUnit,
    //     Fv_SetPumpCoeffData he_input.cpp:12515)
    for ( i = 0; i < IsNPump; i++ )
    {
        TIsPump* pu = &IsPump[i];
        TIsPumpIn* in = NULL;
        double cap, h0, h1, q1, h2, q2;
        for ( m = 0; m < NpumpIn; m++ ) if ( PumpIn[m].link == pu->swmm ) in = &PumpIn[m];
        if ( !in ) return netError(ERR_ISEWER_PUMP, Link[pu->swmm].ID);
        pu->type = in->type;
        pu->parallel = in->extra + 1;
        cap = in->v[0];  h0 = in->v[1];  h1 = in->v[2];  q1 = in->v[3];
        h2 = in->v[4];   q2 = in->v[5];
        pu->speed = in->v[6];
        pu->pattern = -1;
        if ( in->pat[0] )
        {
            pu->pattern = findPattern(in->pat);
            if ( pu->pattern < 0 ) return netError(ERR_NAME, in->pat);
        }
        pu->staticFlow = in->qs;
        pu->baseFlow = in->qd <= 0.0 ? 1.0 : in->qd;
        if ( pu->type == IS_PUMP_CAPACITY )
        {
            h0 = h1 = q1 = h2 = q2 = 0.0;
            pu->speed = 1.0;
        }
        else
        {
            cap = 0.0;
            if ( pu->type == IS_PUMP_POINT )
            {
                h0 = 1.33334*h1;    // shut-off head
                h2 = 0.0;           // max. flow head
                q2 = 2.0*q1;        // max. flow
            }
        }
        pu->capacity = cap / IsFlowFactor;
        pu->head0 = h0 / lenF;  pu->head1 = h1 / lenF;  pu->head2 = h2 / lenF;
        pu->flow1 = q1 / IsFlowFactor;  pu->flow2 = q2 / IsFlowFactor;
        pu->a = pu->b = 0.0;
        if ( pu->type != IS_PUMP_CAPACITY )
        {
            double deltaH2 = pu->head0 - pu->head2;
            double deltaH1 = pu->head0 - pu->head1;
            pu->b = log(deltaH2/deltaH1)/log(pu->flow2/pu->flow1);
            pu->a = (pu->head0 - pu->head1)/pow((double)pu->parallel*pu->flow1, pu->b);
            if ( pu->a < 0.0 || pu->b <= 0.0 || pu->b > 20.0 )
                return netError(ERR_ISEWER_PUMP, Link[pu->swmm].ID);
        }
    }
    for ( i = 0; i < IsNStation; i++ )             // Fv_CheckParallelPumpV2
    {
        int ncap = 0;
        for ( m = 0; m < IsStation[i].npump; m++ )
            if ( IsPump[IsStation[i].pump[m]].type == IS_PUMP_CAPACITY ) ncap++;
        if ( ncap > 0 && ncap < IsStation[i].npump )
            return netError(ERR_ISEWER_PUMP, Link[IsPump[IsStation[i].pump[0]].swmm].ID);
    }

    // --- siphons (Fv_SetSiphonData he_input.cpp:6240-6717)
    if ( buildSiphons() ) return ErrorCode;

    // --- peaking curves (he_input.cpp:12269: peak-flow curves only scaled)
    IsPeakCurve = IsCoverCurve = -1;
    if ( IsPeakCurveId[0] )
    {
        IsPeakCurve = findCurve(IsPeakCurveId);
        if ( IsPeakCurve < 0 ) return netError(ERR_NAME, IsPeakCurveId);
        {
            TIsCurve* c = &IsCurves[IsPeakCurve];
            double s = 1.0 / IsFlowFactor;
            for ( i = 0; i < c->n; i++ ) { c->x[i] *= s;  c->y[i] *= s; }
        }
    }
    if ( IsCoverCurveId[0] )
    {
        IsCoverCurve = findCurve(IsCoverCurveId);
        if ( IsCoverCurve < 0 ) return netError(ERR_NAME, IsCoverCurveId);
    }

    // --- wet well level-volume curves (he_input.cpp:12218: x / length
    //     factor, y / volume factor) and their volumes (Fv_SetWetWellVolumeData)
    if ( UnitSystem != US )
    {
        int c, *done = (int*) calloc(IsNcurves + 1, sizeof(int));
        if ( !done ) return netError(ERR_MEMORY, "");
        for ( i = 0; i < IsNWell; i++ )
        {
            c = IsWell[i].curve;
            if ( c < 0 || done[c] ) continue;
            done[c] = 1;
            for ( m = 0; m < IsCurves[c].n; m++ )
            {
                IsCurves[c].x[m] *= 1.0 / lenF;
                IsCurves[c].y[m] *= 1.0 / IS_M3_PER_FT3;
            }
        }
        FREE(done);
    }
    for ( i = 0; i < IsNWell; i++ )
    {
        TIsWetWell* w = &IsWell[i];
        if ( w->curve < 0 ) continue;
        w->minVolume  = isewer_wellVolumeFromLevel(i, w->minLevel);
        w->maxVolume  = isewer_wellVolumeFromLevel(i, w->maxLevel);
        w->initVolume = isewer_wellVolumeFromLevel(i, w->initLevel);
    }

    // --- flow tolerances and minimum pipe segment (he_input.cpp:12320-12339)
    IsMinSegment /= lenF;
    IsFlowTol = IsFlowTolIn / IsFlowFactor;
    IsIterTol = IS_ITER_FLOW_TOL;
    if ( IsIterTol > IsFlowTol ) IsIterTol = IsFlowTol;

    // --- full flow of each gravity pipe x barrels (HH_PipeCapacity)
    for ( i = 0; i < IsNP; i++ ) IsPipe[i].fullFlow = isewer_fullFlowBarrel(i) * IsPipe[i].count;
    return 0;
}

static int buildSiphons(void)
//  Force mains grouped into siphons: Fv_GetSiphonSource, Fv_SetWorkSiphonData,
//  Fv_SetSiphonResource, Fv_SetSiphonElementData, Fv_SetSiphonForceMainData,
//  Fv_SetStationSiphonData (he_input.cpp:6240-7790). Force mains are indexed
//  0..IsNF-1 here (link IsNP + f).
{
    int i, j, nf = IsNF, group = -1;
    int *type, *flag, *index, *order;

    IsNSiphon = 0;
    IsSiphon = NULL;
    if ( nf <= 0 )
    {
        if ( IsNStation > 0 )                 // a station must feed a siphon
            return netError(ERR_ISEWER_PUMP, Link[IsPump[IsStation[0].pump[0]].swmm].ID);
        return 0;
    }
    type = (int*) calloc(nf, sizeof(int));
    flag = (int*) calloc(nf, sizeof(int));
    index = (int*) calloc(nf, sizeof(int));
    order = (int*) calloc(nf, sizeof(int));
    if ( !type || !flag || !index || !order )
    {
        FREE(type); FREE(flag); FREE(index); FREE(order);
        return netError(ERR_MEMORY, "");
    }
    for ( i = 0; i < nf; i++ ) { type[i] = -1;  index[i] = -1; }

    // --- split / parallel force mains are invalid (Fv_CheckSplitForceMain)
    for ( i = 0; i < nf; i++ )
        for ( j = 0; j < nf; j++ )
            if ( i != j && IsPipe[IsNP+i].frm == IsPipe[IsNP+j].frm )
            {
                FREE(type); FREE(flag); FREE(index); FREE(order);
                return netError(ERR_ISEWER_LINK, Link[IsPipe[IsNP+j].swmm].ID);
            }

    // --- leaders (Fv_GetSiphonSource)
    for ( i = 0; i < nf; i++ )
        if ( IsNode[IsPipe[IsNP+i].frm].kind == IS_NODE_LOAD ) type[i] = IS_SIPHON_STANDARD;
    for ( i = 0; i < nf; i++ )
    {
        int from = IsPipe[IsNP+i].frm, cand = 1;
        if ( IsNode[from].kind != IS_NODE_CHAMBER ) continue;
        for ( j = 0; j < nf; j++ )
            if ( i != j && IsPipe[IsNP+j].to == from ) { cand = 0;  break; }
        if ( cand ) type[i] = IS_SIPHON_STATION;
    }

    // --- chains (Fv_SetWorkSiphonData)
    for ( i = 0; i < nf; i++ )
    {
        int ord = 0, toNode;
        if ( type[i] < 0 ) continue;
        group++;
        flag[i] = 1;  index[i] = group;  order[i] = ord;
        toNode = IsPipe[IsNP+i].to;
        ord++;
        for ( j = 0; j < nf; j++ )
        {
            int from;
            if ( j == i || flag[j] ) continue;
            from = IsPipe[IsNP+j].frm;
            if ( from == toNode )
            {
                flag[j] = 1;  index[j] = group;  order[j] = ord;
                if ( IsNode[from].kind == IS_NODE_LOAD ) break;
                ord++;
                toNode = IsPipe[IsNP+j].to;
                j = -1;
            }
        }
    }
    for ( i = 0; i < nf; i++ )                 // a force main outside every siphon
        if ( index[i] < 0 )                     // (e.g. leaving a wet well) has no
        {                                       // flow or head rule in the engine
            int id = IsPipe[IsNP+i].swmm;
            FREE(type); FREE(flag); FREE(index); FREE(order);
            return netError(ERR_ISEWER_LINK, Link[id].ID);
        }
    IsNSiphon = 0;
    for ( i = 0; i < nf; i++ ) IsNSiphon = MAX(IsNSiphon, index[i] + 1);
    if ( IsNSiphon <= 0 )
    {
        FREE(type); FREE(flag); FREE(index); FREE(order);
        return netError(ERR_ISEWER_LINK, Link[IsPipe[IsNP].swmm].ID);
    }
    IsSiphon = (TIsSiphon*) calloc(IsNSiphon, sizeof(TIsSiphon));
    if ( !IsSiphon ) { FREE(type); FREE(flag); FREE(index); FREE(order); return netError(ERR_MEMORY, ""); }

    // --- elements (Fv_SetSiphonResource, Fv_SetSiphonElementData)
    {
        int siphonType = IS_SIPHON_STANDARD;    // QUIRK: carried between siphons
        for ( i = 0; i < IsNSiphon; i++ )
        {
            int cnt = 0, e, f;
            for ( f = 0; f < nf; f++ ) if ( index[f] == i ) cnt++;
            IsSiphon[i].fm = (int*) calloc(cnt + 1, sizeof(int));
            if ( !IsSiphon[i].fm ) { FREE(type); FREE(flag); FREE(index); FREE(order); return netError(ERR_MEMORY, ""); }
            IsSiphon[i].nfm = cnt;
            IsSiphon[i].station = -1;
            for ( e = 0; e < cnt; e++ )
                for ( f = 0; f < nf; f++ )
                    if ( index[f] == i && order[f] == e )
                    {
                        IsSiphon[i].fm[e] = IsNP + f;
                        if ( e == 0 && type[f] >= 0 ) siphonType = type[f];
                    }
            IsSiphon[i].type = siphonType;
        }
    }
    FREE(type); FREE(flag); FREE(index); FREE(order);

    // --- force main siphon data (Fv_SetSiphonForceMainData)
    for ( i = 0; i < IsNSiphon; i++ )
        for ( j = 0; j < IsSiphon[i].nfm; j++ )
        {
            TIsPipe* p = &IsPipe[IsSiphon[i].fm[j]];
            p->siphon = i;
            p->siphonType = IsSiphon[i].type;
        }

    // --- stations and their siphons (Fv_SetStationSiphonData he_input.cpp:7695)
    for ( i = 0; i < IsNStation; i++ )
    {
        int s = -1;
        for ( j = 0; j < IsNSiphon; j++ )
            if ( IsPipe[IsSiphon[j].fm[0]].frm == IsStation[i].to ) s = j;
        if ( s < 0 ) return netError(ERR_ISEWER_PUMP, Link[IsPump[IsStation[i].pump[0]].swmm].ID);
        IsStation[i].siphon = s;
        IsSiphon[s].station = i;
        IsSiphon[s].type = IS_SIPHON_STATION;
    }
    for ( i = 0; i < IsNSiphon; i++ )
    {
        // a chamber-led siphon no station feeds has no flow rule in the engine
        if ( IsSiphon[i].type == IS_SIPHON_STATION && IsSiphon[i].station < 0 )
            return netError(ERR_ISEWER_LINK, Link[IsPipe[IsSiphon[i].fm[0]].swmm].ID);
        for ( j = 0; j < IsSiphon[i].nfm; j++ )
            IsPipe[IsSiphon[i].fm[j]].siphonType = IsSiphon[i].type;
    }

    // --- Fv_CheckStationData he_input.cpp:7800: a station discharges into a
    //     junction chamber, and every curve pump's shut-off head must exceed
    //     the lift from the wet well's initial level to the siphon tail invert
    for ( i = 0; i < IsNStation; i++ )
    {
        TIsStation* st = &IsStation[i];
        TIsSiphon* sp = &IsSiphon[st->siphon];
        double headDown = IsPipe[sp->fm[sp->nfm - 1]].toInv;
        if ( IsNode[st->to].kind != IS_NODE_CHAMBER )
            return netError(ERR_ISEWER_PUMP, Link[IsPump[st->pump[0]].swmm].ID);
        for ( j = 0; j < st->npump; j++ )
        {
            TIsPump* pu = &IsPump[st->pump[j]];
            TIsWetWell* w = &IsWell[pu->wetwell];
            if ( pu->type == IS_PUMP_CAPACITY ) continue;
            if ( pu->head0 <= headDown - (w->initLevel + w->bottom) )
                return netError(ERR_ISEWER_LIFT, Link[pu->swmm].ID);
        }
    }
    return 0;
}

static int buildTopology(void)
//  Fv_SetLinkConnection, Fv_SetFlowSplitLoadNode, Fv_ManholeOrder,
//  Fv_LinkOrder (he_topology.cpp, he_input.cpp:5680). Topology links are the
//  gravity pipes, force mains and station links.
{
    int  i, j, k, l, m, N = IsNN, L = IsNT;
    int  *indeg, *vn, *vl, *cnt, *pool;
    int  norder, done, nvisit, nto;

    indeg = (int*) calloc(N + 1, sizeof(int));
    vn    = (int*) calloc(N + 1, sizeof(int));
    vl    = (int*) calloc(L + 1, sizeof(int));
    cnt   = (int*) calloc(N + 1, sizeof(int));
    if ( !indeg || !vn || !vl || !cnt )
    {
        FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
        return netError(ERR_MEMORY, "");
    }
    pool = IsAllocList;

    // --- connecting links in link order
    for ( i = 0; i < L; i++ ) { cnt[IsPipe[i].frm]++;  cnt[IsPipe[i].to]++; }
    for ( j = 0; j < N; j++ )
    {
        IsNode[j].conn = pool;  pool += cnt[j];
        IsNode[j].down = pool;  pool += cnt[j];
        IsNode[j].up   = pool;  pool += cnt[j];
        IsNode[j].nconn = IsNode[j].ndown = IsNode[j].nup = 0;
    }
    for ( i = 0; i < L; i++ )
    {
        TIsNode* a = &IsNode[IsPipe[i].frm];
        TIsNode* b = &IsNode[IsPipe[i].to];
        a->conn[a->nconn++] = i;
        b->conn[b->nconn++] = i;
    }

    // --- split manholes (loading manholes only) and dead ends
    for ( j = 0; j < N; j++ )
    {
        int nout = 0;
        char* id = Node[IsNode[j].swmm].ID;
        for ( m = 0; m < IsNode[j].nconn; m++ )
            if ( IsPipe[IsNode[j].conn[m]].frm == j ) nout++;
        IsNode[j].split = (IsNode[j].kind == IS_NODE_LOAD && nout > 1);
        if ( IsNode[j].kind != IS_NODE_OUTLET && nout == 0 )
        {
            FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
            return netError(ERR_ISEWER_DEADEND, id);
        }
        if ( IsNode[j].kind == IS_NODE_OUTLET && nout > 0 )
        {
            FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
            return netError(ERR_OUTFALL, id);
        }
        if ( IsNode[j].kind != IS_NODE_LOAD && nout > 1 )
        {
            FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
            return netError(ERR_ISEWER_NODE, id);
        }
    }

    // --- a force main may not leave a split manhole (he_input.cpp:5947)
    for ( i = IsNP; i < IsNPM; i++ )
        if ( IsNode[IsPipe[i].frm].split )
        {
            FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
            return netError(ERR_ISEWER_LINK, Link[IsPipe[i].swmm].ID);
        }

    // --- Fv_ManholeOrder he_topology.cpp:1377
    for ( j = 0; j < N; j++ )
        for ( m = 0; m < IsNode[j].nconn; m++ )
            if ( IsPipe[IsNode[j].conn[m]].to == j ) indeg[j]++;
    norder = 0;
    for (;;)
    {
        done = 0;
        for ( i = 0; i < N; i++ )
        {
            if ( vn[i] ) continue;
            k = i;
            while ( k >= 0 && indeg[k] == 0 )
            {
                done = 1;
                vn[k] = 1;
                IsOrder[norder++] = k;
                nvisit = 0;  nto = -1;
                for ( m = 0; m < IsNode[k].nconn; m++ )
                {
                    l = IsNode[k].conn[m];
                    if ( vl[l] ) continue;
                    nto = IsPipe[l].to;
                    if ( vn[nto] )
                    {
                        FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
                        return netError(ERR_LOOP, linkId(l));
                    }
                    vl[l] = 1;
                    if ( indeg[nto] > 0 ) { nvisit++;  indeg[nto]--; }
                }
                if ( nvisit == 1 ) k = nto;
                else break;
            }
        }
        if ( norder == N || !done ) break;
    }
    IsNorder = norder;
    if ( norder < N )
    {
        for ( j = 0; j < N; j++ ) if ( !vn[j] ) break;
        FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
        return netError(ERR_LOOP, j < N ? Node[IsNode[j].swmm].ID : "");
    }

    // --- Fv_LinkOrder he_topology.cpp:1525
    for ( i = 0; i < L; i++ ) vl[i] = 0;
    for ( m = 0; m < norder; m++ )
    {
        TIsNode* nd;
        k = IsOrder[m];
        nd = &IsNode[k];
        if ( nd->kind == IS_NODE_OUTLET ) continue;
        if ( nd->nconn == 1 )
        {
            l = nd->conn[0];
            if ( !vl[l] ) { nd->down[nd->ndown++] = l;  vl[l] = 1; }
        }
        else if ( nd->nconn > 1 )
        {
            if ( !nd->split )
            {
                int dl = -1, q;
                for ( q = 0; q < nd->nconn; q++ )
                {
                    l = nd->conn[q];
                    if ( !vl[l] && IsPipe[l].frm == k ) { dl = l; break; }
                }
                if ( dl >= 0 ) { nd->down[nd->ndown++] = dl;  vl[dl] = 1; }
                for ( q = 0; q < nd->nconn; q++ )
                    if ( nd->conn[q] != dl ) nd->up[nd->nup++] = nd->conn[q];
            }
            else
            {
                int q;
                for ( q = 0; q < nd->nconn; q++ )
                {
                    l = nd->conn[q];
                    if ( !vl[l] ) { nd->down[nd->ndown++] = l;  vl[l] = 1; }
                    else nd->up[nd->nup++] = l;
                }
            }
        }
    }

    // --- links used by the HGL sweep: conn minus down, in conn order
    for ( j = 0; j < N; j++ )
    {
        TIsNode* nd = &IsNode[j];
        int q, r, isDown;
        nd->uph = pool;  nd->nuph = 0;
        for ( q = 0; q < nd->nconn; q++ )
        {
            l = nd->conn[q];
            isDown = 0;
            for ( r = 0; r < nd->ndown; r++ ) if ( nd->down[r] == l ) isDown = 1;
            if ( !isDown ) nd->uph[nd->nuph++] = l;
        }
        pool += nd->nconn;
    }

    // --- loads of each node, in input order
    for ( j = 0; j < N; j++ ) cnt[j] = 0;
    for ( i = 0; i < IsNloads; i++ ) cnt[IsLoads[i].node]++;
    for ( j = 0; j < N; j++ ) { IsNode[j].load = pool;  pool += cnt[j];  IsNode[j].nload = 0; }
    for ( i = 0; i < IsNloads; i++ )
    {
        TIsNode* nd = &IsNode[IsLoads[i].node];
        nd->load[nd->nload++] = i;
    }

    FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
    return 0;
}

//=============================================================================
//  Steady state (he_staticrun.cpp, he_hydrun.cpp)
//=============================================================================

static double staticGetLoading(int k, double* up, double* pk, double* cv,
                               double* cc, double* inf)
//  HHU_StaticGetLoading he_hydrun.cpp:984; unpeak includes infiltration of
//  the node's downstream pipes.
{
    int m;
    double infil = 0.0;
    *up = *pk = *cv = *cc = 0.0;
    for ( m = 0; m < IsNode[k].nload; m++ )
    {
        TIsLoad* ld = &IsLoads[IsNode[k].load[m]];
        if      ( ld->type == IS_UNPEAK ) *up += IsLoadMult * ld->base;
        else if ( ld->type == IS_PEAK )   *pk += IsLoadMult * ld->base;
        else if ( ld->type == IS_COVER )
        {
            *cc += 0.001 * ld->cover;           // F_GetPerThousandUnit
            *cv += IsLoadMult * ld->base;
        }
    }
    for ( m = 0; m < IsNode[k].ndown; m++ ) infil += isewer_infiltration(IsNode[k].down[m]);
    *up += infil;
    *inf = infil;
    return 0.0;
}

static double peakFlow(double base)
//  HHU_StaticGetPeakFlow he_hydrun.cpp:1423
{
    if ( IsPeakCurve < 0 ) return IsPeakK * pow(base, IsPeakP);
    return curveY(IsPeakCurve, base);
}

static double coverFlow(double cflow, double cunit)
//  HHU_StaticGetCoverFlow he_hydrun.cpp:1458
{
    double varC, factor;
    if ( cunit < IS_ZERO_TOL ) return 0.0;
    if ( IsCoverCurve < 0 )
    {
        varC = pow(cunit * IsPeakE, IsPeakC);
        factor = (IsPeakA / (IsPeakB + varC)) + IsPeakD;
        return factor * cflow;
    }
    return cflow * curveY(IsCoverCurve, cunit * 1000);
}

static double staticTotal(int k)
//  HHU_StaticGetTotalLoading he_hydrun.cpp:933
{
    double up, pk, cv, cc, inf;
    staticGetLoading(k, &up, &pk, &cv, &cc, &inf);
    return up + peakFlow(pk) + coverFlow(cv, cc) + 0.0;
}

//-----------------------------------------------------------------------------
//  Force mains, siphons and pump stations (steady state)
//-----------------------------------------------------------------------------

static double fmHeadloss(int l, double flow)
//  HD_ForceMainData::F_GetHeadloss he_data.cpp:6339
{
    TIsPipe* p = &IsPipe[l];
    double headloss = 0.0, varN, varD;
    if ( p->siphonType == IS_SIPHON_STATION || p->friction >= IS_MAX_MANNING )
    {
        double frictionPower = pow(p->friction, 1.852);
        double diameterPower = pow(p->dia, 4.871);
        double flowPower     = pow(flow, 1.852);
        varN = 4.727*p->length*flowPower;
        varD = frictionPower*diameterPower;
        if ( fabs(varD) > IS_SYS_ZERO ) headloss = varN/varD;
    }
    else
    {
        varN = 4.66*p->friction*p->friction*p->length*flow*flow;
        varD = pow(p->dia, 5.33);
        if ( fabs(varD) > IS_SYS_ZERO ) headloss = varN/varD;
    }
    return headloss;
}

static void fmHydraulics(int l, double flow, double* velocity, double* headloss)
//  HD_ForceMainData::F_GetHydraulicData he_data.cpp:6443
{
    TIsPipe* p = &IsPipe[l];
    double area;
    *velocity = 0.0;
    *headloss = 0.0;
    if ( p->length <= 0.0 || p->dia <= 0.0 || p->friction <= 0.0 || p->siphon < 0 ) return;
    area = 0.25*IS_PI_DEF*p->dia*p->dia;
    if ( fabs(area) > IS_SYS_ZERO ) *velocity = flow/area;
    *headloss = fmHeadloss(l, flow);
}

static double siphonHeadloss(int s, double q)
//  HD_SiphonData::F_GetSiphonHeadloss he_data.cpp:8114
{
    int i;
    double h = 0.0;
    for ( i = 0; i < IsSiphon[s].nfm; i++ ) h += fmHeadloss(IsSiphon[s].fm[i], q);
    return h;
}

static double siphonTailHead(int s)
//  HD_SiphonData::F_GetSiphonTailHead he_data.cpp:8082
{
    TIsPipe* p = &IsPipe[IsSiphon[s].fm[IsSiphon[s].nfm - 1]];
    return p->toInv + p->dia;
}

static double siphonNewtonFactor(int s, double q)
//  HD_SiphonData::F_GetSiphonNewtonFactor he_data.cpp:8180
{
    int i;
    double factor = 0.0, constant = 4.727*1.852, flowPower = pow(q, 0.852);
    for ( i = 0; i < IsSiphon[s].nfm; i++ )
    {
        TIsPipe* p = &IsPipe[IsSiphon[s].fm[i]];
        double frictionPower = pow(p->friction, 1.852);
        double diameterPower = pow(p->dia, 4.871);
        double varN = constant*p->length*flowPower;
        double varD = frictionPower*diameterPower;
        if ( fabs(varD) > IS_SYS_ZERO ) factor += varN/varD;
    }
    return factor;
}


// --- force main and siphon helpers shared with isewer_eps.c
double isewer_siphonHeadloss(int s, double q)   { return siphonHeadloss(s, q); }
double isewer_siphonTailHead(int s)             { return siphonTailHead(s); }
double isewer_siphonNewtonFactor(int s, double q) { return siphonNewtonFactor(s, q); }
void   isewer_fmHydraulics(int l, double q, double* v, double* h) { fmHydraulics(l, q, v, h); }

static int    *SVisit;                    // links visited by the flow pass
static double *SPumpHead0, *SPumpA, *SPumpB, *SPumpPB, *SPumpOn, *SPumpSpeed;
static double  SUpdatedPumpFlow;          // hm_updatedPumpFlow

static int isPointPump(int p) { return IsPump[p].type != IS_PUMP_CAPACITY; }

static void standardSiphonFlow(int s, double inflow)
//  Fc_StaticGetStandardSiphonFlow he_hydrun.cpp:4273
{
    int i;
    for ( i = 0; i < IsSiphon[s].nfm; i++ )
    {
        int l = IsSiphon[s].fm[i];
        SQ[l] = inflow;
        fmHydraulics(l, inflow, &SVel[l], &SFmLoss[l - IsNP]);
    }
}

static void stationSiphonFlow(int s)
//  Fc_StaticGetStationSiphonFlow he_hydrun.cpp:4313. QUIRK: one unit of each
//  pump (no parallel count), read from the pump element flows.
{
    TIsStation* st = &IsStation[IsSiphon[s].station];
    double stationFlow = 0.0;
    int i;
    for ( i = 0; i < st->npump; i++ ) stationFlow += SQ[IsNPM + st->pump[i]];
    for ( i = 0; i < IsSiphon[s].nfm; i++ )
    {
        int l = IsSiphon[s].fm[i];
        SQ[l] = stationFlow;
        SUP[l] = stationFlow;
        fmHydraulics(l, stationFlow, &SVel[l], &SFmLoss[l - IsNP]);
    }
}

static void siphonFlow(int l, double inflow, double up, double pk, double cv,
                       double cc, double inf)
//  Fc_StaticGetSiphonFlow he_hydrun.cpp:4127 (siphon flow option Upstream)
{
    int s, i;
    if ( SVisit[l] ) return;
    s = IsPipe[l].siphon;
    if ( IsSiphon[s].type == IS_SIPHON_STANDARD ) standardSiphonFlow(s, inflow);
    else stationSiphonFlow(s);
    for ( i = 0; i < IsSiphon[s].nfm; i++ )
    {
        int f = IsSiphon[s].fm[i];
        if ( IsSiphon[s].type == IS_SIPHON_STATION )
        {
            if ( IsPumpContinuity == 0 )
            {
                SQ[f] = inflow;
                SUP[f] = up;  SPK[f] = pk;  SCV[f] = cv;  SCC[f] = cc;  SINF[f] = inf;
            }
            else
            {
                SUP[f] = SQ[f];
                SPK[f] = SCV[f] = SCC[f] = SINF[f] = 0.0;
            }
        }
        else
        {
            SUP[f] = up;  SPK[f] = pk;  SCV[f] = cv;  SCC[f] = cc;  SINF[f] = inf;
        }
    }
    for ( i = 0; i < IsSiphon[s].nfm; i++ ) SVisit[IsSiphon[s].fm[i]] = 1;
}

static void topSiphonLink(int l, double inflow, double up, double pk, double cv,
                          double cc, double inf)
//  HHU_StaticTopSiphonLink he_hydrun.cpp:863 (force main leaving a top manhole)
{
    int s, i;
    if ( SVisit[l] ) return;
    s = IsPipe[l].siphon;
    if ( s < 0 || IsSiphon[s].type != IS_SIPHON_STANDARD ) return;
    standardSiphonFlow(s, inflow);
    for ( i = 0; i < IsSiphon[s].nfm; i++ )
    {
        int f = IsSiphon[s].fm[i];
        SUP[f] = up;  SPK[f] = pk;  SCV[f] = cv;  SCC[f] = cc;  SINF[f] = inf;
    }
    for ( i = 0; i < IsSiphon[s].nfm; i++ ) SVisit[IsSiphon[s].fm[i]] = 1;
}

static double pumpFlowBisection(int pumpIndex, double siphonQ)
//  Fc_PumpFlowBisection he_hydrun.cpp:5039
{
    TIsPump* pu = &IsPump[pumpIndex];
    TIsStation* st = &IsStation[pu->station];
    TIsWetWell* w = &IsWell[pu->wetwell];
    int s = st->siphon, i;
    double wellHead = w->initLevel + w->bottom;
    double tailHead = siphonTailHead(s);
    double loss, functX, functX1, functX11, speedCoef, adjA;
    for (;;)
    {
        loss = siphonHeadloss(s, siphonQ);
        functX = 0.0;  functX11 = 0.0;
        for ( i = 0; i < st->npump; i++ )
        {
            functX1 = 0.0;
            if ( SPumpOn[i] > 0 && isPointPump(st->pump[i]) )
            {
                functX = wellHead - tailHead - loss + pow(SPumpSpeed[i], 2.0)*SPumpHead0[i];
                if ( functX < 0.0 ) break;
                speedCoef = pow(SPumpSpeed[i], SPumpB[i] - 2.0);
                adjA = SPumpA[i]/speedCoef;
                functX1 = pow(functX/adjA, 1.0/SPumpB[i]);
            }
            functX11 += functX1;
        }
        if ( functX < 0.0 ) { siphonQ = 0.99*siphonQ;  continue; }
        break;
    }
    SUpdatedPumpFlow = siphonQ;
    return functX11 - siphonQ;
}

static int pointPumpFlow(int pumpIndex, double* pumpFlow)
//  Fc_StaticGetPointPumpFlowV2 he_hydrun.cpp:4698: Newton on the station
//  flow, bisection if Newton fails. Returns 1 on a fatal failure.
{
    TIsPump* pu0 = &IsPump[pumpIndex];
    TIsStation* st = &IsStation[pu0->station];
    TIsWetWell* w = &IsWell[pu0->wetwell];
    int s = st->siphon, n = st->npump, i, iter = 0, maxIter = 150;
    int parallelCount = pu0->parallel;
    double wellHead = w->initLevel + w->bottom;
    double tailHead = siphonTailHead(s);
    double totalPumpHead = 0.0, limiting, sQ, loss = 0.0, newQ = 0.0, deltaFlow = 1000.0;
    double speedCoef = 0.0, adjA = 0.0, prevFunctX = 0.0, prevQ = 0.0, tol = 0.0001;
    double Q = 0.0, functX;

    for ( i = 0; i < n; i++ )
    {
        int p = st->pump[i];
        if ( SPumpOn[i] > 0 )
        {
            if ( isPointPump(p) ) Q += parallelCount*IsPump[p].flow1;
            if ( SPumpSpeed[i] < 0.0 ) SPumpSpeed[i] = 1.0;
        }
    }
    for ( i = 0; i < n; i++ ) totalPumpHead += pow(SPumpSpeed[i], 2.0)*SPumpHead0[i];
    limiting = totalPumpHead - (tailHead - wellHead);
    if ( limiting <= 0.0 )
    {
        for ( i = 0; i < n; i++ )
            if ( SPumpOn[i] > 0 && isPointPump(st->pump[i]) )
            {
                SQ[IsNPM + st->pump[i]] = 0.0;
                SPumpHead[st->pump[i]] = 0.0;
            }
        *pumpFlow = Q;
        return 0;
    }
    sQ = Q;
    while ( deltaFlow > tol && iter <= maxIter )
    {
        double functX11 = 0.0, functX12, functX22 = 0.0, functX23;
        iter++;
        for (;;)                                    // tryAgain
        {
            loss = siphonHeadloss(s, sQ);
            if ( loss > limiting ) { sQ = sQ/2.0;  Q = sQ;  continue; }
            break;
        }
        functX = 0.0;
        for ( i = 0; i < n; i++ )
        {
            double functX1 = 0.0, functX21 = 0.0, functP, functX2;
            int p = st->pump[i];
            if ( SPumpOn[i] > 0 && isPointPump(p) )
            {
                functX = wellHead - tailHead - loss + pow(SPumpSpeed[i], 2.0)*SPumpHead0[i];
                if ( functX < 0.0 ) break;
                speedCoef = pow(SPumpSpeed[i], SPumpB[i] - 2.0);
                adjA = SPumpA[i]/speedCoef;
                functX1 = pow(functX/adjA, 1.0/SPumpB[i]);
                functP = pow(functX/adjA, SPumpPB[i]);
                functX2 = (1.0/SPumpB[i])*functP;
                functX21 = functX2*(siphonNewtonFactor(s, sQ)/adjA);
            }
            functX11 += functX1;
            functX22 += functX21;
        }
        if ( functX >= 0.0 )
        {
            functX12 = Q - functX11;
            functX23 = 1 + functX22;
            if ( fabs(functX23) <= IsFlowTol ) newQ = Q;
            else newQ = Q - functX12/functX23;
            if ( newQ < 0.001 )
            {
                for ( i = 0; i < n; i++ )
                    if ( SPumpOn[i] > 0 && isPointPump(st->pump[i]) )
                    {
                        SQ[IsNPM + st->pump[i]] = 0.0;
                        SPumpHead[st->pump[i]] = 0.0;
                    }
                *pumpFlow = Q;
                return 0;
            }
            prevFunctX = functX;
            deltaFlow = fabs((newQ - Q)/newQ);
            prevQ = Q;
            Q = newQ;
            sQ = Q;
        }
        else
        {
            if ( Q <= 0.001 )
            {
                for ( i = 0; i < n; i++ )
                    if ( SPumpOn[i] > 0 && isPointPump(st->pump[i]) )
                    {
                        SQ[IsNPM + st->pump[i]] = 0.0;
                        SPumpHead[st->pump[i]] = 0.0;
                    }
                *pumpFlow = Q;
                return 0;
            }
            if ( prevFunctX >= 0.0 )
            {
                Q = 0.5*(Q + prevQ);
                prevQ = sQ;
                deltaFlow = fabs((Q - prevQ)/Q);
                sQ = Q;
            }
            else
            {
                Q = 0.5*Q;
                prevQ = sQ;
                deltaFlow = fabs((Q - prevQ)/Q);
                sQ = Q;
            }
            prevFunctX = functX;
        }
    }
    if ( deltaFlow > tol )                          // bisection
    {
        int bIter = 0;
#ifdef ISEWER_TRACE
        fprintf(stderr, "ISEWER_TRACE bisection pump %s\n", Link[pu0->swmm].ID);
#endif
        double lowerQ = IsFlowTol, upperQ = pu0->flow2, lowerF, upperF, avgQ, avgF;
        lowerF = pumpFlowBisection(pumpIndex, lowerQ);  lowerQ = SUpdatedPumpFlow;
        upperF = pumpFlowBisection(pumpIndex, upperQ);  upperQ = SUpdatedPumpFlow;
        avgQ = 0.5*(lowerQ + upperQ);
        avgF = pumpFlowBisection(pumpIndex, avgQ);      avgQ = SUpdatedPumpFlow;
        do
        {
            if ( avgF*upperF <= 0.0 )
            {
                lowerQ = avgQ;
                lowerF = pumpFlowBisection(pumpIndex, lowerQ);  lowerQ = SUpdatedPumpFlow;
            }
            else if ( avgF*lowerF <= 0.0 )
            {
                upperQ = avgQ;
                upperF = pumpFlowBisection(pumpIndex, upperQ);  upperQ = SUpdatedPumpFlow;
            }
            else
            {
                upperQ = 1.01*upperQ;
                upperF = pumpFlowBisection(pumpIndex, upperQ);  upperQ = SUpdatedPumpFlow;
            }
            avgQ = 0.5*(lowerQ + upperQ);
            avgF = pumpFlowBisection(pumpIndex, avgQ);  avgQ = SUpdatedPumpFlow;
            bIter++;
        } while ( (fabs(avgF)/avgQ) >= 0.00001 && bIter < 200 );
        if ( bIter >= 100 ) return 1;
        loss = siphonHeadloss(s, avgQ);
    }
    functX = 0.0;
    for ( i = 0; i < n; i++ )                       // Compute Output Data
    {
        int p = st->pump[i];
        if ( SPumpOn[i] > 0 && isPointPump(p) )
        {
            int e = IsNPM + p;
            functX = wellHead - tailHead - loss + pow(SPumpSpeed[i], 2.0)*SPumpHead0[i];
            speedCoef = pow(SPumpSpeed[i], SPumpB[i] - 2.0);
            adjA = SPumpA[i]/speedCoef;
            SQ[e] = pow(functX/adjA, 1.0/SPumpB[i]);
            SPumpHead[p] = pow(SPumpSpeed[i], 2.0)*SPumpHead0[i] - adjA*pow(fabs(SQ[e]), SPumpB[i]);
            SPumpHead[p] = MAX(0.0, SPumpHead[p]);
        }
    }
    *pumpFlow = Q;
    return 0;
}

static int stationFlow(int link, double up, double pk, double cv, double cc,
                       double inf)
//  Fc_StaticGetStationFlowV2 he_hydrun.cpp:4365. Returns 1 on a fatal failure.
{
    TIsStation* st = &IsStation[IsPipe[link].station];
    int n = st->npump, i, numOff = 0;
    TIsWetWell* w = &IsWell[st->wetwell];
    double totalCap = 0.0, totalInflow, wellHead, bottomElev;

    for ( i = 0; i < n; i++ )
    {
        SPumpHead0[i] = SPumpA[i] = SPumpB[i] = SPumpPB[i] = 0.0;
        SPumpOn[i] = 1;  SPumpSpeed[i] = 1.0;
    }
    totalInflow = up + pk + cv + inf + 0.0;
    wellHead = w->initLevel + w->bottom;
    bottomElev = w->bottom;
    if ( wellHead <= bottomElev && totalInflow <= IsFlowTol )
        for ( i = 0; i < n; i++ ) SPumpOn[i] = 0;
    for ( i = 0; i < n; i++ )                      // pumps off by usage or speed
    {
        int p = st->pump[i], e = IsNPM + p;
        double usage = 1.0;
        SPumpSpeed[i] = IsPump[p].speed;
        if ( SPumpSpeed[i] < IS_SYS_ZERO ) usage = 0.0;
        if ( usage < IS_SYS_ZERO )
        {
            SQ[e] = SUP[e] = SPK[e] = SCV[e] = SCC[e] = SINF[e] = 0.0;
            SPumpHead[p] = 0.0;
            SPumpOn[i] = 0;
            numOff += 1;
        }
    }
    if ( numOff < n )
    {
        for ( i = 0; i < n; i++ )
        {
            int p = st->pump[i], e = IsNPM + p;
            if ( SPumpOn[i] > 0 && IsPump[p].type == IS_PUMP_CAPACITY )
            {
                double totalFlow;
                SQ[e] = IsPump[p].capacity;
                totalCap += IsPump[p].parallel*SQ[e];
                totalFlow = IsPump[p].parallel*SQ[e];
                SPumpHead[p] = siphonTailHead(st->siphon) + siphonHeadloss(st->siphon, totalFlow) - wellHead;
                SPumpHead[p] = MAX(0.0, SPumpHead[p]);
            }
        }
        if ( isPointPump(st->pump[n - 1]) )
        {
            double pumpQ = 0.0;
            for ( i = 0; i < n; i++ )
            {
                int p = st->pump[i];
                if ( SPumpOn[i] > 0 )
                {
                    SPumpHead0[i] = IsPump[p].head0;
                    SPumpA[i] = IsPump[p].a;
                    SPumpB[i] = IsPump[p].b;
                    SPumpPB[i] = (1.0/SPumpB[i]) - 1.0;
                }
            }
            if ( pointPumpFlow(st->pump[0], &pumpQ) ) return 1;
            totalCap = pumpQ;
        }
    }
    if ( IsPumpContinuity == 0 )
    {
        SUP[link] = up;  SPK[link] = pk;  SCV[link] = cv;  SCC[link] = cc;  SINF[link] = inf;
    }
    else
    {
        SUP[link] = totalCap;
        SPK[link] = SCV[link] = SCC[link] = SINF[link] = 0.0;
    }
    SVisit[link] = 1;
    return 0;
}

static void steadySolve(void)
//  Fv_StaticRun he_staticrun.cpp:344, then Fv_GetPipeHydraulic,
//  Fv_GetCriticalDepth, Fv_GetPipeFlowDepth, Fv_StaticFilterResult, He_RunHead.
{
    int i, k, m, N = IsNN, LE = IsNLE, maxPumps = 1;
    int *visit;
    double *Q;

    for ( i = 0; i < IsNStation; i++ ) maxPumps = MAX(maxPumps, IsStation[i].npump);
    SQ  = (double*) calloc(LE + 1, sizeof(double));
    SUP = (double*) calloc(LE + 1, sizeof(double));
    SPK = (double*) calloc(LE + 1, sizeof(double));
    SCV = (double*) calloc(LE + 1, sizeof(double));
    SCC = (double*) calloc(LE + 1, sizeof(double));
    SINF = (double*) calloc(LE + 1, sizeof(double));
    SDepth = (double*) calloc(LE + 1, sizeof(double));
    SVel = (double*) calloc(LE + 1, sizeof(double));
    SCrit = (double*) calloc(LE + 1, sizeof(double));
    SFull = (int*) calloc(LE + 1, sizeof(int));
    SFtype = (int*) calloc(LE + 1, sizeof(int));
    SPartial = (int*) calloc(LE + 1, sizeof(int));
    SAdjFull = (int*) calloc(LE + 1, sizeof(int));
    SStatic = (double*) calloc(N + 1, sizeof(double));
    SHead = (double*) calloc(N + 1, sizeof(double));
    SInc = (double*) calloc(N + 1, sizeof(double));
    SFmLoss = (double*) calloc(IsNF + 1, sizeof(double));
    SPumpHead = (double*) calloc(IsNPump + 1, sizeof(double));
    SWellHead = (double*) calloc(IsNWell + 1, sizeof(double));
    SPumpHead0 = (double*) calloc(maxPumps, sizeof(double));
    SPumpA = (double*) calloc(maxPumps, sizeof(double));
    SPumpB = (double*) calloc(maxPumps, sizeof(double));
    SPumpPB = (double*) calloc(maxPumps, sizeof(double));
    SPumpOn = (double*) calloc(maxPumps, sizeof(double));
    SPumpSpeed = (double*) calloc(maxPumps, sizeof(double));
    visit = (int*) calloc(LE + 1, sizeof(int));
    if ( !SQ || !SUP || !SPK || !SCV || !SCC || !SINF || !SDepth || !SVel ||
         !SCrit || !SFull || !SFtype || !SPartial || !SAdjFull || !SStatic ||
         !SHead || !SInc || !SFmLoss || !SPumpHead || !SWellHead || !SPumpHead0 ||
         !SPumpA || !SPumpB || !SPumpPB || !SPumpOn || !SPumpSpeed || !visit )
    {
        FREE(visit);
        report_writeErrorMsg(ERR_MEMORY, "");
        return;
    }
    Q = SQ;
    SVisit = visit;

    // --- flows, upstream to downstream
    for ( m = 0; m < IsNorder; m++ )
    {
        TIsNode* nd;
        double up, pk, cv, cc, inf, inflow;
        k = IsOrder[m];
        nd = &IsNode[k];
        if ( nd->nconn == 1 && nd->kind == IS_NODE_LOAD )
        {
            // Fv_GetSingleSourceNodeFlow he_staticrun.cpp:423
            int l = nd->conn[0];
            SStatic[k] = staticTotal(k);
            staticGetLoading(k, &SUP[l], &SPK[l], &SCV[l], &SCC[l], &SINF[l]);
            Q[l] = SUP[l] + peakFlow(SPK[l]) + coverFlow(SCV[l], SCC[l]) + 0.0;
            if ( IsPipe[l].kind == IS_LINK_FM )
                topSiphonLink(l, Q[l], SUP[l], SPK[l], SCV[l], SCC[l], SINF[l]);
            visit[l] = 1;
            continue;
        }

        // Fv_GetConnectNodeFlow he_staticrun.cpp:518
        up = pk = cv = cc = inf = 0.0;
        if ( nd->kind == IS_NODE_LOAD )
        {
            SStatic[k] = staticTotal(k);
            staticGetLoading(k, &up, &pk, &cv, &cc, &inf);
        }
        for ( i = 0; i < nd->nconn; i++ )         // HHU_StaticAddPrecedeFlow
        {
            int l = nd->conn[i];
            if ( visit[l] )
            {
                up += SUP[l];  pk += SPK[l];  cv += SCV[l];  cc += SCC[l];  inf += SINF[l];
            }
        }
        inflow = up + peakFlow(pk) + coverFlow(cv, cc) + 0.0;
        if ( fabs(inflow) < IsFlowTol && nd->kind != IS_NODE_WETWELL &&
             nd->kind != IS_NODE_CHAMBER )
        {
            // QUIRK: the loop starts at the second downstream link (he_staticrun.cpp:609)
            for ( i = 1; i < nd->ndown; i++ )
            {
                int l = nd->down[i];
                Q[l] = SUP[l] = SPK[l] = SCV[l] = SINF[l] = 0.0;
                visit[l] = 1;
            }
            continue;
        }
        if ( nd->ndown == 1 )                     // HHU_StaticSingleDownLink
        {
            int l = nd->down[0];
            if ( IsPipe[l].kind == IS_LINK_PIPE )
            {
                Q[l] = inflow;  SUP[l] = up;  SPK[l] = pk;  SCV[l] = cv;  SCC[l] = cc;
                SINF[l] = inf;
                visit[l] = 1;
            }
            else if ( IsPipe[l].kind == IS_LINK_FM )
                siphonFlow(l, inflow, up, pk, cv, cc, inf);
            else if ( stationFlow(l, up, pk, cv, cc, inf) )
            {
                report_writeErrorMsg(ERR_ISEWER_PUMPFLOW, linkId(l));
                FREE(visit);
                return;
            }
        }
        else if ( nd->ndown > 1 )
        {
            splitLinks(k, inflow, Q, visit);       // HHU_StaticFlowSplitLinks
            for ( i = 0; i < nd->ndown; i++ )
            {
                int l = nd->down[i];
                double ratio = Q[l] / inflow;
                SUP[l] = ratio * up;  SPK[l] = ratio * pk;  SCV[l] = ratio * cv;
                SCC[l] = ratio * cc;  SINF[l] = ratio * inf;
                visit[l] = 1;
            }
        }
    }
    FREE(visit);
    SVisit = NULL;

    // --- gravity pipe hydraulics (he_staticrun.cpp:681-904, gravity pipes only)
    for ( i = 0; i < IsNP; i++ )                   // Fv_GetPipeHydraulic
        SFull[i] = (Q[i] > IsFlowTol && Q[i] >= IsPipe[i].fullFlow);
    for ( i = 0; i < IsNP; i++ )                   // Fv_GetCriticalDepth
    {
        if ( fabs(Q[i]) < IsFlowTol ) SCrit[i] = 0.0;
        else
        {
            double q = MIN(Q[i], IsPipe[i].fullFlow) / IsPipe[i].count;
            SCrit[i] = isewer_criticalDepth(i, q, IsPipe[i].dia);
        }
    }
    for ( i = 0; i < IsNP; i++ )                   // Fv_GetPipeFlowDepth
    {
        double q = Q[i] / IsPipe[i].count;
        SDepth[i] = SFull[i] ? IsPipe[i].dia : isewer_flowDepth(i, q);
        if ( SDepth[i] > 0.0 )
            SVel[i] = isewer_pipeVelocity(SDepth[i], 0.5 * IsPipe[i].dia,
                                          MIN(q, Q[i] / IsPipe[i].count));
        else SVel[i] = 0.0;
    }
    for ( i = 0; i < IsNP; i++ )                   // Fv_StaticFilterResult
    {
        if ( !(SDepth[i] / IsPipe[i].dia >= IS_DD_TOL) )
        {
            SDepth[i] = 0.0;
            SVel[i] = 0.0;
        }
    }

    // --- heads (Fv_StaticGetHead he_staticrun.cpp:905): wet wells at their
    //     initial level
    for ( i = 0; i < IsNWell; i++ ) SWellHead[i] = IsWell[i].initLevel + IsWell[i].bottom;
    isewer_setHeadPressureData(SWellHead, SPumpHead);
    if ( isewer_runHead(Q, SDepth, SHead, SFtype, SPartial, SAdjFull, 1) && !ErrorCode )
        report_writeErrorMsg(ERR_ISEWER_LINK, "(head sweep)");
}

static void splitLinks(int k, double inflow, double* Q, int* visit)
//  HHU_StaticFlowSplitLinks he_hydrun.cpp:477 + Fc_GetDownstreamSplitData,
//  HHU_CompPercentSplitFlow he_hydrun.cpp:1744 (percent type)
{
    TIsNode* nd = &IsNode[k];
    int *hyd, *pct, nhyd = 0, npct = 0, m;
    double s = 0.0, tot = 0.0, ratio;

    hyd = (int*) calloc(nd->ndown + 1, sizeof(int));
    pct = (int*) calloc(nd->ndown + 1, sizeof(int));
    if ( !hyd || !pct ) { FREE(hyd); FREE(pct); report_writeErrorMsg(ERR_MEMORY, ""); return; }
    for ( m = 0; m < nd->ndown; m++ )
    {
        int l = nd->down[m];
        if ( visit[l] ) continue;
        Q[l] = 0.0;
        if ( IsPipe[l].split == IS_SPLIT_HYD ) hyd[nhyd++] = l;
        else pct[npct++] = l;
    }
    if ( npct > 0 )
    {
        int scale;
        for ( m = 0; m < npct; m++ ) tot += IsPipe[pct[m]].splitPct;
        ratio = 1.0;
        if ( tot == 0.0 ) ratio = 0.0;
        else if ( tot < 100.0 ) ratio = (nhyd > 0) ? 1.0 : 100.0 / tot;
        else if ( tot > 100.0 ) ratio = 100.0 / tot;
        scale = !(tot < 100.0 && nhyd > 0);
        for ( m = 0; m < npct; m++ )
        {
            int l = pct[m];
            double w = IsPipe[l].splitPct;
            if ( scale ) w = w * ratio;
            Q[l] = 0.01 * w * inflow;
            s += Q[l];
            visit[l] = 1;
        }
    }
    if ( nhyd == 1 )
    {
        Q[hyd[0]] = inflow - s;
        visit[hyd[0]] = 1;
    }
    else if ( nhyd > 1 ) multiHydSplit(hyd, nhyd, inflow - s, Q, visit);
    FREE(hyd);  FREE(pct);
}

static void multiHydSplit(int* hyd, int nhyd, double qh, double* Q, int* visit)
//  Fc_StaticMultiHydSplitFlow he_hydrun.cpp:690: a common water level found
//  by bisection.  QUIRK: the stop test uses the previous level's flows.
{
    int m, allFull;
    double minInv, maxCrown, wl, out, lo, hi, d, D, tf;
    TIsPipe* p;
    long iter = 0;

    if ( qh < IsFlowTol )
    {
        for ( m = 0; m < nhyd; m++ ) { Q[hyd[m]] = 0.0;  visit[hyd[m]] = 1; }
        return;
    }
    p = &IsPipe[hyd[0]];
    minInv = p->fromInv;  maxCrown = minInv + p->dia;
    for ( m = 1; m < nhyd; m++ )
    {
        p = &IsPipe[hyd[m]];
        minInv = MIN(minInv, p->fromInv);
        maxCrown = MAX(maxCrown, p->fromInv + p->dia);
    }
    wl = 0.5 * (minInv + maxCrown);
    out = 0.0;
    for ( m = 0; m < nhyd; m++ )
    {
        int l = hyd[m];
        p = &IsPipe[l];
        d = wl - p->fromInv;
        if ( d <= 0.0 ) Q[l] = 0.0;
        else
        {
            d = MIN(d, p->dia);
            Q[l] = p->count * isewer_linkFlow(l, d, -1.0, -1.0);
        }
        out += Q[l];
    }
    lo = minInv;  hi = maxCrown;
    if ( out > qh ) hi = wl;
    if ( out < qh ) lo = wl;
    wl = 0.5 * (hi + lo);
    while ( fabs(out - qh) / out > IsIterTol )
    {
        out = 0.0;
        allFull = 1;
        for ( m = 0; m < nhyd; m++ )
        {
            int l = hyd[m];
            p = &IsPipe[l];
            D = p->dia;
            d = wl - p->fromInv;
            if ( d <= 0.0 ) { Q[l] = 0.0;  allFull = 0; }
            else
            {
                Q[l] = p->count * isewer_linkFlow(l, MIN(d, D), -1.0, -1.0);
                if ( d < (D - IS_ZERO_TOL) ) allFull = 0;
            }
            out += Q[l];
        }
        if ( allFull )
        {
            tf = 0.0;
            for ( m = 0; m < nhyd; m++ ) tf += IsPipe[hyd[m]].fullFlow;
            if ( tf <= 0.0 ) break;
            for ( m = 0; m < nhyd; m++ ) Q[hyd[m]] = qh * IsPipe[hyd[m]].fullFlow / tf;
            break;
        }
        if ( out > qh ) hi = wl;
        if ( out < qh ) lo = wl;
        wl = 0.5 * (hi + lo);
        if ( ++iter > 100000 ) break;           // safety only
    }
    for ( m = 0; m < nhyd; m++ ) visit[hyd[m]] = 1;
}

//=============================================================================
//  HGL sweep (He_RunHead, he_runhead.cpp) - shared by steady state and EPS
//=============================================================================

static const double* HQ;      // link flows (cfs), link element indexes
static const double* HY;      // gravity pipe depths (ft)
static const double* HWell;   // wet well heads (ft)
static const double* HPumpH;  // pump head increases (ft)
static double* HHyd;          // working node heads
static double* HOut;          // output node heads
static int*    HLevel;        // upstream links done per node
static int*    HFtype;
static int*    HPartial;
static int*    HAdjFull;
static int*    HSplitOrder;   // order in which split nodes were reached
static int*    HSplitNode;    // node for each split order
static int*    HFmFlag;       // force main heads done (mv_pForceMainFlag)
static int     HCounter;
static int     HError;

void isewer_setHeadPressureData(const double* wellHead, const double* pumpHead)
//  F_SetWetWellHead, F_SetPumpHead he_runhead.cpp:270, 483
{
    HWell = wellHead;
    HPumpH = pumpHead;
}

static int fully(int k) { return HLevel[k] == IsNode[k].nuph; }

static void minorLoss(int k, double* h)
//  Fv_AdjustNodeMinorLoss he_runhead.cpp:2258: full-pipe velocity of the
//  total flow of each flowing downstream link; K from the SWMM inlet loss.
//  Force mains use Q/(3.14 r^2).
{
    TIsNode* nd = &IsNode[k];
    int m;
    for ( m = 0; m < nd->ndown; m++ )
    {
        int l = nd->down[m];
        if ( HQ[l] > IS_FLOW_TOL )
        {
            double D = IsPipe[l].dia, v;
            if ( IsPipe[l].kind == IS_LINK_PIPE )
                v = isewer_pipeVelocity(D, 0.5 * D, HQ[l]);
            else if ( IsPipe[l].kind == IS_LINK_FM )
            {
                double r = 0.5 * D;
                v = HQ[l] / (3.14 * r * r);
            }
            else continue;
            *h += fabs(IsPipe[l].kLoss) * ((v * v) / (2.0 * 32.174));
        }
    }
}

static void adjustNode(int k, double h)
//  Fv_AdjustNodeHead he_runhead.cpp:1910
{
    TIsNode* nd = &IsNode[k];
    int m;
    double dh = -IS_BIG;
    switch ( nd->kind )
    {
      case IS_NODE_LOAD:
        // Fv_AdjustLoadNodeDownstream he_runhead.cpp:2152
        // QUIRK: h = max(h, dh) runs for every link, with dh kept from the
        // last flowing one.
        for ( m = 0; m < nd->ndown; m++ )
        {
            int l = nd->down[m];
            if ( HQ[l] > IS_FLOW_TOL )
            {
                if ( IsPipe[l].kind == IS_LINK_PIPE ) dh = HY[l] + IsPipe[l].fromInv;
                else if ( IsPipe[l].kind == IS_LINK_FM ) dh = IsPipe[l].dia + IsPipe[l].fromInv;
            }
            h = MAX(h, dh);
        }
        minorLoss(k, &h);
        if ( nd->split ) HHyd[k] = MAX(HHyd[k], h);
        else             HHyd[k] = h;
        HOut[k] = MIN(HHyd[k], nd->cap);   // unlocked: cap = rim
        break;

      case IS_NODE_CHAMBER:
        minorLoss(k, &h);
        HOut[k] = h;
        break;

      case IS_NODE_WETWELL:
        HOut[k] = HWell[nd->wetwell];      // head is not carried upstream
        break;

      default:
        HError = 1;
        return;
    }
    HHyd[k] = HOut[k];
}

static double adjustPipe(int l, double dh, double uh)
//  Fv_AdjustPipeHead he_runhead.cpp:1591.
//  QUIRK: backwater below the crown adds full-length full-pipe friction loss.
{
    TIsPipe* p = &IsPipe[l];
    double L_, D, s, ti, downRef, upH, deltaLevel, angle, refLength, dH, crossing;
    if ( HQ[l] >= (p->fullFlow - IS_FLOW_TOL) ) return uh;
    L_ = p->length;  D = p->dia;  s = p->slope;  ti = p->toInv;
    downRef = ti + D;
    if ( HQ[l] > IS_FLOW_TOL )
    {
        if ( dh > uh && dh < downRef )
        {
            upH = dh + isewer_headloss(l, HQ[l]);
            if ( upH > uh )
            {
                uh = upH;
                HPartial[l] = 1;
                if ( uh >= p->fromInv + D ) HAdjFull[l] = 1;
            }
        }
        if ( dh <= downRef + IS_DEPTH_TOL ) return uh;
        deltaLevel = dh - (ti + HY[l]);
        angle = atan(s);
        refLength = deltaLevel / sin(angle);
        dH = 0.0;
        if ( refLength > L_ ) dH = isewer_headloss(l, HQ[l]);
        crossing = dh + dH;
        if ( dh > ti ) HPartial[l] = 1;
        if ( crossing > uh )
        {
            uh = crossing;
            HPartial[l] = 1;
            if ( uh >= p->fromInv + D ) HAdjFull[l] = 1;
        }
    }
    return uh;
}

static void pipeHead(int l, int uk, double dh)
//  Fv_GetPipeHead he_runhead.cpp:1404 (circular)
//  QUIRK: zero flow "by invert" reads hm_pLoadNode[upstream node index]; the
//  node and load node indexes agree for load manholes, which come first.
{
    TIsPipe* p = &IsPipe[l];
    double uh, dH, ref, ui, di, D, dd;
    if ( p->slopeType == IS_SLOPE_STANDARD )
    {
        if ( HQ[l] < p->fullFlow )
        {
            if ( HQ[l] <= IS_FLOW_TOL ) uh = MAX(dh, IsNode[uk].bottomRef);  // by invert
            else uh = HY[l] + p->fromInv;
            uh = adjustPipe(l, dh, uh);
        }
        else
        {
            dH = isewer_headloss(l, HQ[l]);
            ref = p->dia + p->toInv;
            uh = MAX(dh, ref) + dH;
        }
    }
    else
    {
        ui = p->fromInv;  di = p->toInv;
        if ( HQ[l] > IS_FLOW_TOL )
        {
            D = p->dia;
            dH = (HY[l] >= D) ? isewer_headloss(l, HQ[l]) : 0.0;
            ref = HY[l] + di;
            dd = MAX(dh, ref);
            uh = MAX(ui + HY[l], dd + dH);
            dh = dd;
            uh = adjustPipe(l, dh, uh);
        }
        else
        {
            uh = ui;
            dh = MAX(dh, di);
            if ( dh > di ) HPartial[l] = 1;
            if ( dh > di + p->dia ) HAdjFull[l] = 1;
        }
    }
    adjustNode(uk, uh);
    uh = HHyd[uk];
    HFtype[l] = (uh >= (p->fromInv + p->dia + IS_DEPTH_TOL)) ? 1 : 0;
}

static void forceMainHead(int l, double dh)
//  Fv_GetForceMainHead he_runhead.cpp:1722: heads along the whole siphon
//  chain, from its tail up.
{
    int s, i, nfm;
    double ref = dh, stationFlow = 0.0, linkFlow;
    if ( HFmFlag[l - IsNP] ) return;
    s = IsPipe[l].siphon;
    if ( s < 0 ) { HError = 1;  return; }
    nfm = IsSiphon[s].nfm;
    if ( IsSiphon[s].type == IS_SIPHON_STATION )
    {
        TIsStation* st = &IsStation[IsSiphon[s].station];
        double stationHead = HPumpH[st->pump[0]];
        double siphonHead, down;
        stationFlow = HQ[l];               // flow of the force main given
        siphonHead = siphonHeadloss(s, stationFlow);
        down = HWell[st->wetwell] + stationHead - siphonHead;
        if ( down > dh )
        {
            adjustNode(IsPipe[IsSiphon[s].fm[nfm - 1]].to, down);
            ref = down;
        }
    }
    linkFlow = stationFlow;
    for ( i = nfm - 1; i >= 0; i-- )
    {
        int f = IsSiphon[s].fm[i], uk = IsPipe[f].frm;
        double uh;
        if ( IsSiphon[s].type == IS_SIPHON_STANDARD ) linkFlow = HQ[f];
        uh = ref + fmHeadloss(f, linkFlow);
        if ( linkFlow > IS_FLOW_TOL ) uh = MAX(uh, IsPipe[f].fromInv + IsPipe[f].dia);
        else uh = IsPipe[f].fromInv;       // zero flow by invert
        adjustNode(uk, uh);
        ref = HHyd[uk];
        HFmFlag[f - IsNP] = 1;
    }
}

static void stationHead(int l, double dh)
//  Fv_GetStationHead he_runhead.cpp:1871
{
    TIsStation* st = &IsStation[IsPipe[l].station];
    adjustNode(st->frm, dh - HPumpH[st->pump[0]]);
}

static void nodeHead(int l, int uk, double dh)
//  Fv_GetNodeHead he_runhead.cpp:1353
{
    switch ( IsPipe[l].kind )
    {
      case IS_LINK_PIPE: pipeHead(l, uk, dh);  break;
      case IS_LINK_FM:   forceMainHead(l, dh); break;
      default:           stationHead(l, dh);   break;
    }
}

static void depthSearch(int k)
//  Fv_DepthSearch he_runhead.cpp:1162. A split node is suspended here and
//  resumed level by level from isewer_runHead.
{
    TIsNode* nd = &IsNode[k];
    int m;
    if ( fully(k) ) return;
    if ( nd->split )
    {
        if ( HSplitOrder[k] < 0 )
        {
            HSplitOrder[k] = HCounter;
            HSplitNode[HCounter] = k;
            HCounter++;
        }
        return;
    }
    for ( m = 0; m < nd->nuph; m++ )
    {
        int l = nd->uph[m];
        int nk = IsPipe[l].frm;
        nodeHead(l, nk, HHyd[k]);
        HLevel[k]++;
        depthSearch(nk);
    }
}

static void splitDepthSearch(int k)
//  Fv_FlowSplitDepthSearch he_runhead.cpp:1239
{
    TIsNode* nd = &IsNode[k];
    int m;
    if ( fully(k) ) return;
    for ( m = 0; m < nd->nuph; m++ )
    {
        int l = nd->uph[m];
        int nk = IsPipe[l].frm;
        nodeHead(l, nk, HHyd[k]);
        HLevel[k]++;
        depthSearch(nk);
    }
}

int isewer_runHead(const double* Q, const double* y, double* head,
                   int* ftype, int* partial, int* adjFull, int resetFtype)
//  He_RunHead::F_Run he_runhead.cpp:554. Heads from each outlet upstream.
//  Returns 1 if the engine's sweep would fail (F_Run false).
//  QUIRK: nested flow splits can leave a stale head upstream of a split.
//  QUIRK: in EPS the flow type flags are never reset between report times
//  (resetFtype = 0); a pipe the sweep does not reach keeps its old flag.
{
    int k, i, N = IsNN;
    int nsplit = 0, remain, start, end, tmpStart;
    double ff;

    HQ = Q;  HY = y;  HOut = head;
    HFtype = ftype;  HPartial = partial;  HAdjFull = adjFull;
    HError = 0;
    HHyd = (double*) calloc(N + 1, sizeof(double));
    HLevel = (int*) calloc(N + 1, sizeof(int));
    HSplitOrder = (int*) calloc(N + 1, sizeof(int));
    HSplitNode = (int*) calloc(N + 1, sizeof(int));
    HFmFlag = (int*) calloc(IsNF + 1, sizeof(int));
    if ( !HHyd || !HLevel || !HSplitOrder || !HSplitNode || !HFmFlag )
    {
        FREE(HHyd); FREE(HLevel); FREE(HSplitOrder); FREE(HSplitNode); FREE(HFmFlag);
        report_writeErrorMsg(ERR_MEMORY, "");
        return 1;
    }
    for ( k = 0; k < N; k++ )
    {
        HHyd[k] = -IS_BIG;
        HOut[k] = 0.0;            // never initialised in C++; unreached nodes
        HSplitOrder[k] = -1;
        if ( IsNode[k].split ) nsplit++;
    }
    for ( i = 0; i < IsNP; i++ )
    {
        if ( resetFtype ) HFtype[i] = 0;
        HPartial[i] = 0;  HAdjFull[i] = 0;
    }
    HCounter = 0;

    for ( k = 0; k < N && !HError; k++ )
    {
        TIsNode* nd = &IsNode[k];
        int m;
        if ( nd->kind != IS_NODE_OUTLET ) continue;
        // Fv_GetOutletHead he_runhead.cpp:1296
        ff = -IS_BIG;
        for ( m = 0; m < nd->nuph; m++ )
        {
            int l = nd->uph[m];
            if ( IsPipe[l].kind != IS_LINK_PIPE ) { HError = 1;  break; }
            ff = MAX(ff, HY[l] + IsPipe[l].toInv);
        }
        HHyd[k] = MAX(nd->submerge, ff);
        HOut[k] = HHyd[k];
        depthSearch(k);
    }

    // --- split nodes, level by level (QUIRK: remain counts levels)
    remain = nsplit;
    start = 0;
    while ( remain > 0 && !HError )
    {
        tmpStart = HCounter;
        end = HCounter;
        if ( start >= end ) break;
        for ( i = start; i < end; i++ ) splitDepthSearch(HSplitNode[i]);
        remain--;
        start = tmpStart;
    }
    FREE(HHyd); FREE(HLevel); FREE(HSplitOrder); FREE(HSplitNode); FREE(HFmFlag);
    return HError;
}

//=============================================================================
//  SWMM state, loads and execution
//=============================================================================

static void setLinkState(int j, double q, double y)
{
    int k = Link[j].subIndex;
    double a;
    if ( Link[j].type != CONDUIT )
    {
        Link[j].newFlow = q;
        Link[j].newDepth = 0.0;
        Link[j].newVolume = 0.0;
        return;
    }
    a = (y > 0.0) ? xsect_getAofY(&Link[j].xsect, y) : 0.0;
    Link[j].newFlow = q;
    Link[j].newDepth = y;
    Conduit[k].q1 = Conduit[k].q2 = q / Conduit[k].barrels;
    Conduit[k].a1 = Conduit[k].a2 = a;
    Link[j].newVolume = a * link_getLength(j) * Conduit[k].barrels;
}

static void setNodeState(int j, double head)
{
    double y = head - Node[j].invertElev;
    if ( y < 0.0 ) y = 0.0;
    Node[j].newDepth = y;
    Node[j].newVolume = (Node[j].type == STORAGE) ? node_getVolume(j, y) : 0.0;
    Node[j].overflow = 0.0;
}

static double steadyLinkFlow(int j, double* depth)
//  InfoSewer steady flow (cfs) and depth (ft) of SWMM link j. A force main
//  runs full. A pump carries its element flow; a capacity pump's element
//  flow is one unit (station total = parallel x capacity), a curve pump's is
//  the flow of all its parallel units (head = head0 - a*Q^b, a built on
//  parallel x flow1).
{
    int l = IsLinkOfSwmm[j];
    if ( l < IsNP ) { *depth = SDepth[l];  return MAX(SQ[l], 0.0); }
    if ( l < IsNPM ) { *depth = IsPipe[l].dia;  return MAX(SQ[l], 0.0); }
    *depth = 0.0;
    {
        TIsPump* pu = &IsPump[l - IsNPM];
        double q = MAX(SQ[l], 0.0);
        if ( pu->type == IS_PUMP_CAPACITY ) q *= pu->parallel;
        return q;
    }
}

static void setSteadyState(void)
//  Steady state: SWMM carries the InfoSewer solution unchanged at every step.
//  Each node receives a lateral inflow equal to InfoSewer's flow increment
//  there (outflow - inflow, negative where peaking lowers a combined flow),
//  so SWMM's continuity closes on the peaked flows.
{
    int i, j;
    double y, q;
    for ( j = 0; j < Nobjects[NODE]; j++ ) SInc[j] = 0.0;
    for ( i = 0; i < Nobjects[LINK]; i++ )
    {
        q = steadyLinkFlow(i, &y);
        setLinkState(i, q, y);
        if ( Node[Link[i].node1].type != OUTFALL ) SInc[Link[i].node1] += q;
        if ( Node[Link[i].node2].type != OUTFALL ) SInc[Link[i].node2] -= q;
    }
    for ( j = 0; j < Nobjects[NODE]; j++ ) setNodeState(j, SHead[IsNodeOfSwmm[j]]);
}

void isewer_addLoads(DateTime currentDate)
//  Called by addSystemInflows (routing.c) after dry weather inflows.
{
    int j;
    double q;
    if ( !IsNode ) return;
    if ( IsAnalysis == IS_EPS )
    {
        // loads at the end of the step, where isewer_execute sets the state
        (void)currentDate;
        isewer_epsLoads(NewRoutingTime / 1000.0);
        return;
    }
    for ( j = 0; j < Nobjects[NODE]; j++ )
    {
        q = SInc ? SInc[j] : 0.0;
        if ( q == 0.0 ) continue;
        Node[j].newLatFlow += q;
        massbal_addInflowFlow(DRY_WEATHER_INFLOW, q);
    }
}

int isewer_execute(int links[], double tStep)
//  Called by flowrout_execute: sets the InfoSewer state for the end of the
//  current step and the node inflow/outflow totals SWMM's mass balance uses.
{
    int i, j;
    double y;
    (void)links;  (void)tStep;
    if ( !IsNode ) return 0;
    if ( IsAnalysis == IS_EPS )
    {
        // NewRoutingTime is already the end of this step (routing.c)
        isewer_epsSetState(NewRoutingTime / 1000.0);
        isewer_epsNodeFlows();
        return 1;
    }
    for ( i = 0; i < Nobjects[LINK]; i++ )
    {
        double q = steadyLinkFlow(i, &y);
        setLinkState(i, q, y);
    }
    for ( j = 0; j < Nobjects[NODE]; j++ ) setNodeState(j, SHead[IsNodeOfSwmm[j]]);
    for ( i = 0; i < Nobjects[LINK]; i++ )
    {
        Node[Link[i].node1].outflow += Link[i].newFlow;
        Node[Link[i].node2].inflow  += Link[i].newFlow;
    }
    return 1;
}

double isewer_getVelocity(int j, double flow, double depth)
//  Velocity as InfoSewer computes it, used by link_getVelocity for SWMM's
//  .out and report statistics. Gravity pipes: per barrel flow over the
//  3.1416 circle area (HHU_GetPipeVelocity). Force mains: flow over the
//  3.1415926 full area (F_GetHydraulicData).
{
    int k = Link[j].subIndex, l = IsLinkOfSwmm ? IsLinkOfSwmm[j] : j;
    double r = 0.5 * Link[j].xsect.yFull;
    if ( Link[j].xsect.type == FORCE_MAIN )
    {
        double D = Link[j].xsect.yFull;
        if ( IsAnalysis == IS_STEADY && SVel ) return SVel[l];
        return flow / (0.25 * IS_PI_DEF * D * D);
    }
    if ( depth <= 0.0 ) return 0.0;
    if ( depth > Link[j].xsect.yFull ) depth = Link[j].xsect.yFull;
    return isewer_pipeVelocity(depth, r, flow / Conduit[k].barrels);
}

//=============================================================================
//  Steady-state report and results file
//=============================================================================

typedef struct
{
    double flow, unpeak, peak, cover, infil, flowtype, full, rescap, vel, dD, qQ,
           depth, crit, froude, covcnt, backwater, adjdepth, adjvel, partial;
}  TIsPipeOut;

static void pipeOutputs(int i, TIsPipeOut* o)
//  Fv_PipeData he_staticrun.cpp:1244-1460
{
    TIsPipe* p = &IsPipe[i];
    double Qc = MAX(SQ[i], 0.0), D = p->dia, d = SDepth[i], T, A, uh, dh, dse, a;
    o->flow = Qc;
    o->unpeak = SUP[i] - SINF[i];
    o->peak = SPK[i];
    o->cover = SCV[i];
    o->infil = SINF[i];
    o->flowtype = SFtype[i] == 1 ? 1.0 : -1.0;
    o->full = p->fullFlow;
    o->rescap = MAX(p->fullFlow - Qc, 0.0);
    T = isewer_wetWidth(d, D / 2.0);
    A = isewer_wetArea(d, D / 2.0);
    if ( T <= 0.0 && d > 0.0 ) o->froude = SVel[i] / sqrt(32.2 * d);
    else if ( A > 0.0 ) o->froude = SVel[i] / sqrt(32.2 * (A / T));
    else o->froude = 0.0;
    o->vel = SVel[i];
    o->dD = d / D;
    o->qQ = Qc / isewer_linkFlow(i, D, -1.0, -1.0);
    o->depth = d;
    o->crit = SCrit[i];
    o->covcnt = 1000 * SCC[i];
    uh = SHead[p->frm];  dh = SHead[p->to];
    dse = d + p->toInv;
    if ( dh > dse )
    {
        a = 0.5 * ((dh - p->toInv) + (uh - p->fromInv));
        a = MIN(D, a);
        o->adjdepth = a;
        o->adjvel = isewer_pipeVelocity(a, D / 2, Qc);
        o->backwater = 1.0;
    }
    else
    {
        o->adjdepth = d;  o->adjvel = SVel[i];  o->backwater = 0.0;
    }
    o->partial = SPartial[i] ? 1.0 : -1.0;
}

static double nodeBase(int k)
{
    int m;
    double b = 0.0;
    for ( m = 0; m < IsNode[k].nload; m++ ) b += IsLoads[IsNode[k].load[m]].base;
    return b;
}

static int nodeJump(int k, double* fr)
//  he_staticrun.cpp:1326-1369: supercritical upstream, subcritical downstream
//  (gravity pipes only)
{
    int m, upSuper = 0;
    for ( m = 0; m < IsNode[k].nup; m++ )
    {
        int l = IsNode[k].up[m];
        if ( IsPipe[l].kind == IS_LINK_PIPE && fr[l] > 1 ) { upSuper = 1;  break; }
    }
    if ( !upSuper ) return 0;
    for ( m = 0; m < IsNode[k].ndown; m++ )
    {
        int l = IsNode[k].down[m];
        if ( IsPipe[l].kind == IS_LINK_PIPE && fr[l] < 1 ) return 1;
    }
    return 0;
}

static double outletFlow(int k)
//  Fv_GetOutletFlow: sums the UNclipped flows of the links into the outlet
{
    int i;
    double q = 0.0;
    for ( i = 0; i < IsNPM; i++ ) if ( IsPipe[i].to == k ) q += SQ[i];
    return q;
}

static double pumpOut(int p, int what)
//  Fv_PumpData he_staticrun.cpp:1681: 0 flow, 1 head, 2 power (hp),
//  3 usage, 4 speed
{
    double q = SQ[IsNPM + p], h = SPumpHead[p];
    int on = fabs(q) > 1.0e-7;
    switch ( what )
    {
      case 0:  return q;
      case 1:  return h;
      case 2:  return on ? MAX(0.0, q * h * (62.4 / 550.0)) : 0.0;
      case 3:  return on ? (double)IsPump[p].parallel : 0.0;
      default: return on ? IsPump[p].speed : 0.0;
    }
}

#define NODE_ID(k)  Node[IsNode[k].swmm].ID
#define LINK_ID(l)  Link[IsPipe[l].swmm].ID

static void writeSteadyReport(void)
{
    int i, k, N = IsNN;
    double ff = IsFlowFactor, lf = UCF(LENGTH);
    TIsPipeOut o;
    FILE* f = Frpt.file;
    if ( !f || RptFlags.disabled ) return;

    fprintf(f,
"\n  ***********************************"
"\n  InfoSewer Steady-State Pipe Results"
"\n  ***********************************"
"\n  Flows in %s; peaked by type; d/D from Manning normal depth.\n",
    FlowUnitWords[FlowUnits]);
    fprintf(f,
"\n  ------------------------------------------------------------------------------------------------------------"
"\n                         Total     Unpeak       Peak      Cover      Infil   Velocity       d/D        q/Q   Surch"
"\n  Pipe                    Flow       Flow       Flow       Flow       Flow     %-6s                            "
"\n  ------------------------------------------------------------------------------------------------------------",
    UnitSystem == US ? "ft/s" : "m/s");
    for ( i = 0; i < IsNP; i++ )
    {
        pipeOutputs(i, &o);
        fprintf(f, "\n  %-16s %10.4f %10.4f %10.4f %10.4f %10.4f %10.3f %9.3f %10.3f   %s",
            LINK_ID(i), o.flow * ff, o.unpeak * ff, o.peak * ff, o.cover * ff,
            o.infil * ff, o.vel * lf, o.dD, o.qQ, o.flowtype > 0 ? "YES" : "no");
    }
    if ( IsNF > 0 )
    {
        fprintf(f,
"\n\n  *****************************************"
"\n  InfoSewer Steady-State Force Main Results"
"\n  *****************************************\n");
        fprintf(f,
"\n  ---------------------------------------------------------------------------------------"
"\n                         Total     Unpeak       Peak      Cover      Infil   Velocity  Headloss"
"\n  Force Main              Flow       Flow       Flow       Flow       Flow     %-6s     %-4s"
"\n  ---------------------------------------------------------------------------------------",
        UnitSystem == US ? "ft/s" : "m/s", UnitSystem == US ? "ft" : "m");
        for ( i = IsNP; i < IsNPM; i++ )
            fprintf(f, "\n  %-16s %10.4f %10.4f %10.4f %10.4f %10.4f %10.3f %9.3f",
                LINK_ID(i), SQ[i] * ff, (SUP[i] - SINF[i]) * ff, SPK[i] * ff,
                SCV[i] * ff, SINF[i] * ff, SVel[i] * lf, SFmLoss[i - IsNP] * lf);
    }
    if ( IsNPump > 0 )
    {
        fprintf(f,
"\n\n  ***********************************"
"\n  InfoSewer Steady-State Pump Results"
"\n  ***********************************"
"\n  Curve pump flow is for all parallel units; capacity pump flow is per unit.\n");
        fprintf(f,
"\n  -----------------------------------------------------------------------"
"\n                         Flow       Head      Power    Units      Speed"
"\n  Pump               %9s     %6s         hp"
"\n  -----------------------------------------------------------------------",
        FlowUnitWords[FlowUnits], UnitSystem == US ? "ft" : "m");
        for ( i = 0; i < IsNPump; i++ )
            fprintf(f, "\n  %-16s %10.4f %10.3f %10.3f %8.0f %10.3f", Link[IsPump[i].swmm].ID,
                pumpOut(i, 0) * ff, pumpOut(i, 1) * lf, pumpOut(i, 2), pumpOut(i, 3),
                pumpOut(i, 4));
    }
    fprintf(f,
"\n\n  **************************************"
"\n  InfoSewer Steady-State Manhole Results"
"\n  **************************************\n");
    fprintf(f,
"\n  -------------------------------------------------------------------"
"\n                         Base Load    Total Load        HGL     Surch"
"\n  Manhole                %9s     %9s     %6s"
"\n  -------------------------------------------------------------------",
    FlowUnitWords[FlowUnits], FlowUnitWords[FlowUnits], UnitSystem == US ? "ft" : "m");
    for ( k = 0; k < N; k++ )
    {
        if ( IsNode[k].kind != IS_NODE_LOAD ) continue;
        fprintf(f, "\n  %-16s %13.4f %13.4f %10.3f     %s", NODE_ID(k),
            nodeBase(k) * ff, SStatic[k] * ff, SHead[k] * lf,
            SHead[k] >= IsNode[k].rim ? "YES" : "no");
    }
    for ( k = 0; k < N; k++ )
    {
        const char* what;
        if ( IsNode[k].kind == IS_NODE_CHAMBER ) what = "chamber";
        else if ( IsNode[k].kind == IS_NODE_WETWELL ) what = "wet well";
        else continue;
        fprintf(f, "\n  %-16s %13s %13s %10.3f     %s", NODE_ID(k), what, "",
            SHead[k] * lf, "");
    }
    fprintf(f, "\n");
}

static void writeSteadyCsv(void)
//  Every InfoSewer steady-state output array at full precision, in the
//  engine's internal units (cfs, ft) and its output order, for verification.
{
    int i, k, N = IsNN;
    double *fr;
    TIsPipeOut o;
    FILE* f = fopen(IsResultsFile, "wt");
    if ( !f ) { report_writeErrorMsg(ERR_RPT_FILE, IsResultsFile); return; }
    fr = (double*) calloc(IsNLE + 1, sizeof(double));
    fprintf(f, "kind,id,flow,unpeak,peak,cover,infil,storm,flowtype,full,rescap,vel,"
               "dD,qQ,depth,crit,froude,covcnt,backwater,adjdepth,adjvel,partial\n");
    for ( i = 0; i < IsNP; i++ )
    {
        pipeOutputs(i, &o);
        if ( fr ) fr[i] = o.froude;
        fprintf(f, "pipe,%s,%.17g,%.17g,%.17g,%.17g,%.17g,0,%.17g,%.17g,%.17g,%.17g,"
                   "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            LINK_ID(i), o.flow, o.unpeak, o.peak, o.cover, o.infil, o.flowtype,
            o.full, o.rescap, o.vel, o.dD, o.qQ, o.depth, o.crit, o.froude,
            o.covcnt, o.backwater, o.adjdepth, o.adjvel, o.partial);
    }
    fprintf(f, "kind,id,base,static,nstorm,head,surch,jump\n");
    for ( k = 0; k < N; k++ )
    {
        if ( IsNode[k].kind != IS_NODE_LOAD ) continue;
        fprintf(f, "node,%s,%.17g,%.17g,0,%.17g,%.17g,%.17g\n", NODE_ID(k),
            nodeBase(k), SStatic[k], SHead[k],
            SHead[k] >= IsNode[k].rim ? 1.0 : 0.0,
            fr ? (double)nodeJump(k, fr) : 0.0);
    }
    fprintf(f, "kind,id,chead\n");
    for ( k = 0; k < N; k++ )
        if ( IsNode[k].kind == IS_NODE_CHAMBER )
            fprintf(f, "chamber,%s,%.17g\n", NODE_ID(k), SHead[k]);
    fprintf(f, "kind,id,oflow,ohead\n");
    for ( k = 0; k < N; k++ )
        if ( IsNode[k].kind == IS_NODE_OUTLET )
            fprintf(f, "outlet,%s,%.17g,%.17g\n", NODE_ID(k), outletFlow(k), SHead[k]);
    fprintf(f, "kind,id,whead\n");
    for ( k = 0; k < N; k++ )
        if ( IsNode[k].kind == IS_NODE_WETWELL )
            fprintf(f, "wetwell,%s,%.17g\n", NODE_ID(k), SHead[k]);
    fprintf(f, "kind,id,fflow,funpeak,fpeak,fcover,finfil,fstorm,fvel,floss\n");
    for ( i = IsNP; i < IsNPM; i++ )
        fprintf(f, "fm,%s,%.17g,%.17g,%.17g,%.17g,%.17g,0,%.17g,%.17g\n", LINK_ID(i),
            SQ[i], SUP[i] - SINF[i], SPK[i], SCV[i], SINF[i], SVel[i], SFmLoss[i - IsNP]);
    fprintf(f, "kind,id,pflow,phead,ppower,pusage,pspeed\n");
    for ( i = 0; i < IsNPump; i++ )
        fprintf(f, "pump,%s,%.17g,%.17g,%.17g,%.17g,%.17g\n", Link[IsPump[i].swmm].ID,
            pumpOut(i, 0), pumpOut(i, 1), pumpOut(i, 2), pumpOut(i, 3), pumpOut(i, 4));
    FREE(fr);
    fclose(f);
}

//=============================================================================

static void freeAll(void)
{
    int p;
    FREE(SQ); FREE(SUP); FREE(SPK); FREE(SCV); FREE(SCC); FREE(SINF);
    FREE(SDepth); FREE(SVel); FREE(SCrit); FREE(SFull); FREE(SFtype);
    FREE(SPartial); FREE(SAdjFull); FREE(SStatic); FREE(SHead); FREE(SInc);
    FREE(SFmLoss); FREE(SPumpHead); FREE(SWellHead);
    FREE(SPumpHead0); FREE(SPumpA); FREE(SPumpB); FREE(SPumpPB); FREE(SPumpOn);
    FREE(SPumpSpeed);
    if ( IsStation ) for ( p = 0; p < IsNStation; p++ ) FREE(IsStation[p].pump);
    if ( IsSiphon )  for ( p = 0; p < IsNSiphon; p++ )  FREE(IsSiphon[p].fm);
    FREE(IsStation); FREE(IsSiphon); FREE(IsPump); FREE(IsWell);
    FREE(IsNode); FREE(IsPipe); FREE(IsOrder); FREE(IsAllocList); FREE(IsRim);
    FREE(IsFullDepth0); FREE(IsLinkInv); FREE(IsNodeOfSwmm); FREE(IsLinkOfSwmm);
    for ( p = 0; p < IsNpats; p++ ) FREE(IsPats[p].f);
    for ( p = 0; p < IsNcurves; p++ ) { FREE(IsCurves[p].x); FREE(IsCurves[p].y); }
    FREE(IsPats); FREE(IsCurves); FREE(IsLoads); FREE(SplitIn); FREE(InfilIn);
    FREE(ChamberIn); FREE(WellIn); FREE(PumpIn);
    IsNpats = MaxPats = IsNcurves = MaxCurves = IsNloads = MaxLoads = 0;
    NsplitIn = MaxSplitIn = NinfilIn = MaxInfilIn = 0;
    NchamberIn = MaxChamberIn = NwellIn = MaxWellIn = NpumpIn = MaxPumpIn = 0;
    IsNN = IsNP = IsNF = IsNPM = IsNPump = IsNWell = IsNStation = IsNSiphon = 0;
    IsNT = IsNLE = 0;
}

void isewer_close(void)
{
    if ( IsAnalysis == IS_EPS ) isewer_epsClose();
    freeAll();
}
