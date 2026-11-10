//-----------------------------------------------------------------------------
//   isewer.c
//
//   Project:  swmm5-infosewer (EPA SWMM 5.2.4 + InfoSewer gravity routing)
//   Purpose:  InfoSewer's gravity-sewer calculations as a SWMM5 flow routing
//             option:  FLOW_ROUTING INFOSEWER
//               INFOSEWER_ANALYSIS STEADY  - InfoSewer steady state (loads by
//                                            type, peaking, normal depth, HGL)
//               INFOSEWER_ANALYSIS EPS     - InfoSewer extended period routing
//                                            (isewer_eps.c)
//
//   This file: input sections and options, the InfoSewer network built from
//   SWMM's conduits, junctions and outfalls, the steady-state solver, the HGL
//   sweep shared with EPS, and the SWMM state and report hooks.
//
//   The calculations follow the InfoSewer engine (HeEngPkg, MWH Soft/Innovyze,
//   1996-2011) as written up in spec/steady.md and spec/eps.md, which cite the
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
char    IsResultsFile[MAXFNAME+1];  // optional full-precision CSV

TIsLoad*    IsLoads;   int IsNloads;   static int MaxLoads;
TIsPattern* IsPats;    int IsNpats;    static int MaxPats;
static TIsCurve* IsCurves; static int IsNcurves, MaxCurves;

TIsNode* IsNode;
TIsPipe* IsPipe;
int*     IsOrder;
static int  IsNorder;
static double* IsRim;        // junction rim saved before SWMM raises full depth

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
static int    buildTopology(void);
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
    IsResultsFile[0] = '\0';
    IsLoads = NULL;  IsNloads = MaxLoads = 0;
    IsPats = NULL;   IsNpats = MaxPats = 0;
    IsCurves = NULL; IsNcurves = MaxCurves = 0;
    SplitIn = NULL;  NsplitIn = MaxSplitIn = 0;
    InfilIn = NULL;  NinfilIn = MaxInfilIn = 0;
    IsNode = NULL;  IsPipe = NULL;  IsOrder = NULL;  IsRim = NULL;
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
    if ( match(k, "RESULTS_FILE") )
    {
        sstrncpy(IsResultsFile, s2, MAXFNAME);
        return 0;
    }
    if ( match(k, "PATTERN_STEP") )
    {
        DateTime t;
        if ( getDouble(s2, &x) ) { IsPatternStep = x; return 0; }
        if ( !datetime_strToTime(s2, &t) ) return error_setInpError(ERR_DATETIME, s2);
        IsPatternStep = floor(t * SECperDAY + 0.5);
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
//  manhole heads at the rim the user gave, so keep it.
{
    int j;
    FREE(IsRim);
    IsRim = (double*) calloc(Nobjects[NODE] + 1, sizeof(double));
    if ( !IsRim ) { report_writeErrorMsg(ERR_MEMORY, ""); return; }
    for ( j = 0; j < Nobjects[NODE]; j++ )
        IsRim[j] = Node[j].invertElev + Node[j].fullDepth;
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

static int buildNetwork(void)
//  InfoSewer nodes and pipes from SWMM junctions, outfalls and conduits.
//  Mirrors he_input.cpp: Fv_ReadPipe (slope, bottom reference),
//  Fv_UnitConversion, and the option conversions of he_input.cpp:12143-12339.
{
    int    i, j, k, n1, n2, N = Nobjects[NODE], L = Nobjects[LINK];
    double lenF, diaF;

    if ( UnitSystem == US ) { lenF = 1.0;    diaF = 12.0; }
    else                    { lenF = 0.3048; diaF = 304.8; }

    IsNode = (TIsNode*) calloc(N + 1, sizeof(TIsNode));
    IsPipe = (TIsPipe*) calloc(L + 1, sizeof(TIsPipe));
    IsOrder = (int*) calloc(N + 1, sizeof(int));
    IsAllocList = (int*) calloc(8 * (L + 1) + IsNloads + 1, sizeof(int));  // conn, down, up, uph: 2L each
    if ( !IsNode || !IsPipe || !IsOrder || !IsAllocList ) return netError(ERR_MEMORY, "");

    // --- nodes
    for ( j = 0; j < N; j++ )
    {
        TIsNode* nd = &IsNode[j];
        if ( Node[j].type == JUNCTION ) nd->isOutlet = 0;
        else if ( Node[j].type == OUTFALL ) nd->isOutlet = 1;
        else return netError(ERR_ISEWER_NODE, Node[j].ID);
        nd->rim = IsRim ? IsRim[j] : Node[j].invertElev + Node[j].fullDepth;
        nd->cap = nd->rim + Node[j].surDepth;
        nd->bottomRef = IS_BIG;
        nd->submerge = -IS_BIG;
        if ( nd->isOutlet )
        {
            k = Node[j].subIndex;
            if ( Outfall[k].type == FIXED_OUTFALL ) nd->submerge = Outfall[k].fixedStage;
        }
    }

    // --- pipes
    for ( i = 0; i < L; i++ )
    {
        TIsPipe* p = &IsPipe[i];
        if ( Link[i].type != CONDUIT ) return netError(ERR_ISEWER_LINK, Link[i].ID);
        if ( Link[i].xsect.type != CIRCULAR ) return netError(ERR_ISEWER_LINK, Link[i].ID);
        k = Link[i].subIndex;
        n1 = Link[i].node1;  n2 = Link[i].node2;
        p->frm = n1;  p->to = n2;
        p->fromInv = Node[n1].invertElev + Link[i].offset1;
        p->toInv   = Node[n2].invertElev + Link[i].offset2;
        p->length  = Conduit[k].length;
        p->dia     = Link[i].xsect.yFull;
        p->n       = Conduit[k].roughness;
        p->count   = Conduit[k].barrels;
        p->kLoss   = fabs(Link[i].cLossInlet);
        p->split   = IS_SPLIT_HYD;
        p->splitPct = 0.0;
        p->infType = IS_INF_NONE;  p->infRate = 0.0;  p->infCount = 0.0;
        p->infPat  = -1;

        // Fv_UpdateSlopeData he_data.cpp:5600
        p->slope = (p->fromInv - p->toInv) / p->length;
        if ( fabs(p->slope) < IS_SYS_ZERO_FLOAT ) p->slopeType = IS_SLOPE_FLAT;
        else if ( p->slope > 0.0 )                p->slopeType = IS_SLOPE_STANDARD;
        else                                      p->slopeType = IS_SLOPE_ADVERSE;
        p->adverse = p->fromInv < p->toInv;       // Fv_CheckPipeInvert

        // bottom reference levels (he_input.cpp:5256-5270)
        if ( !IsNode[n1].isOutlet )
            IsNode[n1].bottomRef = MIN(IsNode[n1].bottomRef, p->fromInv);
        if ( !IsNode[n2].isOutlet )
            IsNode[n2].bottomRef = MIN(IsNode[n2].bottomRef, p->toInv);
    }

    // --- percent splits (Fv_ReadFlowSplit; F_SetSplitFlowPercent)
    for ( i = 0; i < NsplitIn; i++ )
    {
        TIsPipe* p = &IsPipe[SplitIn[i].link];
        p->splitPct = SplitIn[i].pct;
        if ( p->splitPct > 0.0 ) p->split = IS_SPLIT_PCT;
    }

    // --- infiltration (Fv_ReadInfiltrateRate; HD_PipeData::F_ConvertUnit)
    for ( i = 0; i < NinfilIn; i++ )
    {
        TIsPipe* p = &IsPipe[InfilIn[i].link];
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
        IsLoads[i].base = IsLoads[i].baseIn / IsFlowFactor;
        if ( IsLoads[i].pattern[0] )
        {
            IsLoads[i].pat = findPattern(IsLoads[i].pattern);
            if ( IsLoads[i].pat < 0 ) return netError(ERR_NAME, IsLoads[i].pattern);
        }
        if ( IsNode[IsLoads[i].node].isOutlet )
            return netError(ERR_ISEWER_NODE, Node[IsLoads[i].node].ID);
    }

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

    // --- flow tolerances (he_input.cpp:12321, 12339)
    IsFlowTol = IsFlowTolIn / IsFlowFactor;
    IsIterTol = IS_ITER_FLOW_TOL;
    if ( IsIterTol > IsFlowTol ) IsIterTol = IsFlowTol;

    // --- full flow of each pipe x barrels (HH_PipeCapacity)
    for ( i = 0; i < L; i++ ) IsPipe[i].fullFlow = isewer_fullFlowBarrel(i) * IsPipe[i].count;
    return 0;
}

static int buildTopology(void)
//  Fv_SetLinkConnection, Fv_SetFlowSplitLoadNode, Fv_ManholeOrder,
//  Fv_LinkOrder (he_topology.cpp, he_input.cpp:5680).
{
    int  i, j, k, l, m, N = Nobjects[NODE], L = Nobjects[LINK];
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

    // --- connecting links in link input order
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

    // --- split manholes
    for ( j = 0; j < N; j++ )
    {
        int nout = 0;
        for ( m = 0; m < IsNode[j].nconn; m++ )
            if ( IsPipe[IsNode[j].conn[m]].frm == j ) nout++;
        IsNode[j].split = nout > 1;
        if ( !IsNode[j].isOutlet && nout == 0 )
        {
            FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
            return netError(ERR_ISEWER_DEADEND, Node[j].ID);
        }
        if ( IsNode[j].isOutlet && nout > 0 )
        {
            FREE(indeg); FREE(vn); FREE(vl); FREE(cnt);
            return netError(ERR_OUTFALL, Node[j].ID);
        }
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
                        return netError(ERR_LOOP, Link[l].ID);
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
        return netError(ERR_LOOP, j < N ? Node[j].ID : "");
    }

    // --- Fv_LinkOrder he_topology.cpp:1525
    for ( i = 0; i < L; i++ ) vl[i] = 0;
    for ( m = 0; m < norder; m++ )
    {
        TIsNode* nd;
        k = IsOrder[m];
        nd = &IsNode[k];
        if ( nd->isOutlet ) continue;
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

static void steadySolve(void)
//  Fv_StaticRun he_staticrun.cpp:344, then Fv_GetPipeHydraulic,
//  Fv_GetCriticalDepth, Fv_GetPipeFlowDepth, Fv_StaticFilterResult, He_RunHead.
{
    int i, k, m, N = Nobjects[NODE], L = Nobjects[LINK];
    int *visit;
    double *Q;

    SQ  = (double*) calloc(L + 1, sizeof(double));
    SUP = (double*) calloc(L + 1, sizeof(double));
    SPK = (double*) calloc(L + 1, sizeof(double));
    SCV = (double*) calloc(L + 1, sizeof(double));
    SCC = (double*) calloc(L + 1, sizeof(double));
    SINF = (double*) calloc(L + 1, sizeof(double));
    SDepth = (double*) calloc(L + 1, sizeof(double));
    SVel = (double*) calloc(L + 1, sizeof(double));
    SCrit = (double*) calloc(L + 1, sizeof(double));
    SFull = (int*) calloc(L + 1, sizeof(int));
    SFtype = (int*) calloc(L + 1, sizeof(int));
    SPartial = (int*) calloc(L + 1, sizeof(int));
    SAdjFull = (int*) calloc(L + 1, sizeof(int));
    SStatic = (double*) calloc(N + 1, sizeof(double));
    SHead = (double*) calloc(N + 1, sizeof(double));
    SInc = (double*) calloc(N + 1, sizeof(double));
    visit = (int*) calloc(L + 1, sizeof(int));
    if ( !SQ || !SUP || !SPK || !SCV || !SCC || !SINF || !SDepth || !SVel ||
         !SCrit || !SFull || !SFtype || !SPartial || !SAdjFull || !SStatic ||
         !SHead || !SInc || !visit )
    {
        FREE(visit);
        report_writeErrorMsg(ERR_MEMORY, "");
        return;
    }
    Q = SQ;

    // --- flows, upstream to downstream
    for ( m = 0; m < IsNorder; m++ )
    {
        TIsNode* nd;
        double up, pk, cv, cc, inf, inflow;
        k = IsOrder[m];
        nd = &IsNode[k];
        if ( nd->nconn == 1 && !nd->isOutlet )
        {
            // Fv_GetSingleSourceNodeFlow he_staticrun.cpp:423
            int l = nd->conn[0];
            SStatic[k] = staticTotal(k);
            staticGetLoading(k, &SUP[l], &SPK[l], &SCV[l], &SCC[l], &SINF[l]);
            Q[l] = SUP[l] + peakFlow(SPK[l]) + coverFlow(SCV[l], SCC[l]) + 0.0;
            visit[l] = 1;
            continue;
        }

        // Fv_GetConnectNodeFlow he_staticrun.cpp:518
        up = pk = cv = cc = inf = 0.0;
        if ( !nd->isOutlet )
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
        if ( fabs(inflow) < IsFlowTol )
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
            Q[l] = inflow;  SUP[l] = up;  SPK[l] = pk;  SCV[l] = cv;  SCC[l] = cc;
            SINF[l] = inf;
            visit[l] = 1;
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

    // --- Fv_GetPipeHydraulic he_staticrun.cpp:681
    for ( i = 0; i < L; i++ )
        SFull[i] = (Q[i] > IsFlowTol && Q[i] >= IsPipe[i].fullFlow);

    // --- Fv_GetCriticalDepth he_staticrun.cpp:732
    for ( i = 0; i < L; i++ )
    {
        if ( fabs(Q[i]) < IsFlowTol ) SCrit[i] = 0.0;
        else
        {
            double q = MIN(Q[i], IsPipe[i].fullFlow) / IsPipe[i].count;
            SCrit[i] = isewer_criticalDepth(i, q, IsPipe[i].dia);
        }
    }

    // --- Fv_GetPipeFlowDepth he_staticrun.cpp:776
    for ( i = 0; i < L; i++ )
    {
        double q = Q[i] / IsPipe[i].count;
        SDepth[i] = SFull[i] ? IsPipe[i].dia : isewer_flowDepth(i, q);
        if ( SDepth[i] > 0.0 )
            SVel[i] = isewer_pipeVelocity(SDepth[i], 0.5 * IsPipe[i].dia,
                                          MIN(q, Q[i] / IsPipe[i].count));
        else SVel[i] = 0.0;
    }

    // --- Fv_StaticFilterResult he_staticrun.cpp:878
    for ( i = 0; i < L; i++ )
    {
        if ( !(SDepth[i] / IsPipe[i].dia >= IS_DD_TOL) )
        {
            SDepth[i] = 0.0;
            SVel[i] = 0.0;
        }
    }

    // --- heads
    isewer_runHead(Q, SDepth, SHead, SFtype, SPartial, SAdjFull, 1);
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

static const double* HQ;      // link flows (cfs)
static const double* HY;      // link depths (ft)
static double* HHyd;          // working node heads
static double* HOut;          // output node heads
static int*    HLevel;        // upstream pipes done per node
static int*    HFtype;
static int*    HPartial;
static int*    HAdjFull;
static int*    HSplitOrder;   // order in which split nodes were reached
static int*    HSplitNode;    // node for each split order
static int     HCounter;

static int fully(int k) { return HLevel[k] == IsNode[k].nuph; }

static void adjustNode(int k, double h)
//  Fv_AdjustNodeHead he_runhead.cpp:1910 (load manhole / flow split)
{
    TIsNode* nd = &IsNode[k];
    int m;
    double dh = -IS_BIG;
    // Fv_AdjustLoadNodeDownstream he_runhead.cpp:2152
    // QUIRK: h = max(h, dh) runs for every pipe, with dh kept from the last
    // flowing one.
    for ( m = 0; m < nd->ndown; m++ )
    {
        int l = nd->down[m];
        if ( HQ[l] > IS_FLOW_TOL ) dh = HY[l] + IsPipe[l].fromInv;
        h = MAX(h, dh);
    }
    // Fv_AdjustNodeMinorLoss he_runhead.cpp:2258: full-pipe velocity of the
    // total flow of each flowing downstream pipe; K from the SWMM inlet loss
    for ( m = 0; m < nd->ndown; m++ )
    {
        int l = nd->down[m];
        if ( HQ[l] > IS_FLOW_TOL )
        {
            double D = IsPipe[l].dia;
            double v = isewer_pipeVelocity(D, 0.5 * D, HQ[l]);
            h += IsPipe[l].kLoss * ((v * v) / (2.0 * 32.174));
        }
    }
    if ( nd->split ) HHyd[k] = MAX(HHyd[k], h);
    else             HHyd[k] = h;
    HOut[k] = MIN(HHyd[k], nd->cap);   // unlocked: cap = rim
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
        pipeHead(l, nk, HHyd[k]);               // Fv_GetNodeHead
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
        pipeHead(l, nk, HHyd[k]);
        HLevel[k]++;
        depthSearch(nk);
    }
}

void isewer_runHead(const double* Q, const double* y, double* head,
                    int* ftype, int* partial, int* adjFull, int resetFtype)
//  He_RunHead::F_Run he_runhead.cpp:554. Heads from each outlet upstream.
//  QUIRK: nested flow splits can leave a stale head upstream of a split.
//  QUIRK: in EPS the flow type flags are never reset between report times
//  (resetFtype = 0); a pipe the sweep does not reach keeps its old flag.
{
    int k, i, N = Nobjects[NODE], L = Nobjects[LINK];
    int nsplit = 0, remain, start, end, tmpStart;
    double ff;

    HQ = Q;  HY = y;  HOut = head;
    HFtype = ftype;  HPartial = partial;  HAdjFull = adjFull;
    HHyd = (double*) calloc(N + 1, sizeof(double));
    HLevel = (int*) calloc(N + 1, sizeof(int));
    HSplitOrder = (int*) calloc(N + 1, sizeof(int));
    HSplitNode = (int*) calloc(N + 1, sizeof(int));
    if ( !HHyd || !HLevel || !HSplitOrder || !HSplitNode )
    {
        FREE(HHyd); FREE(HLevel); FREE(HSplitOrder); FREE(HSplitNode);
        report_writeErrorMsg(ERR_MEMORY, "");
        return;
    }
    for ( k = 0; k < N; k++ )
    {
        HHyd[k] = -IS_BIG;
        HOut[k] = 0.0;            // never initialised in C++; unreached nodes
        HSplitOrder[k] = -1;
        if ( IsNode[k].split ) nsplit++;
    }
    for ( i = 0; i < L; i++ )
    {
        if ( resetFtype ) HFtype[i] = 0;
        HPartial[i] = 0;  HAdjFull[i] = 0;
    }
    HCounter = 0;

    for ( k = 0; k < N; k++ )
    {
        TIsNode* nd = &IsNode[k];
        int m;
        if ( !nd->isOutlet ) continue;
        // Fv_GetOutletHead he_runhead.cpp:1296
        ff = -IS_BIG;
        for ( m = 0; m < nd->nuph; m++ )
        {
            int l = nd->uph[m];
            ff = MAX(ff, HY[l] + IsPipe[l].toInv);
        }
        HHyd[k] = MAX(nd->submerge, ff);
        HOut[k] = HHyd[k];
        depthSearch(k);
    }

    // --- split nodes, level by level (QUIRK: remain counts levels)
    remain = nsplit;
    start = 0;
    while ( remain > 0 )
    {
        tmpStart = HCounter;
        end = HCounter;
        if ( start >= end ) break;
        for ( i = start; i < end; i++ ) splitDepthSearch(HSplitNode[i]);
        remain--;
        start = tmpStart;
    }
    FREE(HHyd); FREE(HLevel); FREE(HSplitOrder); FREE(HSplitNode);
}

//=============================================================================
//  SWMM state, loads and execution
//=============================================================================

static void setLinkState(int j, double q, double y)
{
    int k = Link[j].subIndex;
    double a = (y > 0.0) ? xsect_getAofY(&Link[j].xsect, y) : 0.0;
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
    Node[j].newVolume = 0.0;
    Node[j].overflow = 0.0;
}

static void setSteadyState(void)
//  Steady state: SWMM carries the InfoSewer solution unchanged at every step.
//  Each manhole receives a lateral inflow equal to InfoSewer's flow increment
//  there (outflow - inflow, negative where peaking lowers a combined flow),
//  so SWMM's continuity closes on the peaked flows.
{
    int i, j;
    for ( j = 0; j < Nobjects[NODE]; j++ ) SInc[j] = 0.0;
    for ( i = 0; i < Nobjects[LINK]; i++ )
    {
        double q = MAX(SQ[i], 0.0);
        setLinkState(i, q, SDepth[i]);
        if ( !IsNode[IsPipe[i].frm].isOutlet ) SInc[IsPipe[i].frm] += q;
        if ( !IsNode[IsPipe[i].to].isOutlet )  SInc[IsPipe[i].to]  -= q;
    }
    for ( j = 0; j < Nobjects[NODE]; j++ ) setNodeState(j, SHead[j]);
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
    (void)links;  (void)tStep;
    if ( !IsNode ) return 0;
    if ( IsAnalysis == IS_EPS )
    {
        // NewRoutingTime is already the end of this step (routing.c)
        isewer_epsSetState(NewRoutingTime / 1000.0);
        isewer_epsNodeFlows();
        return 1;
    }
    else
    {
        for ( i = 0; i < Nobjects[LINK]; i++ ) setLinkState(i, MAX(SQ[i], 0.0), SDepth[i]);
        for ( j = 0; j < Nobjects[NODE]; j++ ) setNodeState(j, SHead[j]);
    }
    for ( i = 0; i < Nobjects[LINK]; i++ )
    {
        Node[Link[i].node1].outflow += Link[i].newFlow;
        Node[Link[i].node2].inflow  += Link[i].newFlow;
    }
    return 1;
}

double isewer_getVelocity(int j, double flow, double depth)
//  Velocity as InfoSewer computes it (per barrel flow over the 3.1416 circle
//  area; HHU_GetPipeVelocity), used by link_getVelocity for SWMM's .out and
//  report statistics.
{
    int k = Link[j].subIndex;
    double r = 0.5 * Link[j].xsect.yFull;
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
{
    int m, upSuper = 0;
    for ( m = 0; m < IsNode[k].nup; m++ ) if ( fr[IsNode[k].up[m]] > 1 ) upSuper = 1;
    if ( !upSuper ) return 0;
    for ( m = 0; m < IsNode[k].ndown; m++ ) if ( fr[IsNode[k].down[m]] < 1 ) return 1;
    return 0;
}

static void writeSteadyReport(void)
{
    int i, k, L = Nobjects[LINK], N = Nobjects[NODE];
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
    for ( i = 0; i < L; i++ )
    {
        pipeOutputs(i, &o);
        fprintf(f, "\n  %-16s %10.4f %10.4f %10.4f %10.4f %10.4f %10.3f %9.3f %10.3f   %s",
            Link[i].ID, o.flow * ff, o.unpeak * ff, o.peak * ff, o.cover * ff,
            o.infil * ff, o.vel * lf, o.dD, o.qQ, o.flowtype > 0 ? "YES" : "no");
    }
    fprintf(f,
"\n\n  ************************************"
"\n  InfoSewer Steady-State Manhole Results"
"\n  ************************************\n");
    fprintf(f,
"\n  -------------------------------------------------------------------"
"\n                         Base Load    Total Load        HGL     Surch"
"\n  Manhole                %9s     %9s     %6s"
"\n  -------------------------------------------------------------------",
    FlowUnitWords[FlowUnits], FlowUnitWords[FlowUnits], UnitSystem == US ? "ft" : "m");
    for ( k = 0; k < N; k++ )
    {
        if ( IsNode[k].isOutlet ) continue;
        fprintf(f, "\n  %-16s %13.4f %13.4f %10.3f     %s", Node[k].ID,
            nodeBase(k) * ff, SStatic[k] * ff, SHead[k] * lf,
            SHead[k] >= IsNode[k].rim ? "YES" : "no");
    }
    fprintf(f, "\n");
}

static void writeSteadyCsv(void)
//  Every InfoSewer steady-state output array at full precision, in the
//  engine's internal units (cfs, ft), for verification.
{
    int i, k, L = Nobjects[LINK], N = Nobjects[NODE];
    double *fr;
    TIsPipeOut o;
    FILE* f = fopen(IsResultsFile, "wt");
    if ( !f ) { report_writeErrorMsg(ERR_RPT_FILE, IsResultsFile); return; }
    fr = (double*) calloc(L + 1, sizeof(double));
    fprintf(f, "kind,id,flow,unpeak,peak,cover,infil,storm,flowtype,full,rescap,vel,"
               "dD,qQ,depth,crit,froude,covcnt,backwater,adjdepth,adjvel,partial\n");
    for ( i = 0; i < L; i++ )
    {
        pipeOutputs(i, &o);
        if ( fr ) fr[i] = o.froude;
        fprintf(f, "pipe,%s,%.17g,%.17g,%.17g,%.17g,%.17g,0,%.17g,%.17g,%.17g,%.17g,"
                   "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            Link[i].ID, o.flow, o.unpeak, o.peak, o.cover, o.infil, o.flowtype,
            o.full, o.rescap, o.vel, o.dD, o.qQ, o.depth, o.crit, o.froude,
            o.covcnt, o.backwater, o.adjdepth, o.adjvel, o.partial);
    }
    fprintf(f, "kind,id,base,static,nstorm,head,surch,jump\n");
    for ( k = 0; k < N; k++ )
    {
        if ( IsNode[k].isOutlet ) continue;
        fprintf(f, "node,%s,%.17g,%.17g,0,%.17g,%.17g,%.17g\n", Node[k].ID,
            nodeBase(k), SStatic[k], SHead[k],
            SHead[k] >= IsNode[k].rim ? 1.0 : 0.0,
            fr ? (double)nodeJump(k, fr) : 0.0);
    }
    fprintf(f, "kind,id,oflow,ohead\n");
    for ( k = 0; k < N; k++ )
    {
        double q = 0.0;
        if ( !IsNode[k].isOutlet ) continue;
        // Fv_GetOutletFlow: sums the UNclipped link flows
        for ( i = 0; i < L; i++ ) if ( IsPipe[i].to == k ) q += SQ[i];
        fprintf(f, "outlet,%s,%.17g,%.17g\n", Node[k].ID, q, SHead[k]);
    }
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
    FREE(IsNode); FREE(IsPipe); FREE(IsOrder); FREE(IsAllocList); FREE(IsRim);
    for ( p = 0; p < IsNpats; p++ ) FREE(IsPats[p].f);
    for ( p = 0; p < IsNcurves; p++ ) { FREE(IsCurves[p].x); FREE(IsCurves[p].y); }
    FREE(IsPats); FREE(IsCurves); FREE(IsLoads); FREE(SplitIn); FREE(InfilIn);
    IsNpats = MaxPats = IsNcurves = MaxCurves = IsNloads = MaxLoads = 0;
    NsplitIn = MaxSplitIn = NinfilIn = MaxInfilIn = 0;
}

void isewer_close(void)
{
    if ( IsAnalysis == IS_EPS ) isewer_epsClose();
    freeAll();
}
