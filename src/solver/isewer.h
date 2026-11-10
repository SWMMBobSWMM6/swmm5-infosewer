//-----------------------------------------------------------------------------
//   isewer.h
//
//   Project:  swmm5-infosewer (EPA SWMM 5.2.4 + InfoSewer routing)
//   Purpose:  shared definitions for isewer.c (input, network, steady state,
//             HGL) and isewer_eps.c (extended period routing).
//
//   Written by Claude (AI) for SWMM Enablement, November 2026.
//-----------------------------------------------------------------------------
#ifndef ISEWER_H
#define ISEWER_H

//-----------------------------------------------------------------------------
//  InfoSewer engine constants (he_define.h, stdafx.h). Copied exactly: the
//  normal-depth bisection stops on a 1e-4 relative flow error, so changing
//  1.49, 0.6667 or 3.1416 moves depths by up to ~1e-4.
//-----------------------------------------------------------------------------
#define IS_ZERO_TOL       0.0001       // HMAP_ZeroTolerance
#define IS_DEPTH_TOL      0.0000001    // HMAP_DepthTolerance
#define IS_ITER_FLOW_TOL  0.00001      // HMAP_IterFlowTolerance
#define IS_MAX_MANNING    1.0          // HMAP_MaxManningCoeff
#define IS_WETAREA_TOL    0.00001      // HMAP_WetAreaTolerance
#define IS_DD_TOL         0.01         // HMAP_DepthDiamTolerance
#define IS_FLOW_TOL       0.0001       // HMAP_FlowTolerance (compile-time)
#define IS_SYS_ZERO       1.0e-10      // SYS_Zero
#define IS_SYS_ZERO_FLOAT 1.0e-7       // SYS_ZeroFloat
#define IS_PI_DEF         3.1415926    // 'pi' macro
#define IS_BIG            1.0e37
#define IS_M3_PER_FT3     0.028317     // M3perFT3 (not 0.3048^3)

enum IsAnalysisType { IS_STEADY, IS_EPS };
enum IsLoadType     { IS_UNPEAK, IS_PEAK, IS_COVER };
enum IsSplitType    { IS_SPLIT_HYD, IS_SPLIT_PCT };
enum IsInfilType    { IS_INF_NONE, IS_INF_LEN, IS_INF_LENDIA, IS_INF_SURF,
                      IS_INF_COUNT, IS_INF_PATTERN };
enum IsSlopeType    { IS_SLOPE_STANDARD, IS_SLOPE_FLAT, IS_SLOPE_ADVERSE };
enum IsNodeKind     { IS_NODE_LOAD, IS_NODE_CHAMBER, IS_NODE_OUTLET, IS_NODE_WETWELL };
enum IsLinkKind     { IS_LINK_PIPE, IS_LINK_FM, IS_LINK_STATION };
enum IsPumpType     { IS_PUMP_CAPACITY, IS_PUMP_POINT, IS_PUMP_3POINT };
enum IsSiphonType   { IS_SIPHON_STANDARD, IS_SIPHON_STATION };

typedef struct
{
    int    node;             // SWMM node index
    double baseIn;           // base load as read (project flow units)
    double base;             // base load (cfs) = baseIn / flow factor
    int    type;             // IS_UNPEAK, IS_PEAK, IS_COVER
    char   pattern[MAXMSG+1];
    int    pat;              // pattern index, -1 = none
    double cover;            // coverage unit count
}  TIsLoad;

typedef struct
{
    char    id[MAXMSG+1];
    int     linear;          // 1 = linear interpolation, 0 = stepwise
    double* f;
    int     n, nmax;
}  TIsPattern;

typedef struct
{
    char    id[MAXMSG+1];
    double* x;
    double* y;
    int     n, nmax;
}  TIsCurve;

typedef struct               // one InfoSewer topology link: gravity pipe, force
{                            // main (SWMM conduits) or pump station (virtual)
    int    kind;             // IS_LINK_PIPE, IS_LINK_FM, IS_LINK_STATION
    int    swmm;             // SWMM link index (-1 for a station)
    int    frm, to;          // internal node indexes
    double fromInv, toInv;   // invert elevations (ft)
    double length;           // ft
    double dia;              // ft
    double n;                // Manning n
    int    count;            // parallel barrels
    double slope;
    int    slopeType;
    int    adverse;          // from invert below to invert
    int    split;            // IS_SPLIT_HYD or IS_SPLIT_PCT
    double splitPct;
    int    infType;
    double infRate;          // converted (cfs per ft, per ft2, ...)
    double infCount;
    int    infPat;           // pattern for IS_INF_PATTERN
    double kLoss;            // manhole minor loss K (SWMM inlet loss)
    double fullFlow;         // full flow x count (cfs); -count if adverse
    double friction;         // force main Hazen-Williams C
    int    siphon;           // force main: siphon index
    int    siphonType;       // force main: IS_SIPHON_STANDARD / _STATION
    int    station;          // station link: station index
}  TIsPipe;

typedef struct               // one InfoSewer pump (= SWMM pump)
{
    int    swmm;             // SWMM link index
    int    frm, to;          // internal node indexes (wet well, chamber)
    int    wetwell;          // wet well index
    int    station;          // station index
    int    type;             // IS_PUMP_CAPACITY, _POINT, _3POINT
    int    parallel;         // pumps in parallel (extra + 1)
    double capacity, head0, head1, flow1, head2, flow2;   // cfs, ft
    double speed;
    int    pattern;
    double staticFlow, baseFlow;
    double a, b;             // head = head0 - a*Q^b (parallel embedded)
}  TIsPump;

typedef struct               // one InfoSewer wet well (= SWMM storage node)
{
    int    node;             // internal node index
    int    curve;            // level-volume curve (-1 = cylinder)
    double bottom, minLevel, maxLevel, initLevel, diameter;
    double minVolume, maxVolume, initVolume;
}  TIsWetWell;

typedef struct               // pumps between one wet well and one chamber
{
    int    npump;
    int*   pump;             // pump indexes, input order
    int    link;             // topology link index
    int    wetwell;
    int    frm, to;          // internal node indexes
    int    siphon;           // siphon fed by the station
}  TIsStation;

typedef struct               // chain of force mains
{
    int    nfm;
    int*   fm;               // internal link indexes, upstream first
    int    type;             // IS_SIPHON_STANDARD / IS_SIPHON_STATION
    int    station;
}  TIsSiphon;

typedef struct               // one InfoSewer node (SWMM junction, outfall, storage)
{
    int    kind;             // IS_NODE_LOAD, _CHAMBER, _OUTLET, _WETWELL
    int    swmm;             // SWMM node index
    int    wetwell;          // wet well index
    int    isOutlet;
    double rim;              // ft
    double cap;              // head cap: rim + surcharge depth (ft)
    double bottomRef;        // lowest connected invert (ft)
    double submerge;         // outlet fixed head (ft) or -IS_BIG
    int    split;            // more than one outgoing pipe
    int    nconn, *conn;     // connecting pipes, in link order
    int    ndown, *down;     // downstream pipes (Fv_LinkOrder)
    int    nup,   *up;       // upstream pipes  (Fv_LinkOrder)
    int    nuph,  *uph;      // pipes used by the HGL sweep (conn - down)
    int    nload, *load;     // indexes into IsLoads, in input order
}  TIsNode;

//-----------------------------------------------------------------------------
//  Shared data (defined in isewer.c)
//-----------------------------------------------------------------------------
extern int        IsAnalysis;
extern double     IsLoadMult;
extern double     IsMinSegment, IsMinTravelTime;
extern int        IsMaxSegments, IsFlowAttenuation;
extern double     IsFlowTol, IsIterTol;     // converted to cfs
extern double     IsFlowFactor;             // InfoSewer flow unit factor
extern double     IsPatternStep;            // sec
extern double     IsPumpStep;               // sec (pump station time step)
extern int        IsLinkReport;             // InfoSewer [Report] Link : Yes
extern TIsNode*   IsNode;
extern TIsPipe*   IsPipe;
extern TIsPump*   IsPump;
extern TIsWetWell* IsWell;
extern TIsStation* IsStation;
extern TIsSiphon* IsSiphon;
extern int        IsNN;            // nodes
extern int        IsNP, IsNF;      // gravity pipes, force mains
extern int        IsNPM;           // pipe mains = IsNP + IsNF
extern int        IsNPump, IsNWell, IsNStation, IsNSiphon;
extern int        IsNT;            // topology links = IsNPM + IsNStation
extern int        IsNLE;           // link elements = IsNPM + max(IsNPump, IsNStation)
extern int*       IsNodeOfSwmm;    // SWMM node -> internal
extern int*       IsLinkOfSwmm;    // SWMM link -> internal (pumps: IsNPM + pump)
extern int*       IsOrder;                  // nodes, upstream to downstream
extern TIsLoad*   IsLoads;
extern int        IsNloads;
extern TIsPattern* IsPats;
extern int        IsNpats;
extern char       IsResultsFile[];

//-----------------------------------------------------------------------------
//  Shared functions
//-----------------------------------------------------------------------------
// circular hydraulics (InfoSewer formulas)
double isewer_wetArea(double d, double r);
double isewer_wetPerimeter(double d, double r);
double isewer_wetWidth(double d, double r);
double isewer_linkFlow(int i, double d, double r, double n);
double isewer_fullFlowBarrel(int i);
double isewer_flowDepth(int i, double q);
double isewer_pipeVelocity(double d, double r, double q);
double isewer_headloss(int i, double q);
double isewer_criticalDepth(int i, double q, double D);
double isewer_infiltration(int i);

// wet wells (he_util.cpp)
double isewer_wellVolumeFromLevel(int w, double level);
double isewer_wellLevelFromVolume(int w, double volume);

// force mains and siphons (HD_ForceMainData, HD_SiphonData)
double isewer_siphonHeadloss(int s, double q);
double isewer_siphonTailHead(int s);
double isewer_siphonNewtonFactor(int s, double q);
void   isewer_fmHydraulics(int l, double q, double* velocity, double* headloss);

// HGL sweep (He_RunHead): flows Q[], depths y[] in; heads, flags out
int    isewer_runHead(const double* Q, const double* y, double* head,
                      int* ftype, int* partial, int* adjFull, int resetFtype);
void   isewer_setHeadPressureData(const double* wellHead, const double* pumpHead);

// EPS (isewer_eps.c)
int    isewer_epsInit(void);
void   isewer_epsSetState(double tSec);
void   isewer_epsNodeFlows(void);
void   isewer_epsLoads(double tSec);
void   isewer_epsReport(void);
void   isewer_epsClose(void);

#endif
