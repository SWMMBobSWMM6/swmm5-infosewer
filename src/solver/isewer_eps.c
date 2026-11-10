//-----------------------------------------------------------------------------
//   isewer_eps.c
//
//   Project:  swmm5-infosewer (EPA SWMM 5.2.4 + InfoSewer gravity routing)
//   Purpose:  InfoSewer extended period (EPS, "dynamic") gravity routing.
//
//   InfoSewer EPS is not a network time-stepper: each pipe is routed for the
//   whole run, upstream first, on its own sub-step (spec/eps.md section 1).
//   So the whole run is routed at start-up (isewer_epsInit) and SWMM's routing
//   steps then read the stored report-time results: exactly at report times,
//   linearly interpolated between them.
//
//   Engine functions reproduced (he_dynamicrun.cpp unless noted):
//     Fv_ManholePipeInit, Fv_TimeLoopAttenuate, Fv_SegmentOutflow,
//     Fv_SegmentFlowAttenuate, Fv_DownstreamFlow, Fv_ReportPipeData,
//     Fv_CalPipeReportData, Fv_OutputPipeFlowRange, Fv_SplitFlow,
//     Fv_MultiHydSplitFlow, Fv_SinglePrecedeFlow, HU_GetPatternItemValue
//     (he_util.cpp), Fv_DynamicData / Fv_PipeBackwaterAdjustment (outputs),
//     He_RunHead (isewer.c).
//   Quirks are reproduced on purpose and marked QUIRK.
//
//   Written by Claude (AI) for SWMM Enablement, November 2026.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"
#include "isewer.h"

#define DOUBLE_TIME_TOL  0.0001    // HMAP_DoubleTimeTolerance
#define ZERO_TIME_VALUE  0.001     // he_dynamicrun.cpp:30

//-----------------------------------------------------------------------------
//  Time base (spec 2.1)
//-----------------------------------------------------------------------------
static double StartT;          // InfoSewer simulation start time (sec of day)
static double SimDur;          // duration (sec)
static double EndT;            // StartT + SimDur
static double RptStep;         // report time step (sec)
static double PatStep;         // pattern time step after the split (sec)
static int    PatInc;          // pattern split factor
static double H;               // hydraulic step = PatStep
static int    NHyd;            // number of H steps
static int    NRpt;            // reportStepCount

//-----------------------------------------------------------------------------
//  Expanded patterns (values repeated PatInc times)
//-----------------------------------------------------------------------------
static double** PatV;
static int*     PatN;

//-----------------------------------------------------------------------------
//  Per pipe routing results
//-----------------------------------------------------------------------------
typedef struct
{
    double qUp, qDn, qIn, vel, dep;    // raw report values (spec 6.1); qIn = true inlet flow
}  TRptRaw;

static double**  DF;           // outflow hydrograph per pipe at H steps [NHyd+1]
static TRptRaw** Raw;          // per pipe, per report index [NRpt+2]
static int*      NRaw;
static double*   MaxQ;         // max flow over sub-steps
static double*   MaxT;         // its time (sec, absolute)
static int*      SegN;         // segments
static double*   SegDt;        // sub-step
static double    MbIn, MbOut;  // InfoSewer mass balance terms (ft3)

// work arrays for one pipe
static double *WFormer, *WCurrent, *WInflow, *WOutflow, *WIn0, *WSimT, *WHydT;
static int     WSize;

//-----------------------------------------------------------------------------
//  Report-time outputs (spec 6) [report][pipe or node]
//-----------------------------------------------------------------------------
static double **OFlow, **OVel, **ODepth, **OdD, **OqQ, **OFroude, **OBack,
              **OAdjD, **OAdjV, **OFtype, **OPartial, **OUp, **ODn, **OIn;
static double **OHead, **OLoad, **OStatus, **OJump, **OOutflow;
static double  PipeStorage;    // final pipe storage (ft3)

// current SWMM state
static double *CurQin, *CurQout;  // conduit inlet and outlet flows (cfs)
static double *CurVol;            // water in transit in each conduit (ft3)
static double  CurT;              // time of the current state (sec)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int    allocAll(void);
static double patternValue(int p, double simTime);
static double nodeLoad(int k, double t);
static double nodeLoadTotal(int k, double t);
static double infilFlow(int i, double t);
static double interpHydro(const double* arr, double hydTime);
static int    splitFlows(int k, double** out);
static void   multiHydSplit(int* pipes, int np, double q, double* qout);
static void   pipeInit(int i);
static int    routePipe(int k, int i, int downIdx, double** splitArr);
static void   reportData(int i, double qUp, double* vel, double* dep);
static void   mcCoeffs(int i, double flow0, double deltaX, double* c1,
                       double* c2, double* c3);
static int    computeOutputs(void);
static void   writeCsv(void);
static double** alloc2(int n, int m);
static void   free2(double** a, int n);

//=============================================================================

int isewer_epsInit(void)
//  Routes the whole run (spec 1, eps_ref.py:run) and builds the report-time
//  outputs.
{
    int m, k, i;
    double tod;

    // --- time base (he_input.cpp:820-993)
    tod = StartDateTime - floor(StartDateTime);
    StartT = floor(tod * SECperDAY + 0.5);
    SimDur = floor((EndDateTime - StartDateTime) * SECperDAY + 0.5);
    RptStep = (double)ReportStep;
    PatStep = IsPatternStep;
    PatInc = 1;
    if ( PatStep > RptStep + DOUBLE_TIME_TOL )
    {
        PatInc = (int)(PatStep / RptStep + 0.9995);
        if ( PatInc < 1 ) PatInc = 1;
        PatStep /= PatInc;
    }
    EndT = StartT + SimDur;
    NRpt = (int)(SimDur / RptStep + 1.0 + DOUBLE_TIME_TOL);
    H = PatStep;
    NHyd = (int)(SimDur / H + 1.0 + DOUBLE_TIME_TOL);
    if ( IsMinTravelTime < 1.0 ) IsMinTravelTime = 1.0;     // he_input.cpp:10554
    if ( IsMaxSegments < 1 ) IsMaxSegments = 1000;          // he_input.cpp:10597

    if ( allocAll() ) return ErrorCode;

    // --- route pipes, upstream first
    MbIn = MbOut = 0.0;
    for ( m = 0; m < Nobjects[NODE]; m++ )
    {
        double** splitArr = NULL;
        TIsNode* nd;
        k = IsOrder[m];
        nd = &IsNode[k];
        if ( nd->isOutlet ) continue;
        if ( nd->ndown > 1 && nd->split )
        {
            splitArr = alloc2(nd->ndown, NHyd + 2);
            if ( !splitArr ) { report_writeErrorMsg(ERR_MEMORY, ""); return ErrorCode; }
            splitFlows(k, splitArr);
        }
        for ( i = 0; i < nd->ndown; i++ )
        {
            if ( routePipe(k, nd->down[i], i, splitArr) )
            {
                free2(splitArr, nd->ndown);
                return ErrorCode;
            }
        }
        free2(splitArr, nd->ndown);
    }
    FREE(WFormer); FREE(WCurrent); FREE(WInflow); FREE(WOutflow); FREE(WIn0);
    FREE(WSimT); FREE(WHydT);
    WSize = 0;
    if ( computeOutputs() ) return ErrorCode;
    if ( IsResultsFile[0] ) writeCsv();
    return 0;
}

static double** alloc2(int n, int m)
{
    int i;
    double** a = (double**) calloc(n + 1, sizeof(double*));
    if ( !a ) return NULL;
    for ( i = 0; i < n; i++ )
    {
        a[i] = (double*) calloc(m + 1, sizeof(double));
        if ( !a[i] ) { free2(a, i); return NULL; }
    }
    return a;
}

static void free2(double** a, int n)
{
    int i;
    if ( !a ) return;
    for ( i = 0; i < n; i++ ) FREE(a[i]);
    free(a);
}

static int allocAll(void)
{
    int i, p, L = Nobjects[LINK], N = Nobjects[NODE];

    // --- expanded patterns (Fv_ReadPattern: each value stored PatInc times)
    PatV = (double**) calloc(IsNpats + 1, sizeof(double*));
    PatN = (int*) calloc(IsNpats + 1, sizeof(int));
    if ( !PatV || !PatN ) goto nomem;
    for ( p = 0; p < IsNpats; p++ )
    {
        int j, r;
        PatN[p] = IsPats[p].n * PatInc;
        PatV[p] = (double*) calloc(PatN[p] + 1, sizeof(double));
        if ( !PatV[p] ) goto nomem;
        for ( j = 0; j < IsPats[p].n; j++ )
            for ( r = 0; r < PatInc; r++ ) PatV[p][j * PatInc + r] = IsPats[p].f[j];
    }

    DF = alloc2(L, NHyd + 2);
    Raw = (TRptRaw**) calloc(L + 1, sizeof(TRptRaw*));
    NRaw = (int*) calloc(L + 1, sizeof(int));
    MaxQ = (double*) calloc(L + 1, sizeof(double));
    MaxT = (double*) calloc(L + 1, sizeof(double));
    SegN = (int*) calloc(L + 1, sizeof(int));
    SegDt = (double*) calloc(L + 1, sizeof(double));
    CurQin = (double*) calloc(L + 1, sizeof(double));
    CurQout = (double*) calloc(L + 1, sizeof(double));
    CurVol = (double*) calloc(L + 1, sizeof(double));
    CurT = -1.0;
    if ( !DF || !Raw || !NRaw || !MaxQ || !MaxT || !SegN || !SegDt || !CurQin ||
         !CurQout || !CurVol ) goto nomem;
    for ( i = 0; i < L; i++ )
    {
        Raw[i] = (TRptRaw*) calloc(NRpt + 3, sizeof(TRptRaw));
        if ( !Raw[i] ) goto nomem;
    }
    OFlow = alloc2(NRpt, L);  OVel = alloc2(NRpt, L);  ODepth = alloc2(NRpt, L);
    OdD = alloc2(NRpt, L);    OqQ = alloc2(NRpt, L);   OFroude = alloc2(NRpt, L);
    OBack = alloc2(NRpt, L);  OAdjD = alloc2(NRpt, L); OAdjV = alloc2(NRpt, L);
    OFtype = alloc2(NRpt, L); OPartial = alloc2(NRpt, L);
    OUp = alloc2(NRpt, L);    ODn = alloc2(NRpt, L);   OIn = alloc2(NRpt, L);
    OHead = alloc2(NRpt, N);  OLoad = alloc2(NRpt, N); OStatus = alloc2(NRpt, N);
    OJump = alloc2(NRpt, N);  OOutflow = alloc2(NRpt, N);
    if ( !OFlow || !OVel || !ODepth || !OdD || !OqQ || !OFroude || !OBack ||
         !OAdjD || !OAdjV || !OFtype || !OPartial || !OUp || !ODn || !OIn ||
         !OHead || !OLoad || !OStatus || !OJump || !OOutflow ) goto nomem;
    return 0;
nomem:
    report_writeErrorMsg(ERR_MEMORY, "");
    return ErrorCode ? ErrorCode : ERR_MEMORY;
}

//=============================================================================
//  Loads (spec 3)
//=============================================================================

static double patternValue(int p, double simTime)
//  HU_GetPatternItemValue he_util.cpp:693 for loading patterns.
{
    double *v, ps, value1, value2, deltaTime;
    int size, item, inc, index1, index2, ref;
    if ( p < 0 ) return 1.0;                 // flat default pattern
    v = PatV[p];  size = PatN[p];
    if ( size <= 0 ) return 1.0;
    if ( simTime < 0.0 ) return 0.0;
    item = (int)(simTime / PatStep + DOUBLE_TIME_TOL);
    while ( item >= size ) item -= size;
    if ( item < 0 ) item = 0;
    if ( !IsPats[p].linear ) return v[item];
    inc = PatInc;  ps = PatStep;
    index1 = (int)(simTime / (inc * ps)) * inc;
    while ( index1 >= size ) index1 -= size;
    index2 = index1 + 1;
    if ( index2 >= size ) index2 = 0;
    value1 = v[index1];  value2 = v[index2];
    if ( index2 * ps < simTime )
    {
        ref = MIN(index1, index2) + size;
        while ( (simTime - ref * ps) > (ps * inc) ) ref += size;
        if ( (size - index1) <= 1 ) ref -= 1;
    }
    else ref = index1;
    deltaTime = simTime - ref * ps;
    return value1 + deltaTime * (value2 - value1) / (ps * inc);
}

static double nodeLoad(int k, double t)
//  HD_LoadNodeData::F_GetDynamicTimeFlow: base * (mult * pattern).
//  EPS ignores load type, coverage count and peaking (he_data.cpp:2126).
{
    int m;
    double q = 0.0;
    for ( m = 0; m < IsNode[k].nload; m++ )
    {
        TIsLoad* ld = &IsLoads[IsNode[k].load[m]];
        q += ld->base * (IsLoadMult * patternValue(ld->pat, t));
    }
    return q;
}

static double nodeLoadTotal(int k, double t)
//  HD_LoadNodeData::F_GetDynamicTotalFlow: pattern * mult * base.
{
    int m;
    double q = 0.0;
    for ( m = 0; m < IsNode[k].nload; m++ )
    {
        TIsLoad* ld = &IsLoads[IsNode[k].load[m]];
        q += patternValue(ld->pat, t) * IsLoadMult * ld->base;
    }
    return q;
}

static double infilFlow(int i, double t)
//  HD_PipeData::F_GetInfiltrateFlow he_data.cpp:5420
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
      case IS_INF_PATTERN: f = p->infRate * patternValue(p->infPat, t); break;
      default:             f = 0.0;
    }
    return f * p->count;
}

static double interpHydro(const double* arr, double hydTime)
//  Fv_SinglePrecedeFlow / Fv_SplitPrecedeFlow interpolation (5831-5851).
//  QUIRK: extrapolates past the last H; the 0.001 offset uses the next
//  interval's slope just before a boundary.
{
    int item = (int)(hydTime / H + ZERO_TIME_VALUE);
    double deltaTime;
    if ( item > NHyd - 2 ) item = NHyd - 2;
    if ( item < 0 ) item = 0;
    deltaTime = hydTime - item * H;
    if ( fabs(deltaTime) <= ZERO_TIME_VALUE ) deltaTime = 0.0;
    return arr[item] + deltaTime * (arr[item + 1] - arr[item]) / H;
}

//=============================================================================
//  Flow splits (spec 4.5)
//=============================================================================

static int splitFlows(int k, double** out)
//  Fv_SplitFlow he_dynamicrun.cpp:2010: the split of every H step.
{
    TIsNode* nd = &IsNode[k];
    int *pct, *hyd, *isHyd, npct = 0, nhyd = 0, i, j, m, ti;
    double *qp, *qh, simTime, hydTime, total, inf, s;

    pct = (int*) calloc(nd->ndown + 1, sizeof(int));
    hyd = (int*) calloc(nd->ndown + 1, sizeof(int));
    isHyd = (int*) calloc(nd->ndown + 1, sizeof(int));
    qp = (double*) calloc(nd->ndown + 1, sizeof(double));
    qh = (double*) calloc(nd->ndown + 1, sizeof(double));
    if ( !pct || !hyd || !isHyd || !qp || !qh )
    {
        FREE(pct); FREE(hyd); FREE(isHyd); FREE(qp); FREE(qh);
        report_writeErrorMsg(ERR_MEMORY, "");
        return ErrorCode;
    }
    for ( m = 0; m < nd->ndown; m++ )
    {
        int l = nd->down[m];
        if ( IsPipe[l].split != IS_SPLIT_HYD ) { pct[npct++] = l;  isHyd[m] = 0; }
        else                                   { hyd[nhyd++] = l;  isHyd[m] = 1; }
    }
    ti = 0;
    simTime = StartT;
    hydTime = 0.0;
    do
    {
        total = nodeLoadTotal(k, simTime);
        MbIn += total * H;
        for ( m = 0; m < nd->nuph; m++ ) total += DF[nd->uph[m]][ti];
        // Fv_SplitInfiltrateFlow: percentage pipes first, then hydraulic
        inf = 0.0;
        for ( m = 0; m < npct; m++ ) inf += infilFlow(pct[m], simTime);
        for ( m = 0; m < nhyd; m++ ) inf += infilFlow(hyd[m], simTime);
        total += inf;
        MbIn += inf * H;

        // percentage pipes (not verified against the engine)
        s = 0.0;
        for ( m = 0; m < npct; m++ )
        {
            double percent = IsPipe[pct[m]].splitPct;
            if ( percent < 0.0 ) percent = 0.0;
            qp[m] = 0.01 * total * percent;
            s += qp[m];
        }
        if ( nhyd <= 0 && fabs(s - total) > 1.0e-5 )
            for ( m = 0; m < npct; m++ ) qp[m] = total * (qp[m] / s);

        // hydraulic pipes
        if ( nhyd == 1 ) qh[0] = total - s;
        else if ( nhyd > 1 ) multiHydSplit(hyd, nhyd, total - s, qh);

        for ( m = 0, i = 0, j = 0; m < nd->ndown; m++ )
        {
            if ( isHyd[m] ) out[m][ti] = qh[i++];
            else            out[m][ti] = qp[j++];
        }
        simTime += H;
        hydTime += H;
        ti++;
    } while ( simTime <= EndT && ti <= NHyd );
    FREE(pct); FREE(hyd); FREE(isHyd); FREE(qp); FREE(qh);
    return 0;
}

static void multiHydSplit(int* pipes, int np, double hydCalcFlow, double* q)
//  Fv_MultiHydSplitFlow he_dynamicrun.cpp:2343
{
    int i, allFull;
    double capacity = 0.0, minInvert, maxCrown, waterLevel, outflow, minLevel,
           maxLevel, depth;
    long iter = 0;

    for ( i = 0; i < np; i++ ) q[i] = 0.0;
    if ( hydCalcFlow < IsFlowTol ) return;
    for ( i = 0; i < np; i++ ) capacity += IsPipe[pipes[i]].fullFlow;
    if ( hydCalcFlow >= capacity )
    {
        for ( i = 0; i < np; i++ )
            q[i] = (IsPipe[pipes[i]].fullFlow / capacity) * hydCalcFlow;
        return;
    }
    minInvert = IsPipe[pipes[0]].fromInv;
    maxCrown = minInvert + IsPipe[pipes[0]].dia;
    for ( i = 1; i < np; i++ )
    {
        minInvert = MIN(minInvert, IsPipe[pipes[i]].fromInv);
        maxCrown = MAX(maxCrown, IsPipe[pipes[i]].fromInv + IsPipe[pipes[i]].dia);
    }
    waterLevel = 0.5 * (minInvert + maxCrown);
    outflow = 0.0;
    for ( i = 0; i < np; i++ )
    {
        TIsPipe* p = &IsPipe[pipes[i]];
        depth = waterLevel - p->fromInv;
        if ( depth <= 0.0 ) q[i] = 0.0;
        else
        {
            depth = MIN(depth, p->dia);
            q[i] = p->count * isewer_linkFlow(pipes[i], depth, -1.0, -1.0);
        }
        outflow += q[i];
    }
    minLevel = minInvert;  maxLevel = maxCrown;
    if ( outflow > hydCalcFlow ) maxLevel = waterLevel;
    if ( outflow < hydCalcFlow ) minLevel = waterLevel;
    waterLevel = 0.5 * (maxLevel + minLevel);
    while ( fabs(outflow - hydCalcFlow) / outflow > IsIterTol )
    {
        outflow = 0.0;
        allFull = 1;
        for ( i = 0; i < np; i++ )
        {
            TIsPipe* p = &IsPipe[pipes[i]];
            depth = waterLevel - p->fromInv;
            if ( depth <= 0.0 ) { q[i] = 0.0;  allFull = 0; }
            else
            {
                q[i] = p->count * isewer_linkFlow(pipes[i], MIN(depth, p->dia), -1.0, -1.0);
                if ( depth < (p->dia - IS_ZERO_TOL) ) allFull = 0;
            }
            outflow += q[i];
        }
        if ( allFull )
        {
            for ( i = 0; i < np; i++ ) q[i] = IsPipe[pipes[i]].fullFlow;
            return;
        }
        if ( outflow > hydCalcFlow ) maxLevel = waterLevel;
        if ( outflow < hydCalcFlow ) minLevel = waterLevel;
        waterLevel = 0.5 * (maxLevel + minLevel);
        if ( ++iter > 100000 ) break;            // safety only
    }
}

//=============================================================================
//  Pipe routing (spec 2.2-2.5, 4)
//=============================================================================

static void pipeInit(int i)
//  Fv_ManholePipeInit he_dynamicrun.cpp:4946
{
    TIsPipe* p = &IsPipe[i];
    int segCount;
    double length = p->length, radius, linkArea, deltaX, dt, multiple;

    segCount = (int)(length / IsMinSegment - 0.0001) + 1;
    if ( segCount > IsMaxSegments ) segCount = IsMaxSegments;
    radius = 0.5 * p->dia;
    linkArea = p->count * IS_PI_DEF * radius * radius;
    if ( length > IsMinSegment ) deltaX = length / segCount;
    else { segCount = 1;  deltaX = length; }
    if ( !(fabs(p->slope) > 1.0e-5) ) dt = H * 1.0;          // HU_IsFlatGravityMain
    else dt = (double)(int)((linkArea * deltaX) / p->fullFlow);
    if ( dt > H * 1.0 ) dt = H * 1.0;
    if ( dt < IsMinTravelTime ) dt = IsMinTravelTime;
    multiple = (double)((int)(H / dt) + 1);
    dt = H / multiple;
    if ( dt < 1.0 ) dt = 1.0;
    SegN[i] = segCount;
    SegDt[i] = dt;
    // deltaX kept by the caller through SegDx
}

static void mcCoeffs(int i, double flow0, double deltaX, double* c1, double* c2,
                     double* c3)
//  Fv_SegmentFlowAttenuate coefficient block (he_dynamicrun.cpp:6030-6091).
//  QUIRK: reference flow 6.125 x inflow; Courant number always 1.25.
{
    TIsPipe* p = &IsPipe[i];
    double D = p->dia, n = p->n, flow, fullFlow, depth, radius, area, top,
           celerity, unitq, dT, C, Dn;
    flow = (flow0 < 0.001) ? 0.001 : flow0;
    flow = 3.5 * flow;
    flow = 0.5 * (flow + 2.5 * flow);
    fullFlow = isewer_fullFlowBarrel(i);
    if ( fullFlow < IsFlowTol ) fullFlow = flow;
    if ( (flow >= fullFlow) || p->slopeType != IS_SLOPE_STANDARD )
    {
        flow = 0.25 * fullFlow;
        depth = isewer_flowDepth(i, flow);
        if ( depth >= D )
        {
            depth = 0.25 * D;
            flow = isewer_linkFlow(i, depth, 0.5 * D, n);
        }
    }
    else depth = isewer_flowDepth(i, flow);
    radius = 0.5 * D;
    area = isewer_wetArea(depth, radius);
    top = isewer_wetWidth(depth, radius);
    celerity = 1.25 * (flow / area);
    unitq = flow / top;
    dT = deltaX / (flow / area);
    C = celerity * (dT / deltaX);
    if ( C > 1.25 )
    {
        C = 1.25;
        celerity = C / (dT / deltaX);
    }
    else if ( C < 0.9 )
    {
        C = 0.9;
        celerity = C / (dT / deltaX);
    }
    Dn = (unitq / p->slope) / (celerity * deltaX);
    *c1 = (1 + C - Dn) / (1 + C + Dn);
    *c2 = (-1 + C + Dn) / (1 + C + Dn);
    *c3 = (1 - C + Dn) / (1 + C + Dn);
}

static void reportData(int i, double qUp, double* vel, double* dep)
//  Fv_CalPipeReportData he_dynamicrun.cpp:10419
{
    TIsPipe* p = &IsPipe[i];
    double pipeFlow, depth, v;
    if ( fabs(qUp) < IsFlowTol ) { *vel = 0.0;  *dep = 0.0;  return; }
    pipeFlow = qUp / p->count;
    if ( qUp < p->fullFlow ) depth = isewer_flowDepth(i, pipeFlow);
    else depth = p->dia;
    v = isewer_pipeVelocity(depth, 0.5 * p->dia, pipeFlow);
    if ( depth / p->dia < IS_DD_TOL ) { depth = 0.0;  v = 0.0; }   // HHU_FilterDepthVelocity
    *vel = v;  *dep = depth;
}

static int ensureWork(int n)
{
    if ( n <= WSize ) return 0;
    FREE(WFormer); FREE(WCurrent); FREE(WInflow); FREE(WOutflow); FREE(WIn0);
    FREE(WSimT); FREE(WHydT);
    WFormer = (double*) calloc(n, sizeof(double));
    WCurrent = (double*) calloc(n, sizeof(double));
    WInflow = (double*) calloc(n, sizeof(double));
    WOutflow = (double*) calloc(n, sizeof(double));
    WIn0 = (double*) calloc(n, sizeof(double));
    WSimT = (double*) calloc(n, sizeof(double));
    WHydT = (double*) calloc(n, sizeof(double));
    if ( !WFormer || !WCurrent || !WInflow || !WOutflow || !WIn0 || !WSimT || !WHydT )
    {
        report_writeErrorMsg(ERR_MEMORY, "");
        WSize = 0;
        return ErrorCode;
    }
    WSize = n;
    return 0;
}

static int routePipe(int k, int i, int downIdx, double** splitArr)
//  Fv_ManholePipeInit + Fv_TimeLoopAttenuate for one pipe, whole run
//  (eps_ref.py:route_pipe).
{
    TIsNode* nd = &IsNode[k];
    TIsPipe* p = &IsPipe[i];
    int segCount, K, seg, kk, m, hydroIndex, reportIndex, single, outletPipe,
        attenuate, haveMax;
    double deltaX, dt, simTime, hydTime, nextSimTime, lf0, lf1, c1, c2, c3,
           rmax, rmaxT;
    double *F, *C, *In, *Out, *In0, *DFi = DF[i];

    pipeInit(i);
    segCount = SegN[i];
    dt = SegDt[i];
    if ( p->length > IsMinSegment ) deltaX = p->length / segCount;
    else deltaX = p->length;

    // --- number of sub-steps (accumulated times, spec 2.3)
    K = 0;
    simTime = StartT;  hydTime = 0.0;
    while ( simTime <= (EndT + dt) ) { K++;  simTime = simTime + dt;  hydTime += dt; }
    if ( ensureWork(K + 4) ) return ErrorCode;
    F = WFormer;  C = WCurrent;  In = WInflow;  Out = WOutflow;  In0 = WIn0;
    simTime = StartT;  hydTime = 0.0;
    for ( kk = 1; kk <= K; kk++ )
    {
        WSimT[kk] = simTime;  WHydT[kk] = hydTime;
        simTime = simTime + dt;  hydTime += dt;
    }
    for ( kk = 0; kk < K + 4; kk++ ) F[kk] = C[kk] = In[kk] = Out[kk] = In0[kk] = 0.0;
    for ( kk = 0; kk <= NHyd; kk++ ) DFi[kk] = 0.0;
    DFi[0] = 0.0;
    hydroIndex = 1;
    reportIndex = 0;
    NRaw[i] = 0;
    haveMax = 0;  rmax = 0.0;  rmaxT = 0.0;
    single = (nd->ndown == 1);
    outletPipe = IsNode[p->to].isOutlet;
    attenuate = IsFlowAttenuation >= 1;   // "Flow Attenuation" >= 1e-4 (integer option)
    c1 = c2 = c3 = 0.0;

    for ( seg = 0; seg < segCount; seg++ )
    {
        lf0 = 0.0;
        for ( kk = 1; kk <= K; kk++ )
        {
            simTime = WSimT[kk];  hydTime = WHydT[kk];
            nextSimTime = simTime + dt;
            if ( seg == 0 )
            {
                if ( kk <= 1 ) lf0 = 0.0;
                lf1 = 0.0;
                if ( single )
                {
                    double nextHyd = nextSimTime - StartT, hT;
                    if ( kk <= 1 )
                    {
                        lf1 = nodeLoad(k, nextSimTime);
                        lf0 = nodeLoad(k, simTime);
                        lf0 += infilFlow(i, simTime);
                        lf1 += infilFlow(i, nextSimTime);
                        MbIn += lf0 * dt;
                    }
                    else
                    {
                        lf1 = nodeLoad(k, nextSimTime);
                        lf1 += infilFlow(i, nextSimTime);
                        MbIn += lf1 * dt;
                    }
                    // Fv_SinglePrecedeFlow: upstream pipes in input order
                    hT = simTime - StartT;
                    for ( m = 0; m < nd->nuph; m++ )
                    {
                        int up = nd->uph[m];
                        if ( kk <= 1 ) lf0 += interpHydro(DF[up], hT);
                        lf1 += interpHydro(DF[up], nextHyd);
                    }
                }
                else
                {
                    double* arr = splitArr[downIdx];
                    double hT = simTime - StartT;
                    if ( kk <= 1 ) lf0 = interpHydro(arr, hT);
                    lf1 = interpHydro(arr, nextSimTime - StartT);
                }
                if ( kk <= 1 )
                {
                    In[kk] = lf0;
                    F[kk] = In[kk];
                }
                In[kk + 1] = lf1;
                F[kk + 1] = In[kk + 1];
            }

            // --- segment outflow
            if ( attenuate && p->slope > 0.00001 )
            {
                if ( seg <= 0 && kk <= 10 ) mcCoeffs(i, F[kk], deltaX, &c1, &c2, &c3);
                if ( kk <= 1 ) C[kk] = 0.0;
                C[kk + 1] = c1 * F[kk] + c2 * F[kk + 1] + c3 * C[kk];
                if ( C[kk + 1] <= IsFlowTol ) C[kk + 1] = 0.0;
            }
            else if ( attenuate )
            {
                if ( kk <= 1 ) C[kk] = 0.0;
                C[kk + 1] = F[kk];
                if ( C[kk + 1] <= IsFlowTol ) C[kk + 1] = 0.0;
            }
            else                                    // Fv_SegmentOutflow
            {
                C[kk + 1] = F[kk];
                if ( C[kk + 1] <= IsFlowTol ) C[kk + 1] = 0.0;
            }
            In[kk] = C[kk];                         // QUIRK (5968)
            if ( seg == segCount - 1 )
            {
                Out[kk] = C[kk];
                Out[kk + 1] = C[kk + 1];
            }
            if ( seg <= 0 )
            {
                In[kk] = F[kk];
                In[kk + 1] = F[kk + 1];
            }
            if ( seg == segCount - 1 )
            {
                // Fv_DownstreamFlow (spec 4.4)
                if ( hydroIndex <= NHyd )
                {
                    if ( hydTime > (hydroIndex * H + DOUBLE_TIME_TOL) )
                    {
                        double dF = Out[kk] - DFi[hydroIndex - 1];
                        double dT = hydTime - (hydroIndex - 1) * H;
                        DFi[hydroIndex] = MAX(0.0, DFi[hydroIndex - 1] + H * dF / dT);
                        hydroIndex++;
                    }
                    else if ( fabs(hydTime - hydroIndex * H) < DOUBLE_TIME_TOL )
                    {
                        DFi[hydroIndex] = Out[kk];
                        hydroIndex++;
                    }
                }
                // Fv_ReportPipeData (spec 6.1)
                if ( outletPipe )
                {
                    if ( kk == 1 ) MbOut += Out[kk] * dt;
                    else           MbOut += Out[kk + 1] * dt;
                }
                if ( kk == 1 )
                {
                    TRptRaw* r = &Raw[i][0];
                    r->qUp = In[kk];  r->qDn = Out[kk];
                    r->qIn = (segCount == 1) ? F[kk] : In0[kk];
                    reportData(i, r->qUp, &r->vel, &r->dep);
                    reportIndex++;
                }
                else if ( hydTime >= reportIndex * RptStep && reportIndex <= NRpt )
                {
                    TRptRaw* r = &Raw[i][reportIndex];
                    double deltaTime = reportIndex * RptStep + dt - hydTime;
                    double* src = (segCount == 1) ? F : In0;
                    if ( fabs(deltaTime) <= ZERO_TIME_VALUE ) deltaTime = 0.0;
                    r->qUp = In[kk - 1] + deltaTime * (In[kk] - In[kk - 1]) / dt;
                    r->qDn = Out[kk - 1] + deltaTime * (Out[kk] - Out[kk - 1]) / dt;
                    r->qIn = src[kk - 1] + deltaTime * (src[kk] - src[kk - 1]) / dt;
                    reportData(i, r->qUp, &r->vel, &r->dep);
                    reportIndex++;
                }
                // Fv_OutputPipeFlowRange (spec 6.5)
                {
                    double q = C[kk];
                    if ( fabs(hydTime) < DOUBLE_TIME_TOL || !haveMax )
                    {
                        rmax = q;  rmaxT = simTime;  haveMax = 1;
                    }
                    else if ( q > rmax + IsFlowTol )
                    {
                        rmax = q;  rmaxT = simTime;
                    }
                }
            }
        }
        if ( seg == 0 ) for ( kk = 0; kk < K + 3; kk++ ) In0[kk] = F[kk];
        if ( seg < segCount - 1 )
        {
            // QUIRK: extra one-step shift between segments (5784)
            for ( kk = 1; kk <= K; kk++ ) F[kk] = C[kk - 1];
        }
    }
    NRaw[i] = reportIndex;
    MaxQ[i] = rmax;
    MaxT[i] = rmaxT;
    return 0;
}

//=============================================================================
//  Report-time outputs (spec 6)
//=============================================================================

static int computeOutputs(void)
//  Fv_DynamicData / Fv_ReadPipeMainResultFromFile / Fv_PipeBackwaterAdjustment
//  for every report step (eps_ref.py:outputs).
{
    int r, i, k, m, L = Nobjects[LINK], N = Nobjects[NODE];
    double *up, *dep, *vel, *upAdj, *head;
    int *ftype, *partial, *adjFull;

    up = (double*) calloc(L + 1, sizeof(double));
    dep = (double*) calloc(L + 1, sizeof(double));
    vel = (double*) calloc(L + 1, sizeof(double));
    upAdj = (double*) calloc(L + 1, sizeof(double));
    head = (double*) calloc(N + 1, sizeof(double));
    ftype = (int*) calloc(L + 1, sizeof(int));        // persists across reports
    partial = (int*) calloc(L + 1, sizeof(int));
    adjFull = (int*) calloc(L + 1, sizeof(int));
    if ( !up || !dep || !vel || !upAdj || !head || !ftype || !partial || !adjFull )
    {
        FREE(up); FREE(dep); FREE(vel); FREE(upAdj); FREE(head); FREE(ftype);
        FREE(partial); FREE(adjFull);
        report_writeErrorMsg(ERR_MEMORY, "");
        return ErrorCode;
    }
    PipeStorage = 0.0;
    for ( r = 0; r < NRpt; r++ )
    {
        for ( i = 0; i < L; i++ )
        {
            TRptRaw* w = &Raw[i][MIN(r, NRaw[i] > 0 ? NRaw[i] - 1 : 0)];
            TIsPipe* p = &IsPipe[i];
            double radius = p->dia / 2.0, wd, a, fr;
            up[i] = w->qUp;  dep[i] = w->dep;  vel[i] = w->vel;
            OUp[r][i] = w->qUp;  ODn[r][i] = w->qDn;  OIn[r][i] = w->qIn;
            wd = isewer_wetWidth(dep[i], radius);
            a = isewer_wetArea(dep[i], radius);
            if ( wd <= 0.0 && dep[i] > 0.0 ) fr = vel[i] / sqrt(32.2 * dep[i]);
            else if ( a > 0.0 && wd > 0.0 ) fr = vel[i] / sqrt(32.2 * (a / wd));
            else fr = 0.0;
            OFroude[r][i] = fr;
            if ( r == NRpt - 1 && r > 0 ) PipeStorage += a * p->count * p->length;
        }
        isewer_runHead(up, dep, head, ftype, partial, adjFull, 0);

        // --- Fv_PipeFlow: QUIRK adverse pipes take the sum of the pipes
        //     entering their from-node, in place, in input order
        for ( i = 0; i < L; i++ ) upAdj[i] = up[i];
        for ( i = 0; i < L; i++ )
        {
            TIsPipe* p = &IsPipe[i];
            if ( p->adverse )
            {
                double f = 0.0;
                TIsNode* fn = &IsNode[p->frm];
                for ( m = 0; m < fn->nuph; m++ ) f += upAdj[fn->uph[m]];
                upAdj[i] = f;
                OFlow[r][i] = MAX(0.0, f);
            }
            else OFlow[r][i] = MAX(0.0, up[i]);
        }
        for ( i = 0; i < L; i++ )
        {
            TIsPipe* p = &IsPipe[i];
            double d, ff, v;
            if ( p->adverse )
            {
                d = MIN(p->dia, dep[i]);
                OVel[r][i] = isewer_pipeVelocity(d, 0.5 * p->dia, upAdj[i]);
            }
            else OVel[r][i] = vel[i];
            ODepth[r][i] = MIN(p->dia, dep[i]);
            OdD[r][i] = ODepth[r][i] / p->dia;
            ff = isewer_linkFlow(i, p->dia, -1.0, -1.0);
            v = MAX(0.0, OFlow[r][i] / ff);
            if ( OdD[r][i] >= 1.0 ) v = MAX(v, 1.0);
            OqQ[r][i] = v;
        }
        // --- Fv_PipeBackwaterAdjustment (12760)
        //     QUIRK: adjusted velocity always uses the full-pipe area
        for ( i = 0; i < L; i++ )
        {
            TIsPipe* p = &IsPipe[i];
            double uh = head[p->frm], dh = head[p->to], ad;
            if ( dh > dep[i] + p->toInv )
            {
                ad = 0.5 * ((dh - p->toInv) + (uh - p->fromInv));
                ad = MIN(p->dia, ad);
                ad = MAX(ad, dep[i]);
                OAdjD[r][i] = ad;
                OAdjV[r][i] = isewer_pipeVelocity(p->dia, 0.5 * p->dia, MAX(0.0, upAdj[i]));
                OBack[r][i] = 1.0;
            }
            else
            {
                OAdjD[r][i] = dep[i];
                OAdjV[r][i] = OVel[r][i];
                OBack[r][i] = 0.0;
            }
            OFtype[r][i] = ftype[i] == 1 ? 1.0 : -1.0;
            OPartial[r][i] = partial[i] ? 1.0 : -1.0;
        }
        // --- nodes and outlets
        for ( k = 0; k < N; k++ )
        {
            TIsNode* nd = &IsNode[k];
            OHead[r][k] = head[k];
            if ( nd->isOutlet )
            {
                double q = 0.0;
                for ( m = 0; m < nd->nuph; m++ ) q += ODn[r][nd->uph[m]];
                OOutflow[r][k] = q;
                continue;
            }
            OStatus[r][k] = head[k] >= nd->rim ? 1.0 : 0.0;
            {
                int upSuper = 0, dnSub = 0;
                for ( m = 0; m < nd->nuph; m++ ) if ( OFroude[r][nd->uph[m]] > 1 ) upSuper = 1;
                for ( m = 0; m < nd->ndown; m++ ) if ( OFroude[r][nd->down[m]] < 1 ) dnSub = 1;
                OJump[r][k] = (upSuper && dnSub) ? 1.0 : 0.0;
            }
            OLoad[r][k] = nodeLoadTotal(k, StartT + r * RptStep);
        }
    }
    FREE(up); FREE(dep); FREE(vel); FREE(upAdj); FREE(head); FREE(ftype);
    FREE(partial); FREE(adjFull);
    return 0;
}

//=============================================================================
//  SWMM state
//=============================================================================

void isewer_epsSetState(double tSec)
//  State at tSec (seconds from the start): the report-time result when tSec
//  is a report time, otherwise linear interpolation between the two around it.
//  InfoSewer routes by translation (or Muskingum-Cunge), so the water held in
//  a conduit is what has entered and not yet left: its SWMM volume is the
//  integral of inlet minus outlet flow (trapezoidal, like SWMM's mass
//  balance), starting empty as InfoSewer does. Depth stays InfoSewer's
//  reported normal depth.
{
    int r0, r1, i, j, L = Nobjects[LINK], N = Nobjects[NODE];
    double x, f;
    if ( !OFlow ) return;
    x = tSec / RptStep;
    r0 = (int)floor(x + 1.0e-9);
    f = x - r0;
    if ( fabs(f) < 1.0e-9 ) f = 0.0;
    if ( r0 >= NRpt - 1 ) { r0 = NRpt - 1;  f = 0.0; }
    if ( r0 < 0 ) { r0 = 0;  f = 0.0; }
    r1 = (f > 0.0) ? r0 + 1 : r0;
    for ( i = 0; i < L; i++ )
    {
        int k = Link[i].subIndex;
        double q, y, a, qin, qout;
        if ( f == 0.0 )
        {
            q = OFlow[r0][i];  y = ODepth[r0][i];
            qin = OIn[r0][i];  qout = ODn[r0][i];
        }
        else
        {
            q = (1.0 - f) * OFlow[r0][i] + f * OFlow[r1][i];
            y = (1.0 - f) * ODepth[r0][i] + f * ODepth[r1][i];
            qin = (1.0 - f) * OIn[r0][i] + f * OIn[r1][i];
            qout = (1.0 - f) * ODn[r0][i] + f * ODn[r1][i];
        }
        qin = MAX(qin, 0.0);
        qout = MAX(qout, 0.0);
        if ( CurT < 0.0 ) CurVol[i] = 0.0;
        else CurVol[i] += 0.5 * ((CurQin[i] + qin) - (CurQout[i] + qout)) * (tSec - CurT);
        CurQin[i] = qin;
        CurQout[i] = qout;
        a = (y > 0.0) ? xsect_getAofY(&Link[i].xsect, y) : 0.0;
        Link[i].newFlow = q;
        Link[i].newDepth = y;
        Conduit[k].q1 = CurQin[i] / Conduit[k].barrels;
        Conduit[k].q2 = CurQout[i] / Conduit[k].barrels;
        Conduit[k].a1 = Conduit[k].a2 = a;
        Link[i].newVolume = CurVol[i];
    }
    CurT = tSec;
    for ( j = 0; j < N; j++ )
    {
        double h = (f == 0.0) ? OHead[r0][j] : (1.0 - f) * OHead[r0][j] + f * OHead[r1][j];
        double y = h - Node[j].invertElev;
        if ( y < 0.0 ) y = 0.0;
        Node[j].newDepth = y;
        Node[j].newVolume = 0.0;
        Node[j].overflow = 0.0;
    }
}

void isewer_epsNodeFlows(void)
//  Node inflow/outflow totals for SWMM's mass balance: a conduit takes its
//  true inlet flow from its upstream node and delivers its outlet flow to
//  its downstream node.
{
    int i;
    for ( i = 0; i < Nobjects[LINK]; i++ )
    {
        Node[Link[i].node1].outflow += CurQin[i];
        Node[Link[i].node2].inflow  += CurQout[i];
    }
}

void isewer_epsLoads(double tSec)
//  InfoSewer loads and pipe infiltration at tSec as SWMM lateral inflow
//  (dry weather inflow), each pipe's infiltration at its upstream node.
{
    int k, m;
    double t = StartT + tSec, q;
    if ( !OFlow ) return;
    for ( k = 0; k < Nobjects[NODE]; k++ )
    {
        TIsNode* nd = &IsNode[k];
        if ( nd->isOutlet ) continue;
        q = nodeLoadTotal(k, t);
        for ( m = 0; m < nd->ndown; m++ ) q += infilFlow(nd->down[m], t);
        if ( q == 0.0 ) continue;
        Node[k].newLatFlow += q;
        massbal_addInflowFlow(DRY_WEATHER_INFLOW, q);
    }
}

//=============================================================================
//  Report and results file
//=============================================================================

static void milTime(double t, char* s)
//  HU_GetMilTime: t = hours + 0.0002777, hours and minutes truncated
{
    double h = t / 3600.0 + 0.0002777;
    int hh = (int)h, mm = (int)((h - hh) * 60.0);
    sprintf(s, "%3d:%02d", hh, mm);
}

void isewer_epsReport(void)
{
    int i, L = Nobjects[LINK];
    double ff = IsFlowFactor, lf = UCF(LENGTH), check, stor;
    FILE* f = Frpt.file;
    char ts[32];
    if ( !f || RptFlags.disabled || !OFlow ) return;

    stor = PipeStorage;
    check = MbIn != 0.0 ? fabs(100.0 * (MbIn - (MbOut + stor)) / MbIn) : 0.0;
    fprintf(f,
"\n  ***************************************"
"\n  InfoSewer EPS Maximum Gravity Pipe Flows"
"\n  ***************************************"
"\n  Maximum over every routing sub-step; depth is the normal depth of that flow.\n");
    fprintf(f,
"\n  -----------------------------------------------------------------------------------"
"\n                           Max Flow     Time   Velocity       d/D   Segments   Sub-step"
"\n  Pipe                     %8s    hh:mm     %6s                            sec"
"\n  -----------------------------------------------------------------------------------",
    FlowUnitWords[FlowUnits], UnitSystem == US ? "ft/s" : "m/s");
    for ( i = 0; i < L; i++ )
    {
        TIsPipe* p = &IsPipe[i];
        double q = MaxQ[i] / p->count, d, v;
        d = isewer_flowDepth(i, q);
        v = isewer_pipeVelocity(d, 0.5 * p->dia, q);
        if ( d / p->dia < IS_DD_TOL ) { d = 0.0;  v = 0.0; }
        milTime(MaxT[i], ts);
        fprintf(f, "\n  %-20s %12.4f   %s %10.3f %9.3f %10d %10.1f", Link[i].ID,
            MaxQ[i] * ff, ts, v * lf, d / p->dia, SegN[i], SegDt[i]);
    }
    fprintf(f,
"\n\n  ***************************************"
"\n  InfoSewer EPS Mass Balance (as InfoSewer)"
"\n  ***************************************"
"\n  Inflow volume ............ %14.3f ft3"
"\n  Outflow volume ........... %14.3f ft3"
"\n  Final pipe storage ....... %14.3f ft3"
"\n  Check (%%) ................ %14.3f"
"\n  InfoSewer integrates inflow to the end time plus one sub-step and counts"
"\n  split manholes for one extra hydraulic step; see the README.\n",
    MbIn, MbOut, stor, check);
}

static void writeCsv(void)
//  Every InfoSewer EPS output array at full precision, internal units
//  (cfs, ft), one block per report time, for verification.
{
    int r, i, k, L = Nobjects[LINK], N = Nobjects[NODE];
    FILE* f = fopen(IsResultsFile, "wt");
    if ( !f ) { report_writeErrorMsg(ERR_RPT_FILE, IsResultsFile); return; }
    fprintf(f, "kind,report,id,flow,vel,depth,froude,dD,overflow,backwater,adjdepth,"
               "adjvel,qQ,flowtype,partial,up,dn,in\n");
    for ( r = 0; r < NRpt; r++ )
        for ( i = 0; i < L; i++ )
            fprintf(f, "pipe,%d,%s,%.17g,%.17g,%.17g,%.17g,%.17g,0,%.17g,%.17g,%.17g,"
                       "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                r, Link[i].ID, OFlow[r][i], OVel[r][i], ODepth[r][i], OFroude[r][i],
                OdD[r][i], OBack[r][i], OAdjD[r][i], OAdjV[r][i], OqQ[r][i],
                OFtype[r][i], OPartial[r][i], OUp[r][i], ODn[r][i], OIn[r][i]);
    fprintf(f, "kind,report,id,load,head,status,jump\n");
    for ( r = 0; r < NRpt; r++ )
        for ( k = 0; k < N; k++ )
            if ( !IsNode[k].isOutlet )
                fprintf(f, "node,%d,%s,%.17g,%.17g,%.17g,%.17g\n", r, Node[k].ID,
                    OLoad[r][k], OHead[r][k], OStatus[r][k], OJump[r][k]);
    fprintf(f, "kind,report,id,outflow,outhead\n");
    for ( r = 0; r < NRpt; r++ )
        for ( k = 0; k < N; k++ )
            if ( IsNode[k].isOutlet )
                fprintf(f, "outlet,%d,%s,%.17g,%.17g\n", r, Node[k].ID,
                    OOutflow[r][k], OHead[r][k]);
    fprintf(f, "kind,id,maxflow,maxtime,segments,dt\n");
    for ( i = 0; i < L; i++ )
        fprintf(f, "max,%s,%.17g,%.17g,%d,%.17g\n", Link[i].ID, MaxQ[i], MaxT[i],
            SegN[i], SegDt[i]);
    fprintf(f, "kind,in,out,storage\nmassbal,%.17g,%.17g,%.17g\n", MbIn, MbOut, PipeStorage);
    fclose(f);
}

//=============================================================================

void isewer_epsClose(void)
{
    int p, L = Nobjects[LINK];
    if ( PatV ) for ( p = 0; p < IsNpats; p++ ) FREE(PatV[p]);
    FREE(PatV); FREE(PatN);
    free2(DF, L);  DF = NULL;
    if ( Raw ) for ( p = 0; p < L; p++ ) FREE(Raw[p]);
    FREE(Raw); FREE(NRaw); FREE(MaxQ); FREE(MaxT); FREE(SegN); FREE(SegDt);
    FREE(CurQin); FREE(CurQout); FREE(CurVol);
    free2(OFlow, NRpt); free2(OVel, NRpt); free2(ODepth, NRpt); free2(OdD, NRpt);
    free2(OqQ, NRpt); free2(OFroude, NRpt); free2(OBack, NRpt); free2(OAdjD, NRpt);
    free2(OAdjV, NRpt); free2(OFtype, NRpt); free2(OPartial, NRpt); free2(OUp, NRpt);
    free2(ODn, NRpt); free2(OIn, NRpt); free2(OHead, NRpt); free2(OLoad, NRpt);
    free2(OStatus, NRpt); free2(OJump, NRpt); free2(OOutflow, NRpt);
    OFlow = OVel = ODepth = OdD = OqQ = OFroude = OBack = OAdjD = OAdjV = OFtype =
    OPartial = OUp = ODn = OIn = OHead = OLoad = OStatus = OJump = OOutflow = NULL;
    FREE(WFormer); FREE(WCurrent); FREE(WInflow); FREE(WOutflow); FREE(WIn0);
    FREE(WSimT); FREE(WHydT);
    WSize = 0;
}
