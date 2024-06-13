#include "stdDefinitions.h"
#include "MagFieldOps.h"
#include "mfoGlobals.h"
#include "agpWiegelmannPar.h"
#include "agmRotate3D.h"
#include "TimeTicToc.h"

#include "console_debug.h"
#include "debug_data_trace_win.h"
#include "DebugWrite.h"

static aguTimeTicToc tic;

static double get_max_L_incr(int *N)
{
    return DBL_EPSILON * sqrt(N[0] * N[1] * N[2]);
}

static int proceedDL(double /* dL */, double L, int stepN, double * /* mL */, double *mAv)
{
    if (stepN == 0)
        mAv[stepN] = 0.0;
    else if (stepN < WiegelmannProcdLStdWin)
        mAv[stepN] = L;
    else
    {
        for (int k = 1; k < WiegelmannProcdLStdWin; k++)
            mAv[k-1] = mAv[k];
        mAv[WiegelmannProcdLStdWin-1] = L;

        double s = 0, s2 = 0;
        for (int k = 0; k < WiegelmannProcdLStdWin; k++)
        {
            s += mAv[k];
            s2 += mAv[k]*mAv[k];
        }

        double w = WiegelmannProcdLStdWin;
        double var = sqrt(w/(w-1)*(w*s2/(s*s)-1));
        if (var < WiegelmannProcdLStdVal)
            return 7;
    }

     return 0;
}

uint32_t mfoWiegelmannProcedureCore(CagmVectorField *field, CagmScalarField *weight, 
    CagmVectorField *baseField, CagmVectorField *baseWeight, CagmVectorField *baseField2, CagmVectorField *baseWeight2,
    CagmScalarField *absField, CagmScalarField *absWeight, CagmScalarField *absField2, CagmScalarField *absWeight2,
    CagmScalarField *losField, CagmScalarField *losWeight, CagmScalarField *losField2, CagmScalarField *losWeight2,
    CagmScalarField *bottomWeight,
    double *vcos, int depth, PROTO_mfoWiegelmannCallback callback, double *maxStep)
{
    int N[3];
    field->dimensions(N);

    double *memoryAv = new double[WiegelmannProcdLStdWin];

    CagmVectorField *vF = new CagmVectorField(field);

    int NB[3];
    memcpy(NB, N, 3*sizeof(int));
    NB[0] = 2;
    CagmVectorField boundsx(NB);
    memcpy(NB, N, 3*sizeof(int));
    NB[1] = 2;
    CagmVectorField boundsy(NB);
    memcpy(NB, N, 3*sizeof(int));
    NB[2] = 2;
    CagmVectorField boundsz(NB);
    field->getBounds(&boundsx, &boundsy, &boundsz);
 
    CagmRotate3D rotator(vcos);

    if (WiegelmannProcCondType == 2)
        field->condWeight(WiegelmannProcCondBase, baseField, baseWeight, WiegelmannProcCondBase2, baseField2, baseWeight2,
            WiegelmannProcCondAbs, absField, absWeight, WiegelmannProcCondAbs2, absField2, absWeight2,
            WiegelmannProcCondLOS, losField, losWeight, WiegelmannProcCondLOS2, losField2, losWeight2, &rotator);

    CagpWiegelmann *proc = new CagpWiegelmann(N, CommonThreadsN, WiegelmannDerivStencil
        , field, weight
        , vcos
        , vF
        , baseField, baseWeight
        , absField, absWeight
        , losField, losWeight
        , bottomWeight
        , depth
        , WiegelmannThreadPriority
        );

    DebugWriteData(field, "debug_B", depth);

    proc->max_dL_incr = get_max_L_incr(N);

    int stepN = 0;
    int stop = 0;

    int iterN = 0;
    double L0 = proc->step(iterN);
    iterN++;
    double step0 = sqrt(proc->B2sum[0]/proc->F2max[0])/(N[0]*N[1]*N[2]) * WiegelmannProcStep0;

    DebugWriteData(vF, "debug_vF", depth);

    double L, Lprev = L0;
    proceedDL(0.0, 0.0, 0, nullptr, memoryAv);

    CagmVectorField *prevField = new CagmVectorField(field, true);
    CagmVectorField *prevVF = new CagmVectorField(vF, true);

    double step = step0, dL;
    *maxStep = step0;

    int reason = 0;

    double Lptest = L0;
    while (true)
    {
        field->Copy(prevField);
        vF->Copy(prevVF);
        ///////////////////////
        vF->mult(step);
        field->add(vF);
        field->setBounds(&boundsx, &boundsy, &boundsz);
        ///////////////////////
        L = proc->step(iterN);
        iterN++;

        if (L >= Lptest)
            break;

        Lptest = L;
        step *= 2;
    }
    step0 = step * WiegelmannProcStep0;
    step = step0;
    field->Copy(prevField);
    vF->Copy(prevVF);

    int z_size = N[2];
    int voxels = N[0]*N[1]*N[2];

    proc->step_failed[0] = 0;

    while (true)
    {
        ///////////////////////
        vF->mult(step);
        field->add(vF);
        field->setBounds(&boundsx, &boundsy, &boundsz);
        ///////////////////////
        if (WiegelmannProcCondType == 2)
            field->condWeight(WiegelmannProcCondBase, baseField, baseWeight, WiegelmannProcCondBase2, baseField2, baseWeight2,
                WiegelmannProcCondAbs, absField, absWeight, WiegelmannProcCondAbs2, absField2, absWeight2,
                WiegelmannProcCondLOS, losField, losWeight, WiegelmannProcCondLOS2, losField2, losWeight2, &rotator);

        L = proc->step(iterN);
        iterN++;

        dL = Lprev - L;
        if (dL/Lprev <= - proc->max_dL_incr)
        {
            field->Copy(prevField);
            vF->Copy(prevVF);
            if (step < step0*WiegelmannProcStepLim)
                reason = 1;
            else
                step *= WiegelmannProcStepDecr;
            (proc->step_failed[0])++;
        }
        else
        {
            stepN++;

            proc->setBase(L/L0, step, stepN, depth, iterN, tic.toc(), z_size, voxels);

            Lprev = L;
            prevField->Copy(field);
            prevVF->Copy(vF);

            step *= WiegelmannProcStepIncr;
            if (step > step0*WiegelmannProcStepMax)
                step = step0*WiegelmannProcStepMax;
            if (step > *maxStep)
                *maxStep = step;

            if (stepN > WiegelmannProcMaxSteps)
                reason = 3;
            else
                reason = proceedDL(dL/L0, L/L0, stepN, nullptr, memoryAv);

            if (callback && stepN%WiegelmannProtocolStep == 0)
                callback(step/step0, 0, 0, 0, depth, dL/L0/step, proc, field, &stop);

            proc->step_failed[0] = 0;
            // theta etc. reasons from 8
        }

        if (reason != 0)
        {
            field->Copy(prevField);
            vF->Copy(prevVF);
            break;
        }

        if (stop)
        {
            reason = 4;
            break;
        }
    }

    field->setBounds(&boundsx, &boundsy, &boundsz);
    if (WiegelmannProcCondType == 2)
        field->condWeight(WiegelmannProcCondBase, baseField, baseWeight, WiegelmannProcCondBase2, baseField2, baseWeight2,
            WiegelmannProcCondAbs, absField, absWeight, WiegelmannProcCondAbs2, absField2, absWeight2,
            WiegelmannProcCondLOS, losField, losWeight, WiegelmannProcCondLOS2, losField2, losWeight2, &rotator);

    stop = reason;
    proc->setBase(L/L0, step, stepN, depth, iterN, tic.toc(), z_size, voxels);
    if (callback)
        callback(step/step0, 0, 0, 0, depth, dL/L0, proc, field, &stop);

    delete prevField;
    delete prevVF;
    delete vF;
    delete proc;
    delete [] memoryAv;

    return reason;
}

static int xfloor(double v)
{
    int c = (int)ceil(v);
    if (c - v < 0.5)
        return c;
    else
        return c-1;
}

__declspec( dllexport ) uint32_t mfoWiegelmannProcedure(CagmVectorField *field, CagmScalarField *weight, 
    CagmVectorField *baseField, CagmVectorField *baseWeight, CagmVectorField *baseField2, CagmVectorField *baseWeight2,
    CagmScalarField *absField, CagmScalarField *absWeight, CagmScalarField *absField2, CagmScalarField *absWeight2, 
    CagmScalarField *losField, CagmScalarField *losWeight, CagmScalarField *losField2, CagmScalarField *losWeight2,
    CagmScalarField *bottomWeight,
    double *vcos, PROTO_mfoWiegelmannCallback callback)
{
    DebugWritePars("debug_input", field, weight);
    DebugWriteData(field, "debug_B");
    DebugWriteData(weight, "debug_W");

    console_start();

    if (!baseField || !baseWeight)
        WiegelmannProcCondBase = 0;
    if (!baseField2 || !baseWeight2)
        WiegelmannProcCondBase2 = 0;
    if (!absField || !absWeight)
        WiegelmannProcCondAbs = 0;
    if (!absField2 || !absWeight2)
        WiegelmannProcCondAbs2 = 0;
    if (!losField || !losWeight)
        WiegelmannProcCondLOS = 0;
    if (!losField2 || !losWeight2)
        WiegelmannProcCondLOS2 = 0;
    if (WiegelmannProcCondBase == 0 && WiegelmannProcCondBase2 == 0 &&
        WiegelmannProcCondAbs == 0 && WiegelmannProcCondAbs2 == 0 &&
        WiegelmannProcCondLOS == 0 && WiegelmannProcCondLOS2 == 0)
        WiegelmannProcCondType = 0;

    double maxStep;
    bool bMatr = WiegelmannMatryoshkaUse != 0;
    int depth = 1;
    if (bMatr)
    {
        depth = field->getMatryoshkaDepth(WiegelmannMatryoshkaDeepMinN, WiegelmannMatryoshkaFactor);
        if (depth < 2)
            bMatr = false;
    }

    tic.tic();

    uint32_t dwRC = 0;
    if (!bMatr)
    {
        WiegelmannProcStepIncr = WiegelmannProcStepIncrMain;
        WiegelmannProcStepDecr = WiegelmannProcStepDecrMain;
        WiegelmannProcStepLim = WiegelmannProcStepLimMain;
        WiegelmannProcdLStdVal = WiegelmannProcdLStdValMain;
        WiegelmannProcdLStdWin = WiegelmannProcdLStdWinMain;
        dwRC = mfoWiegelmannProcedureCore(field, weight, baseField, baseWeight, baseField2, baseWeight2,
            absField, absWeight, absField2, absWeight2,
            losField, losWeight, losField2, losWeight2, bottomWeight, vcos, 1, callback, &maxStep);
    }
    else
    {
        double factor = WiegelmannMatryoshkaFactor;
        int *Ns = new int[3*depth];
        double *steps = new double[3*depth];
        field->dimensions(Ns);
        steps[0] = 1.0;
        steps[1] = 1.0;
        steps[2] = 1.0;
        for (int i = 1; i < depth; i++)
        {
            Ns[3*i  ] = xfloor( (Ns[3*(i-1)  ]-1.0)/factor + 1.0 ); 
            Ns[3*i+1] = xfloor( (Ns[3*(i-1)+1]-1.0)/factor + 1.0 ); 
            Ns[3*i+2] = xfloor( (Ns[3*(i-1)+2]-1.0)/factor + 1.0 ); 
        }
        double cf = factor;
        for (int i = 1; i < depth; i++)
        {
            steps[3*i  ] = (Ns[3*i  ]-1.0)/(Ns[0]-1.0);
            steps[3*i+1] = (Ns[3*i+1]-1.0)/(Ns[1]-1.0);
            steps[3*i+2] = (Ns[3*i+2]-1.0)/(Ns[2]-1.0);
            cf *= factor;
        }

        WiegelmannProcStepIncr = WiegelmannProcStepIncrInit;
        WiegelmannProcStepDecr = WiegelmannProcStepDecrInit;
        WiegelmannProcStepLim = WiegelmannProcStepLimInit;
        WiegelmannProcdLStdVal = WiegelmannProcdLStdValInit;
        WiegelmannProcdLStdWin = WiegelmannProcdLStdWinInit;

        int *pN0 = Ns+3*(depth-1);
        CagmVectorField *v0 = new CagmVectorField(pN0, steps+3*(depth-1));

        DebugWriteData(field, "field_init", depth);
        v0->stretch(field);
        DebugWriteData(v0, "field_stretch", depth);

        DebugWriteData(weight, "weight_init", depth);
        CagmScalarField *sW0 = new CagmScalarField(pN0, steps+3*(depth-1));
        sW0->stretch(weight);
        DebugWriteData(sW0, "weight_stretch", depth);

        mfoWiegelmannProcedureCore(v0, sW0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, vcos, depth, callback, &maxStep);
        delete sW0;
        CagmVectorField *v1;
        for (int i = depth-2; i > 0; i--)
        {
            WiegelmannProcStepIncr = WiegelmannProcStepIncrMatr;
            WiegelmannProcStepDecr = WiegelmannProcStepDecrMatr;
            WiegelmannProcStepLim = WiegelmannProcStepLimMatr;
            WiegelmannProcdLStdVal = WiegelmannProcdLStdValMatr;
            WiegelmannProcdLStdWin = WiegelmannProcdLStdWinMatr;
            pN0 = Ns+3*i;
            v1 = new CagmVectorField(pN0, steps+3*i);
            v1->stretch(v0);
            DebugWriteData(v1, "field_stretch", i+1);

            delete v0;
            sW0 = new CagmScalarField(pN0, steps+3*i);
            sW0->stretch(weight);
            DebugWriteData(sW0, "weight_stretch", i+1);
            mfoWiegelmannProcedureCore(v1, sW0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, vcos, i+1, callback, &maxStep);
            delete sW0;
            v0 = v1;
        }

        int NB[3];
        field->dimensions(NB);
        NB[0] = 2;
        CagmVectorField boundsx(NB);
        field->dimensions(NB);
        NB[1] = 2;
        CagmVectorField boundsy(NB);
        field->dimensions(NB);
        NB[2] = 2;
        CagmVectorField boundsz(NB);
        field->getBounds(&boundsx, &boundsy, &boundsz);
        field->stretch(v0);
        field->setBounds(&boundsx, &boundsy, &boundsz);
        DebugWriteData(field, "field_stretch", 1);

        delete v0;

        delete [] Ns;
        delete [] steps;

        WiegelmannProcStepIncr = WiegelmannProcStepIncrMain;
        WiegelmannProcStepDecr = WiegelmannProcStepDecrMain;
        WiegelmannProcStepLim = WiegelmannProcStepLimMain;

        WiegelmannProcdLStdVal = WiegelmannProcdLStdValMain;
        WiegelmannProcdLStdWin = WiegelmannProcdLStdWinMain;

        dwRC = mfoWiegelmannProcedureCore(field, weight, baseField, baseWeight, baseField2, baseWeight2, 
                                          absField, absWeight, absField2, absWeight2, 
                                          losField, losWeight, losField2, losWeight2, bottomWeight, vcos, 1, callback, &maxStep);
    }

    return dwRC;
}
