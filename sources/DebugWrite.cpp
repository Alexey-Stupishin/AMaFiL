#include "stdDefinitions.h"
#include "mfoGlobals.h"
#include "agmScalarField.h"
#include "agmVectorField.h"
#include "binUtilities.h"

#ifdef _WINDOWS
#pragma warning(disable:4996)
#endif

#include "DebugWrite.h"
#include "debug_data_trace_win.h"

void DebugWriteData(CubeXD *v, const char *fname, int depth, int iter)
{
#ifdef _WINDOWS
    if (debug_input)
    {
        char buffer[256], tbuff[32];
        strcpy(buffer, DEBUG_OUT_PATH);
        strcat(buffer, fname);
        if (depth > 0)
        {
            sprintf(tbuff, "_%d", depth);
            strcat(buffer, tbuff);
            if (iter > 0)
            {
                sprintf(tbuff, "_%d", iter);
                strcat(buffer, tbuff);
            }
        }
        strcat(buffer, ".bin");

        debug_data_trace_win *d = new debug_data_trace_win(buffer);
        d->write(v);
        delete d;
    }
#endif
}

void DebugWritePars(const char *fname, CagmVectorField * field, CagmScalarField * w)
{
#ifdef _WINDOWS
    if (debug_input)
    {
        char buffer[256];
        strcpy(buffer, DEBUG_OUT_PATH);
        strcat(buffer, fname);
        strcat(buffer, ".bin");

        FILE *fid = fopen(buffer, "wb");

        CbinDataStruct::WriteHeader(fid);

        CbinDataStruct::Write(fid, &CommonThreadsN, 1, "CommonThreadsN");
        CbinDataStruct::Write(fid, &WiegelmannThreadPriority, 1, "WiegelmannThreadPriority");

        CbinDataStruct::Write(fid, &WiegelmannBoundsCorrection, 1, "WiegelmannBoundsCorrection");

        CbinDataStruct::Write(fid, &WiegelmannWeightType, 1, "WiegelmannWeightType");
        CbinDataStruct::Write(fid, &WiegelmannWeightBound, 1, "WiegelmannWeightBound");
        CbinDataStruct::Write(fid, &WiegelmannWeightDivfree, 1, "WiegelmannWeightDivfree");

        CbinDataStruct::Write(fid, &WiegelmannDerivStencil, 1, "WiegelmannDerivStencil");
        CbinDataStruct::Write(fid, &WiegelmannInversionTolerance, 1, "WiegelmannInversionTolerance");
        CbinDataStruct::Write(fid, &WiegelmannInversionDenom, 1, "WiegelmannInversionDenom");
        
        CbinDataStruct::Write(fid, &WiegelmannProcStep0, 1, "WiegelmannProcStep0");

        CbinDataStruct::Write(fid, &WiegelmannProcStepMax, 1, "WiegelmannProcStepMax");
        CbinDataStruct::Write(fid, &WiegelmannProcMaxSteps, 1, "WiegelmannProcMaxSteps");

        CbinDataStruct::Write(fid, &WiegelmannProcStepIncrInit, 1, "WiegelmannProcStepIncrInit");
        CbinDataStruct::Write(fid, &WiegelmannProcStepIncrMatr, 1, "WiegelmannProcStepIncrMatr");
        CbinDataStruct::Write(fid, &WiegelmannProcStepIncrMain, 1, "WiegelmannProcStepIncrMain");

        CbinDataStruct::Write(fid, &WiegelmannProcStepDecrInit, 1, "WiegelmannProcStepDecrInit");
        CbinDataStruct::Write(fid, &WiegelmannProcStepDecrMatr, 1, "WiegelmannProcStepDecrMatr");
        CbinDataStruct::Write(fid, &WiegelmannProcStepDecrMain, 1, "WiegelmannProcStepDecrMain");

        CbinDataStruct::Write(fid, &WiegelmannProcStepLimInit, 1, "WiegelmannProcStepLimInit");
        CbinDataStruct::Write(fid, &WiegelmannProcStepLimMatr, 1, "WiegelmannProcStepLimMatr");
        CbinDataStruct::Write(fid, &WiegelmannProcStepLimMain, 1, "WiegelmannProcStepLimMain");

        CbinDataStruct::Write(fid, &WiegelmannProcdLStdWinInit, 1, "WiegelmannProcdLStdWinInit");
        CbinDataStruct::Write(fid, &WiegelmannProcdLStdWinMatr, 1, "WiegelmannProcdLStdWinMatr");
        CbinDataStruct::Write(fid, &WiegelmannProcdLStdWinMain, 1, "WiegelmannProcdLStdWinMain");

        CbinDataStruct::Write(fid, &WiegelmannProcdLStdValInit, 1, "WiegelmannProcdLStdValInit");
        CbinDataStruct::Write(fid, &WiegelmannProcdLStdValMatr, 1, "WiegelmannProcdLStdValMatr");
        CbinDataStruct::Write(fid, &WiegelmannProcdLStdValMain, 1, "WiegelmannProcdLStdValMain");

        CbinDataStruct::Write(fid, &WiegelmannMatryoshkaUse, 1, "WiegelmannMatryoshkaUse");
        CbinDataStruct::Write(fid, &WiegelmannMatryoshkaDeepMinN, 1, "WiegelmannMatryoshkaDeepMinN");
        CbinDataStruct::Write(fid, &WiegelmannMatryoshkaFactor, 1, "WiegelmannMatryoshkaFactor");

        CbinDataStruct::Write(fid, &WiegelmannProcCondType, 1, "WiegelmannProcCondType");
        CbinDataStruct::Write(fid, &WiegelmannProcCondAbs, 1, "WiegelmannProcCondAbs");
        CbinDataStruct::Write(fid, &WiegelmannProcCondAbs2, 1, "WiegelmannProcCondAbs2");
        CbinDataStruct::Write(fid, &WiegelmannProcCondLOS2, 1, "WiegelmannProcCondLOS2");
        CbinDataStruct::Write(fid, &WiegelmannProcCondLOS, 1, "WiegelmannProcCondLOS");
        CbinDataStruct::Write(fid, &WiegelmannProcCondBase, 1, "WiegelmannProcCondBase");
        CbinDataStruct::Write(fid, &WiegelmannProcCondBase2, 1, "WiegelmannProcCondBase2");

        int N[3];
        field->dimensions(N);
        CbinDataStruct::Write(fid, N, 3, "N");
        CagmScalarField comp(N);

        field->getComponent(&comp, PLANE_X);
        double v = comp.sum();
        CbinDataStruct::Write(fid, &v, 1, "FieldXSum");
        double *a = field->getAddress(PLANE_X, 5, 7, 0);
        CbinDataStruct::Write(fid, a, 1, "X570");
        a = field->getAddress(PLANE_X, 7, 5, 0);
        CbinDataStruct::Write(fid, a, 1, "X750");

        field->getComponent(&comp, PLANE_Y);
        v = comp.sum();
        CbinDataStruct::Write(fid, &v, 1, "FieldYSum");
        a = field->getAddress(PLANE_Y, 5, 7, 0);
        CbinDataStruct::Write(fid, a, 1, "Y570");
        a = field->getAddress(PLANE_Y, 7, 5, 0);
        CbinDataStruct::Write(fid, a, 1, "Y750");

        field->getComponent(&comp, PLANE_Z);
        v = comp.sum();
        CbinDataStruct::Write(fid, &v, 1, "FieldZSum");
        a = field->getAddress(PLANE_Z, 5, 7, 0);
        CbinDataStruct::Write(fid, a, 1, "Z570");
        a = field->getAddress(PLANE_Z, 7, 5, 0);
        CbinDataStruct::Write(fid, a, 1, "Z750");

        v = w->sum();
        CbinDataStruct::Write(fid, &v, 1, "WeightSum");

        CbinDataStruct::WriteFooter(fid);

        fclose(fid);
    }
#endif
}
