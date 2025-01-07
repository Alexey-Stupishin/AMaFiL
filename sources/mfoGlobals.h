#pragma once

#include "stdDefinitions.h"
#include <map>
#include <string>
#include "agmScalarField.h"
#include "agmVectorFieldOps.h"
#include "TaskQueueProcessor.h"

#define DEBUG_OUT_PATH "g:\\temp\\"

#ifndef d_mfoglobaldefine
#define d_mfoglobaldefineextern extern
#else
#define d_mfoglobaldefineextern
#endif

#define LIB_STATE_LIB_NOT_INIT              0x80010000
#define LIB_STATE_NO_PARAMETER              0x80020000

//d_mfoglobaldefineextern BWndDebug *debugtrace;
//d_mfoglobaldefineextern BWndDebug *debugtrace2;

d_mfoglobaldefineextern int gstcLibraryInit
#ifdef d_mfoglobaldefine
    = 0
#endif
    ;

d_mfoglobaldefineextern std::map<std::string, int> mapInt;
d_mfoglobaldefineextern std::map<std::string, uint64_t> mapuint64_t;
d_mfoglobaldefineextern std::map<std::string, double> mapDouble;

d_mfoglobaldefineextern int WiegelmannWeightType;
d_mfoglobaldefineextern double WiegelmannWeightBound;
d_mfoglobaldefineextern double WiegelmannWeightDivfree;
d_mfoglobaldefineextern int WiegelmannBoundsCorrection;

d_mfoglobaldefineextern double WiegelmannInversionTolerance;
d_mfoglobaldefineextern double WiegelmannInversionDenom;

d_mfoglobaldefineextern int WiegelmannDerivStencil;

d_mfoglobaldefineextern double WiegelmannProcStep0;
d_mfoglobaldefineextern int WiegelmannProcStepMax;
d_mfoglobaldefineextern int WiegelmannProcMaxSteps;

d_mfoglobaldefineextern double WiegelmannProcStepIncr;
d_mfoglobaldefineextern double WiegelmannProcStepIncrInit;
d_mfoglobaldefineextern double WiegelmannProcStepIncrMatr;
d_mfoglobaldefineextern double WiegelmannProcStepIncrMain;

d_mfoglobaldefineextern double WiegelmannProcStepDecr;
d_mfoglobaldefineextern double WiegelmannProcStepDecrInit;
d_mfoglobaldefineextern double WiegelmannProcStepDecrMatr;
d_mfoglobaldefineextern double WiegelmannProcStepDecrMain;

d_mfoglobaldefineextern double WiegelmannProcStepLim;
d_mfoglobaldefineextern double WiegelmannProcStepLimInit;
d_mfoglobaldefineextern double WiegelmannProcStepLimMatr;
d_mfoglobaldefineextern double WiegelmannProcStepLimMain;

d_mfoglobaldefineextern int WiegelmannProcdLStdWin;
d_mfoglobaldefineextern int WiegelmannProcdLStdWinInit;
d_mfoglobaldefineextern int WiegelmannProcdLStdWinMatr;
d_mfoglobaldefineextern int WiegelmannProcdLStdWinMain;

d_mfoglobaldefineextern double WiegelmannProcdLStdVal;
d_mfoglobaldefineextern double WiegelmannProcdLStdValInit;
d_mfoglobaldefineextern double WiegelmannProcdLStdValMatr;
d_mfoglobaldefineextern double WiegelmannProcdLStdValMain;

d_mfoglobaldefineextern int WiegelmannGetMetricsTheta;

d_mfoglobaldefineextern int WiegelmannMatryoshkaUse;
d_mfoglobaldefineextern int WiegelmannMatryoshkaDeepMinN;
d_mfoglobaldefineextern double WiegelmannMatryoshkaFactor;

d_mfoglobaldefineextern int WiegelmannProcCondType;
d_mfoglobaldefineextern int WiegelmannProcCondAbs;
d_mfoglobaldefineextern int WiegelmannProcCondAbs2;
d_mfoglobaldefineextern int WiegelmannProcCondLOS2;
d_mfoglobaldefineextern int WiegelmannProcCondLOS;
d_mfoglobaldefineextern int WiegelmannProcCondBase;
d_mfoglobaldefineextern int WiegelmannProcCondBase2;

d_mfoglobaldefineextern int CommonThreadsN;
d_mfoglobaldefineextern w_priority WiegelmannThreadPriority;

d_mfoglobaldefineextern int WiegelmannProtocolStep;
d_mfoglobaldefineextern int debug_input;

#undef d_mfoglobaldefineextern
