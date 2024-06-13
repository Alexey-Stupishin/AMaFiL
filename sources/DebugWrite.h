#pragma once

#include "stdDefinitions.h"
#include "agmScalarField.h"
#include "agmVectorField.h"
#include "debug_data_trace_win.h"

void DebugWritePars(const char *fname, CagmVectorField * field, CagmScalarField * w);
void DebugWriteData(CubeXD *v, const char *fname, int depth = 0, int iter = 0);
