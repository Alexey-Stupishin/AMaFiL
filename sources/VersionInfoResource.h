#pragma once

#define VIR_QUOTE(a) #a
#define VIR_QUOTE_SUBST(a) VIR_QUOTE(a)

//-------------------------------------------------------------------------
#define VIR_Ver1 4
#define VIR_Ver2 2
#define VIR_Ver3 25
#define VIR_Ver4 326
#define VIR_Revision 26
#define VIR_Year 2025

//-------------------------------------------------------------------------
#define VIR_CompanyName "St. Petersburg State University, Russia"
#define VIR_InternalName "MagFieldOps.dll"
#define VIR_OriginalFilename "WWNLFFFReconstruction.dll"
#define VIR_ProductName "Weighted Wiegelmann NLFFF Reconstruction Library"

#ifdef _WINDOWS
#define VIR_T_REV rev.
#define VIR_T_LIBNAME : Weighted Wiegelmann NLFFF Reconstruction Library
#define VIR_COPYRIGHT Copyright (C) Alexey G. Stupishin (agstup@yandex.ru)
#define VIR_FROM , 2017-
#else
#define VIR_T_REV "rev."
#define VIR_T_LIBNAME ": Weighted Wiegelmann NLFFF Reconstruction Library"
#define VIR_COPYRIGHT "Copyright (C) Alexey G. Stupishin (agstup@yandex.ru)"
#define VIR_FROM ", 2017-"
#endif
