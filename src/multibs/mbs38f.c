
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkgeom.h"

#undef CONST_
#define CONST_

#include "multibs.h"

#include "msgpool.h"

/* /////////////////////////////////////////// */
/* adding B-spline curves */

void mbs_multiAddBSCurvesf ( int ncurves, int spdimen,
                             int degree1, int lastknot1, CONST_ float *knots1,
                             int pitch1, CONST_ float *ctlpoints1,
                             int degree2, int lastknot2, CONST_ float *knots2,
                             int pitch2, CONST_ float *ctlpoints2,
                             int *sumdeg, int *sumlastknot, float *sumknots,
                             int sumpitch, float *sumctlpoints )
{
  void  *sp;
  int   _sumdeg, _sumlastknot, minpitch;
  float *ctlp1, *ctlp2, *_sumknots;

  sp = pkv_GetScratchMemTop ();

  _sumdeg = 0;
  if ( !mbs_FindBSCommonKnotSequencef ( &_sumdeg, &_sumlastknot, &_sumknots,
             2, degree1, lastknot1, knots1, degree2, lastknot2, knots2 ) ) {
    pkv_SignalError ( LIB_MULTIBS, 34, ERRMSG_0 );
    exit ( 1 );
  }
  minpitch = spdimen*(_sumlastknot-_sumdeg);
  ctlp1 = pkv_GetScratchMemf ( ncurves*minpitch );
  ctlp2 = pkv_GetScratchMemf ( ncurves*minpitch );
  if ( !ctlp1 || !ctlp2 ) {
    pkv_SignalError ( LIB_MULTIBS, 35, ERRMSG_0 );
    exit ( 1 );
  }
  if ( !mbs_multiAdjustBSCRepf ( ncurves, spdimen,
           degree1, lastknot1, knots1, pitch1, ctlpoints1,
           _sumdeg, _sumlastknot, _sumknots, minpitch, ctlp1 ) ) {
    pkv_SignalError ( LIB_MULTIBS, 36, ERRMSG_0 );
    exit ( 1 );
  }
  if ( !mbs_multiAdjustBSCRepf ( ncurves, spdimen,
           degree2, lastknot2, knots2, pitch2, ctlpoints2,
           _sumdeg, _sumlastknot, _sumknots, minpitch, ctlp2 ) ) {
    pkv_SignalError ( LIB_MULTIBS, 36, ERRMSG_0 );
    exit ( 1 );
  }
  pkn_AddMatrixf ( ncurves, minpitch,
                   minpitch, ctlp1, minpitch, ctlp2, sumpitch, sumctlpoints );
  *sumdeg = _sumdeg;
  *sumlastknot = _sumlastknot;
  memcpy ( sumknots, _sumknots, (_sumlastknot+1)*sizeof(float) );

  pkv_SetScratchMemTop ( sp );
} /*mbs_multiAddBSCurvesf*/

/* /////////////////////////////////////////// */
/* subtracting B-spline curves */

void mbs_multiSubtractBSCurvesf ( int ncurves, int spdimen,
                             int degree1, int lastknot1, CONST_ float *knots1,
                             int pitch1, CONST_ float *ctlpoints1,
                             int degree2, int lastknot2, CONST_ float *knots2,
                             int pitch2, CONST_ float *ctlpoints2,
                             int *sumdeg, int *sumlastknot, float *sumknots,
                             int sumpitch, float *sumctlpoints )
{
  void  *sp;
  int   _sumdeg, _sumlastknot, minpitch;
  float *ctlp1, *ctlp2, *_sumknots;

  sp = pkv_GetScratchMemTop ();

  _sumdeg = 0;
  if ( !mbs_FindBSCommonKnotSequencef ( &_sumdeg, &_sumlastknot, &_sumknots,
             2, degree1, lastknot1, knots1, degree2, lastknot2, knots2 ) ) {
    pkv_SignalError ( LIB_MULTIBS, 34, ERRMSG_0 );
    exit ( 1 );
  }
  minpitch = spdimen*(_sumlastknot-_sumdeg);
  ctlp1 = pkv_GetScratchMemf ( ncurves*minpitch );
  ctlp2 = pkv_GetScratchMemf ( ncurves*minpitch );
  if ( !ctlp1 || !ctlp2 ) {
    pkv_SignalError ( LIB_MULTIBS, 35, ERRMSG_0 );
    exit ( 1 );
  }
  if ( !mbs_multiAdjustBSCRepf ( ncurves, spdimen,
           degree1, lastknot1, knots1, pitch1, ctlpoints1,
           _sumdeg, _sumlastknot, _sumknots, minpitch, ctlp1 ) ) {
    pkv_SignalError ( LIB_MULTIBS, 36, ERRMSG_0 );
    exit ( 1 );
  }
  if ( !mbs_multiAdjustBSCRepf ( ncurves, spdimen,
           degree2, lastknot2, knots2, pitch2, ctlpoints2,
           _sumdeg, _sumlastknot, _sumknots, minpitch, ctlp2 ) ) {
    pkv_SignalError ( LIB_MULTIBS, 36, ERRMSG_0 );
    exit ( 1 );
  }
  pkn_SubtractMatrixf ( ncurves, minpitch,
                   minpitch, ctlp1, minpitch, ctlp2, sumpitch, sumctlpoints );
  *sumdeg = _sumdeg;
  *sumlastknot = _sumlastknot;
  memcpy ( sumknots, _sumknots, (_sumlastknot+1)*sizeof(float) );

  pkv_SetScratchMemTop ( sp );
} /*mbs_multiSubtractBSCurvesf*/

