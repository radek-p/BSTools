
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
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

/* /////////////////////////////////////////// */
/* adding B-spline curves */

boolean mbs_multiAddBSCurvesd ( int ncurves, int spdimen,
                                int degree1, int lastknot1, CONST_ double *knots1,
                                int pitch1, CONST_ double *ctlpoints1,
                                int degree2, int lastknot2, CONST_ double *knots2,
                                int pitch2, CONST_ double *ctlpoints2,
                                int *sumdeg, int *sumlastknot, double *sumknots,
                                int sumpitch, double *sumctlpoints )
{
  void  *sp;
  int   _sumdeg, _sumlastknot, minpitch;
  double *ctlp1, *ctlp2, *_sumknots;

  sp = pkv_GetScratchMemTop ();

  _sumdeg = 0;
  if ( !mbs_FindBSCommonKnotSequenced ( &_sumdeg, &_sumlastknot, &_sumknots,
             2, degree1, lastknot1, knots1, degree2, lastknot2, knots2 ) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  minpitch = spdimen*(_sumlastknot-_sumdeg);
  ctlp1 = pkv_GetScratchMemd ( ncurves*minpitch );
  ctlp2 = pkv_GetScratchMemd ( ncurves*minpitch );
  if ( !ctlp1 || !ctlp2 ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  if ( !mbs_multiAdjustBSCRepd ( ncurves, spdimen,
           degree1, lastknot1, knots1, pitch1, ctlpoints1,
           _sumdeg, _sumlastknot, _sumknots, minpitch, ctlp1 ) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  if ( !mbs_multiAdjustBSCRepd ( ncurves, spdimen,
           degree2, lastknot2, knots2, pitch2, ctlpoints2,
           _sumdeg, _sumlastknot, _sumknots, minpitch, ctlp2 ) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  pkn_AddMatrixd ( ncurves, minpitch,
                   minpitch, ctlp1, minpitch, ctlp2, sumpitch, sumctlpoints );
  *sumdeg = _sumdeg;
  *sumlastknot = _sumlastknot;
  memcpy ( sumknots, _sumknots, (_sumlastknot+1)*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiAddBSCurvesd*/

/* /////////////////////////////////////////// */
/* subtracting B-spline curves */

boolean mbs_multiSubtractBSCurvesd ( int ncurves, int spdimen,
                             int degree1, int lastknot1, CONST_ double *knots1,
                             int pitch1, CONST_ double *ctlpoints1,
                             int degree2, int lastknot2, CONST_ double *knots2,
                             int pitch2, CONST_ double *ctlpoints2,
                             int *sumdeg, int *sumlastknot, double *sumknots,
                             int sumpitch, double *sumctlpoints )
{
  void  *sp;
  int   _sumdeg, _sumlastknot, minpitch;
  double *ctlp1, *ctlp2, *_sumknots;

  sp = pkv_GetScratchMemTop ();

  _sumdeg = 0;
  if ( !mbs_FindBSCommonKnotSequenced ( &_sumdeg, &_sumlastknot, &_sumknots,
             2, degree1, lastknot1, knots1, degree2, lastknot2, knots2 ) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  minpitch = spdimen*(_sumlastknot-_sumdeg);
  ctlp1 = pkv_GetScratchMemd ( ncurves*minpitch );
  ctlp2 = pkv_GetScratchMemd ( ncurves*minpitch );
  if ( !ctlp1 || !ctlp2 ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  if ( !mbs_multiAdjustBSCRepd ( ncurves, spdimen,
           degree1, lastknot1, knots1, pitch1, ctlpoints1,
           _sumdeg, _sumlastknot, _sumknots, minpitch, ctlp1 ) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  if ( !mbs_multiAdjustBSCRepd ( ncurves, spdimen,
           degree2, lastknot2, knots2, pitch2, ctlpoints2,
           _sumdeg, _sumlastknot, _sumknots, minpitch, ctlp2 ) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  pkn_SubtractMatrixd ( ncurves, minpitch,
                   minpitch, ctlp1, minpitch, ctlp2, sumpitch, sumctlpoints );
  *sumdeg = _sumdeg;
  *sumlastknot = _sumlastknot;
  memcpy ( sumknots, _sumknots, (_sumlastknot+1)*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiSubtractBSCurvesd*/

