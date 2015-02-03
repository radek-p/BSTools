
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"

#undef CONST_
#define CONST_

#include "multibs.h"

boolean _mbs_BCHornerDer3Pf ( int degreeu, int degreev, int spdimen,
                              CONST_ float *ctlpoints,
                              float u, float v,
                              float *p, float *pu, float *pv,
                              float *puu, float *puv, float *pvv,
                              float *puuu, float *puuv, float *puvv, float *pvvv,
                              float *workspace )
{
  float *dfa, *dfb, *dfc, *dfd;
  float s, t;

  if ( degreeu < 3 || degreev < 3 ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    return false;
  }

  dfa = workspace;
  dfb = &dfa[16*spdimen];
  dfc = &dfb[12*spdimen];
  dfd = &dfc[8*spdimen];
  workspace = &dfd[4*spdimen];

  if ( !_mbs_FindBezPatchDiagFormf ( degreeu, degreev, spdimen, ctlpoints,
                                     3, 3, u, v, dfa, workspace ) )
    return false;
  s = 1.0-u;  t = 1.0-v;

  pkn_MatrixLinCombf ( 4, 3*spdimen, 4*spdimen, dfa, t,
                       4*spdimen, &dfa[spdimen], v, 3*spdimen, dfb );
  pkn_MatrixLinCombf ( 4, 2*spdimen, 3*spdimen, dfb, t,
                       3*spdimen, &dfb[spdimen], v, 2*spdimen, dfc );
  pkn_MatrixLinCombf ( 4, spdimen, 2*spdimen, dfc, t,
                       2*spdimen, &dfc[spdimen], v, spdimen, dfd );
  pkn_SubtractMatrixf ( 1, spdimen, 0, &dfd[3*spdimen], 0, dfd, 0, puuu );
  pkn_SubtractMatrixf ( 1, spdimen, 0, &dfd[2*spdimen], 0, &dfd[spdimen],
                        0, dfd );
  pkn_AddMatrixMf ( 1, spdimen, 0, puuu, 0, dfd, -3.0, 0, puuu );
  pkn_MultMatrixNumf ( 1, spdimen, 0, puuu,
                       (float)(degreeu*(degreeu-1)*(degreeu-2)), 0, puuu );

  pkn_MatrixLinCombf ( 1, 12*spdimen, 0, dfa, s, 0, &dfa[4*spdimen], u,
                       0, dfb );
  pkn_MatrixLinCombf ( 3, 3*spdimen, 4*spdimen, dfb, t,
                       4*spdimen, &dfb[spdimen], v, 3*spdimen, dfa );
  pkn_MatrixLinCombf ( 3, 2*spdimen, 3*spdimen, dfa, t,
                       3*spdimen, &dfa[spdimen], v, 2*spdimen, dfc );
  pkn_AddMatrixf ( 1, 2*spdimen, 0, dfc, 0, &dfc[4*spdimen], 0, dfc );
  pkn_AddMatrixMf ( 1, 2*spdimen, 0, dfc, 0, &dfc[2*spdimen], -2.0, 0, dfc );
  pkn_SubtractMatrixf ( 1, spdimen, 0, &dfc[spdimen], 0, dfc, 0, puuv );
  pkn_MultMatrixNumf ( 1, spdimen, 0, puuv,
                       (float)(degreeu*(degreeu-1)*degreev), 0, puuv );
  pkn_MatrixLinCombf ( 1, spdimen, 0, dfc, t, 0, &dfc[spdimen], v, 0, puu );
  pkn_MultMatrixNumf ( 1, spdimen, 0, puu,
                       (float)(degreeu*(degreeu-1)), 0, puu );

  pkn_MatrixLinCombf ( 1, 8*spdimen, 0, dfb, s, 0, &dfb[4*spdimen], u,
                       0, dfc );
  pkn_MatrixLinCombf ( 2, 3*spdimen, 4*spdimen, dfc, t,
                       4*spdimen, &dfc[spdimen], v, 3*spdimen, dfa );
  pkn_SubtractMatrixf ( 1, 3*spdimen, 0, &dfa[3*spdimen], 0, dfa, 0, dfa );
  pkn_AddMatrixf ( 1, spdimen, 0, dfa, 0, &dfa[2*spdimen], 0, puvv );
  pkn_AddMatrixMf ( 1, spdimen, 0, puvv, 0, &dfa[spdimen], -2.0, 0, puvv );
  pkn_MultMatrixNumf ( 1, spdimen, 0, puvv,
                       (float)(degreeu*degreev*(degreev-1)), 0, puvv );
  pkn_MatrixLinCombf ( 1, 2*spdimen, 0, dfa, t, 0, &dfa[spdimen], v, 0, dfb );
  pkn_SubtractMatrixf ( 1, spdimen, 0, &dfb[spdimen], 0, dfb, 0, puv );
  pkn_MultMatrixNumf ( 1, spdimen, 0, puv, (float)(degreeu*degreev), 0, puv );
  pkn_MatrixLinCombf ( 1, spdimen, 0, dfb, t, 0, &dfb[spdimen], v, 0, pu );
  pkn_MultMatrixNumf ( 1, spdimen, 0, pu, (float)degreeu, 0, pu );

  pkn_MatrixLinCombf ( 1, 4*spdimen, 0, dfc, s, 0, &dfc[4*spdimen], u,
                       0, dfd );
  pkn_SubtractMatrixf ( 1, spdimen, 0, &dfd[3*spdimen], 0, dfd, 0, pvvv );
  pkn_SubtractMatrixf ( 1, spdimen, 0, &dfd[2*spdimen], 0, &dfd[spdimen],
                        0, dfa );
  pkn_AddMatrixMf ( 1, spdimen, 0, pvvv, 0, dfa, -3.0, 0, pvvv );
  pkn_MultMatrixNumf ( 1, spdimen, 0, pvvv,
                       (float)(degreev*(degreev-1)*(degreev-2)), 0, pvvv );
  pkn_MatrixLinCombf ( 1, 3*spdimen, 0, dfd, t, 0, &dfd[spdimen], v, 0, dfa );
  pkn_AddMatrixf ( 1, spdimen, 0, dfa, 0, &dfa[2*spdimen], 0, pvv );
  pkn_AddMatrixMf ( 1, spdimen, 0, pvv, 0, &dfa[spdimen], -2.0, 0, pvv );
  pkn_MultMatrixNumf ( 1, spdimen, 0, pvv, (float)(degreev*(degreev-1)),
                       0, pvv );
  pkn_MatrixLinCombf ( 1, 2*spdimen, 0, dfa, t, 0, &dfa[spdimen], v, 0, dfb );
  pkn_SubtractMatrixf ( 1, spdimen, 0, &dfb[spdimen], 0, dfb, 0, pv );
  pkn_MultMatrixNumf ( 1, spdimen, 0, pv, (float)degreev, 0, pv );
  pkn_MatrixLinCombf ( 1, spdimen, 0, dfb, t, 0, &dfb[spdimen], v, 0, p );
  return true;
} /*_mbs_BCHornerDer3Pf*/

