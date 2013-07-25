
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

#include "msgpool.h"

/* /////////////////////////////////////////// */

void mbs_FindBezPatchDiagFormd ( int degreeu, int degreev, int spdimen,
                                 CONST_ double *cpoints,
                                 int k, int l, double u, double v,
                                 double *dfcp )
{
  void  *sp;
  double *p;
  int   i, pitch, ptch;

  if ( k > degreeu || l > degreev )
    pkv_SignalError ( LIB_MULTIBS, 51, ERRMSG_1 );
  sp = pkv_GetScratchMemTop ();

  pitch = (degreev+1)*spdimen;
  if ( k < degreeu ) {
    p = pkv_GetScratchMemd ( (k+1)*pitch );
    if ( !p )
      pkv_SignalError ( LIB_MULTIBS, 52, ERRMSG_0 );
    mbs_multiBCHornerd ( degreeu-k, k+1, pitch, pitch, cpoints, u, p );
  }
  else
    p = cpoints;
  if ( l < degreev ) {
    ptch = (l+1)*spdimen;
    for ( i = 0; i <= k; i++ )
      mbs_multiBCHornerd ( degreev-l, l+1, spdimen, spdimen, &p[i*pitch], v,
                           &dfcp[i*ptch] );
  }
  else
    memcpy ( dfcp, p, (k+1)*(l+1)*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
} /*mbs_FindBezPatchDiagFormd*/

void mbs_BCHornerDer3Pd ( int degreeu, int degreev, int spdimen,
                          CONST_ double *ctlpoints,
                          double u, double v,
                          double *p, double *pu, double *pv,
                          double *puu, double *puv, double *pvv,
                          double *puuu, double *puuv, double *puvv, double *pvvv )
{
  void  *sp;
  double *dfa, *dfb, *dfc, *dfd;
  double s, t;

  if ( degreeu < 3 || degreev < 3 )
    pkv_SignalError ( LIB_MULTIBS, 53, ERRMSG_1 );

  sp = pkv_GetScratchMemTop ();
  dfa = pkv_GetScratchMemd ( 16*spdimen );
  dfb = pkv_GetScratchMemd ( 12*spdimen );
  dfc = pkv_GetScratchMemd ( 8*spdimen );
  dfd = pkv_GetScratchMemd ( 4*spdimen );
  if ( !dfa || !dfb || !dfc || !dfd )
    pkv_SignalError ( LIB_MULTIBS, 54, ERRMSG_0 );

  mbs_FindBezPatchDiagFormd ( degreeu, degreev, spdimen, ctlpoints,
                              3, 3, u, v, dfa );
  s = 1.0-u;  t = 1.0-v;

  pkn_MatrixLinCombd ( 4, 3*spdimen, 4*spdimen, dfa, t,
                       4*spdimen, &dfa[spdimen], v, 3*spdimen, dfb );
  pkn_MatrixLinCombd ( 4, 2*spdimen, 3*spdimen, dfb, t,
                       3*spdimen, &dfb[spdimen], v, 2*spdimen, dfc );
  pkn_MatrixLinCombd ( 4, spdimen, 2*spdimen, dfc, t,
                       2*spdimen, &dfc[spdimen], v, spdimen, dfd );
  pkn_SubtractMatrixd ( 1, spdimen, 0, &dfd[3*spdimen], 0, dfd, 0, puuu );
  pkn_SubtractMatrixd ( 1, spdimen, 0, &dfd[2*spdimen], 0, &dfd[spdimen],
                        0, dfd );
  pkn_AddMatrixMd ( 1, spdimen, 0, puuu, 0, dfd, -3.0, 0, puuu );
  pkn_MultMatrixNumd ( 1, spdimen, 0, puuu,
                       (double)(degreeu*(degreeu-1)*(degreeu-2)), 0, puuu );


  pkn_MatrixLinCombd ( 1, 12*spdimen, 0, dfa, s, 0, &dfa[4*spdimen], u,
                       0, dfb );
  pkn_MatrixLinCombd ( 3, 3*spdimen, 4*spdimen, dfb, t,
                       4*spdimen, &dfb[spdimen], v, 3*spdimen, dfa );
  pkn_MatrixLinCombd ( 3, 2*spdimen, 3*spdimen, dfa, t,
                       3*spdimen, &dfa[spdimen], v, 2*spdimen, dfc );
  pkn_AddMatrixd ( 1, 2*spdimen, 0, dfc, 0, &dfc[4*spdimen], 0, dfc );
  pkn_AddMatrixMd ( 1, 2*spdimen, 0, dfc, 0, &dfc[2*spdimen], -2.0, 0, dfc );
  pkn_SubtractMatrixd ( 1, spdimen, 0, &dfc[spdimen], 0, dfc, 0, puuv );
  pkn_MultMatrixNumd ( 1, spdimen, 0, puuv,
                       (double)(degreeu*(degreeu-1)*degreev), 0, puuv );
  pkn_MatrixLinCombd ( 1, spdimen, 0, dfc, t, 0, &dfc[spdimen], v, 0, puu );
  pkn_MultMatrixNumd ( 1, spdimen, 0, puu,
                       (double)(degreeu*(degreeu-1)), 0, puu );


  pkn_MatrixLinCombd ( 1, 8*spdimen, 0, dfb, s, 0, &dfb[4*spdimen], u,
                       0, dfc );
  pkn_MatrixLinCombd ( 2, 3*spdimen, 4*spdimen, dfc, t,
                       4*spdimen, &dfc[spdimen], v, 3*spdimen, dfa );
  pkn_SubtractMatrixd ( 1, 3*spdimen, 0, &dfa[3*spdimen], 0, dfa, 0, dfa );
  pkn_AddMatrixd ( 1, spdimen, 0, dfa, 0, &dfa[2*spdimen], 0, puvv );
  pkn_AddMatrixMd ( 1, spdimen, 0, puvv, 0, &dfa[spdimen], -2.0, 0, puvv );
  pkn_MultMatrixNumd ( 1, spdimen, 0, puvv,
                       (double)(degreeu*degreev*(degreev-1)), 0, puvv );
  pkn_MatrixLinCombd ( 1, 2*spdimen, 0, dfa, t, 0, &dfa[spdimen], v, 0, dfb );
  pkn_SubtractMatrixd ( 1, spdimen, 0, &dfb[spdimen], 0, dfb, 0, puv );
  pkn_MultMatrixNumd ( 1, spdimen, 0, puv, (double)(degreeu*degreev), 0, puv );
  pkn_MatrixLinCombd ( 1, spdimen, 0, dfb, t, 0, &dfb[spdimen], v, 0, pu );
  pkn_MultMatrixNumd ( 1, spdimen, 0, pu, (double)degreeu, 0, pu );


  pkn_MatrixLinCombd ( 1, 4*spdimen, 0, dfc, s, 0, &dfc[4*spdimen], u,
                       0, dfd );
  pkn_SubtractMatrixd ( 1, spdimen, 0, &dfd[3*spdimen], 0, dfd, 0, pvvv );
  pkn_SubtractMatrixd ( 1, spdimen, 0, &dfd[2*spdimen], 0, &dfd[spdimen],
                        0, dfa );
  pkn_AddMatrixMd ( 1, spdimen, 0, pvvv, 0, dfa, -3.0, 0, pvvv );
  pkn_MultMatrixNumd ( 1, spdimen, 0, pvvv,
                       (double)(degreev*(degreev-1)*(degreev-2)), 0, pvvv );
  pkn_MatrixLinCombd ( 1, 3*spdimen, 0, dfd, t, 0, &dfd[spdimen], v, 0, dfa );
  pkn_AddMatrixd ( 1, spdimen, 0, dfa, 0, &dfa[2*spdimen], 0, pvv );
  pkn_AddMatrixMd ( 1, spdimen, 0, pvv, 0, &dfa[spdimen], -2.0, 0, pvv );
  pkn_MultMatrixNumd ( 1, spdimen, 0, pvv, (double)(degreev*(degreev-1)),
                       0, pvv );
  pkn_MatrixLinCombd ( 1, 2*spdimen, 0, dfa, t, 0, &dfa[spdimen], v, 0, dfb );
  pkn_SubtractMatrixd ( 1, spdimen, 0, &dfb[spdimen], 0, dfb, 0, pv );
  pkn_MultMatrixNumd ( 1, spdimen, 0, pv, (double)degreev, 0, pv );
  pkn_MatrixLinCombd ( 1, spdimen, 0, dfb, t, 0, &dfb[spdimen], v, 0, p );

  pkv_SetScratchMemTop ( sp );
} /*mbs_BCHornerDer3Pd*/

