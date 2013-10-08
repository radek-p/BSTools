
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
#include "multibs.h"

boolean mbs_multiBSDegElevClosedd ( int ncurves, int spdimen,
                         int indegree, int inlastknot, const double *inknots,
                         int inpitch, const double *inctlpoints,
                         int deltadeg,
                         int *outdegree, int *outlastknot,
                         double *outknots, int outpitch, double *outctlpoints )
{
  void  *sp;
  double *nlkn, T, u;
  int   deg, lkn, i, r, K;

  sp = pkv_GetScratchMemTop ();
  T = inknots[inlastknot-indegree]-inknots[indegree];
  if ( !mbs_multiBSDegElevd ( ncurves, spdimen, indegree, inlastknot, inknots,
                              inpitch, inctlpoints, deltadeg, outdegree, &lkn,
                              outknots, outpitch, outctlpoints, false ) )
    goto failure;
  deg = *outdegree;
  if ( !(nlkn = pkv_GetScratchMemd ( deg+1 )) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  r = 2*indegree+1;  r = min ( r, inlastknot );
  r = mbs_KnotMultiplicityd ( r, inknots, (u = inknots[indegree]) );
  *outlastknot = lkn += r;
  K = lkn-2*deg;
  for ( i = deg; i >= 0 && i >= deg-r; i-- )
    nlkn[i] = u;
  for ( ; i >= 0; i-- )
    nlkn[i] = min ( outknots[i+K]-T, u );  /* because of rounding errors */
  mbs_multiBSChangeLeftKnotsd ( ncurves, spdimen, deg, outknots,
                                outpitch, outctlpoints, nlkn );
  for ( i = 0; i < deg; i++ )
    pkv_Moved ( ncurves, spdimen, outpitch, spdimen*K, &outctlpoints[i*spdimen] );
  for ( i = deg+1; i+K <= lkn; i++ )
    outknots[i+K] = outknots[i]+T;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBSDegElevClosedd*/

