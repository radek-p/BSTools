
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"

#undef CONST_
#define CONST_
#include "multibs.h"

boolean mbs_FindBSCommonKnotSequenced ( int *degree, int *lastknot,
            double **knots, int nsequences, ... )
{
  void    *sp;
  va_list ap;
  int     i, j, k, l, r;
  int     deg, nki, maxk, degc, lknc;
  double   *knotsc, *kn, a, b;
  int     *knmult;

  sp = pkv_GetScratchMemTop ();
  if ( nsequences < 1 )
    return false;

        /* find the degree and allocate memory */
  degc = lknc = 0;  knotsc = NULL;
  va_start ( ap, nsequences );
  deg = *degree;  /* deg must not be less than the initial value of *degree */
  nki = 0;
  for ( i = 0; i < nsequences; i++ ) {
    degc = va_arg ( ap, int );
    lknc = va_arg ( ap, int );
    knotsc = va_arg ( ap, double* );
    deg = max ( deg, degc );
    nki += mbs_NumKnotIntervalsd ( degc, lknc, knotsc );
  }
  *degree = deg;
  a = knotsc[degc];  b = knotsc[lknc-degc];
  va_end ( ap );

        /* allocate an array and copy the knots */
  *knots = kn = pkv_GetScratchMemd ( (maxk = nki*deg+2) );
  if ( !kn ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  va_start ( ap, nsequences );
  for ( i = k = 0; i < nsequences; i++ ) {
    degc = va_arg ( ap, int );
    lknc = va_arg ( ap, int );
    knotsc = va_arg ( ap, double* );
    for ( j = degc; j <= lknc-degc; j++ ) {
      if ( k > maxk ) {
        PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_6, ERRMSG_6 );
        goto failure;
      }
      kn[k++] = knotsc[j];
    }
        /* verify whether all curves have the same domain */
    if ( knotsc[degc] != a || knotsc[lknc-degc] != b ) {
      va_end ( ap );
      goto failure;
    }
  }
  va_end ( ap );
        /* remove multiple knots, then sort and remove again */
  for ( j = k, i = k = 1; i < j; i++ )
    if ( kn[i] != kn[i-1] )
      kn[k++] = kn[i];
  if ( pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double),
                      0, k, kn ) != SORT_OK )
    goto failure;
  for ( j = k, i = k = 1; i < j; i++ )
    if ( kn[i] != kn[i-1] )
      kn[k++] = kn[i];

        /* now find the proper knot multiplicities */
  if ( !(knmult = pkv_GetScratchMem ( k*sizeof(int) )) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  for ( j = 0; j < k; j++ )
    knmult[j] = 1;
  knmult[0] = knmult[k-1] = deg+1;
  va_start ( ap, nsequences );
  for ( i = 0; i < nsequences; i++ ) {
    degc = va_arg ( ap, int );
    lknc = va_arg ( ap, int );
    knotsc = va_arg ( ap, double* );
    for ( j = l = 0;  j < k;  j++ ) {
      while ( knotsc[l] < kn[j] )
        l ++;
      r = 0;
      while ( l+r <= lknc && r <= degc && knotsc[l+r] == kn[j] )
        r ++;
      if ( r ) {
        r += deg-degc;
        knmult[j] = max ( knmult[j], r );
      }
    }
  }
  va_end ( ap );
        /* and setup the final knot sequence */
  memmove ( &kn[maxk-k], kn, k*sizeof(double) );
  for ( i = j = 0, l = maxk-k;  i < k;  i++, l++ ) {
    for ( r = 0; r < knmult[i]; r++ )
      kn[j++] = kn[l];
  }
  *lastknot = j-1;
  pkv_SetScratchMemTop ( &kn[j] );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_FindBSCommonKnotSequenced*/

boolean mbs_multiAdjustBSCRepd ( int ncurves, int spdimen,
     int indegree, int inlastknot, const double *inknots,
     int inpitch, const double *inctlpoints,
     int outdegree, int outlastknot, CONST_ double *outknots,
     int outpitch, double *outctlpoints )
{
  void  *sp;
  int   kni, lknt, ptch;
  double *knt, *ctlp;

  sp = pkv_GetScratchMemTop ();
  if ( outdegree < indegree )
    return false;

  if ( outdegree > indegree ) {
        /* degree elevation */
    kni = mbs_NumKnotIntervalsd ( indegree, inlastknot, inknots );
    lknt = inlastknot+(kni+1)*(outdegree-indegree);
    ptch = spdimen*(lknt-indegree);
    knt = pkv_GetScratchMemd ( lknt+1 );
    ctlp = pkv_GetScratchMemd ( ncurves*ptch );
    if ( !knt || !ctlp )
      goto failure;
    mbs_multiBSDegElevd ( ncurves, spdimen, indegree, inlastknot, inknots,
                          inpitch, inctlpoints, outdegree-indegree,
                          &outdegree, &lknt, knt, ptch, ctlp, false );
  }
  else {
        /* just copy data */
    lknt = inlastknot;
    ptch = spdimen*(lknt-indegree);
    knt = pkv_GetScratchMemd ( lknt+1 );
    ctlp = pkv_GetScratchMemd ( ncurves*ptch );
    if ( !knt || !ctlp )
      goto failure;
    memcpy ( knt, inknots, (lknt+1)*sizeof(double) );
    pkv_Selectd ( ncurves, ptch, inpitch, ptch, inctlpoints, ctlp );
  }
        /* set the end knots as in the final sequence */
  mbs_multiBSChangeLeftKnotsd ( ncurves, spdimen, outdegree, knt,
                                ptch, ctlp, outknots );
  mbs_multiBSChangeRightKnotsd ( ncurves, spdimen, outdegree, lknt, knt,
                                 ptch, ctlp, &outknots[outlastknot-outdegree] );
        /* it is assumed that now the knot sequence of the curve */
        /* is a subsequence of the final knot sequence */
  mbs_multiOsloInsertKnotsd ( ncurves, spdimen, outdegree,
                              lknt, knt, ptch, ctlp,
                              outlastknot, outknots, outpitch, outctlpoints );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiAdjustBSCRepd*/

