
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

boolean mbs_FindBSCommonKnotSequencef ( int *degree, int *lastknot,
            float **knots, int nsequences, ... )
{
  va_list ap;
  int     i, j, k, l, r;
  int     deg, nki, maxk, degc, lknc;
  float   *knotsc, *kn, a, b;
  int     *knmult;

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
    knotsc = va_arg ( ap, float* );
    deg = max ( deg, degc );
    nki += mbs_NumKnotIntervalsf ( degc, lknc, knotsc );
  }
  *degree = deg;
  a = knotsc[degc];  b = knotsc[lknc-degc];
  va_end ( ap );

        /* allocate an array and copy the knots */
  *knots = kn = pkv_GetScratchMemf ( (maxk = nki*deg+2) );
  if ( !kn )
    return false;
  va_start ( ap, nsequences );
  for ( i = k = 0; i < nsequences; i++ ) {
    degc = va_arg ( ap, int );
    lknc = va_arg ( ap, int );
    knotsc = va_arg ( ap, float* );
    for ( j = degc; j <= lknc-degc; j++ ) {
      if ( k > maxk )
        PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_6, ERRMSG_6 );
      kn[k++] = knotsc[j];
    }
        /* verify whether all curves have the same domain */
    if ( knotsc[degc] != a || knotsc[lknc-degc] != b ) {
      va_end ( ap );
      return false;
    }
  }
  va_end ( ap );
        /* remove multiple knots, then sort and remove again */
  for ( j = k, i = k = 1; i < j; i++ )
    if ( kn[i] != kn[i-1] )
      kn[k++] = kn[i];
  if ( pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(float),
                      0, k, kn ) != SORT_OK )
    return false;
  for ( j = k, i = k = 1; i < j; i++ )
    if ( kn[i] != kn[i-1] )
      kn[k++] = kn[i];

        /* now find the proper knot multiplicities */
  if ( !(knmult = pkv_GetScratchMem ( k*sizeof(int) )) )
    return false;
  for ( j = 0; j < k; j++ )
    knmult[j] = 1;
  knmult[0] = knmult[k-1] = deg+1;
  va_start ( ap, nsequences );
  for ( i = 0; i < nsequences; i++ ) {
    degc = va_arg ( ap, int );
    lknc = va_arg ( ap, int );
    knotsc = va_arg ( ap, float* );
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
  memmove ( &kn[maxk-k], kn, k*sizeof(float) );
  for ( i = j = 0, l = maxk-k;  i < k;  i++, l++ ) {
    for ( r = 0; r < knmult[i]; r++ )
      kn[j++] = kn[l];
  }
  *lastknot = j-1;
  pkv_SetScratchMemTop ( &kn[j] );
  return true;
} /*mbs_FindBSCommonKnotSequencef*/

boolean mbs_multiAdjustBSCRepf ( int ncurves, int spdimen,
     int indegree, int inlastknot, const float *inknots,
     int inpitch, const float *inctlpoints,
     int outdegree, int outlastknot, CONST_ float *outknots,
     int outpitch, float *outctlpoints )
{
  void  *sp;
  int   kni, lknt, ptch;
  float *knt, *ctlp;

  sp = pkv_GetScratchMemTop ();
  if ( outdegree < indegree )
    return false;

  if ( outdegree > indegree ) {
        /* degree elevation */
    kni = mbs_NumKnotIntervalsf ( indegree, inlastknot, inknots );
    lknt = inlastknot+(kni+1)*(outdegree-indegree);
    ptch = spdimen*(lknt-indegree);
    knt = pkv_GetScratchMemf ( lknt+1 );
    ctlp = pkv_GetScratchMemf ( ncurves*ptch );
    if ( !knt || !ctlp )
      goto failure;
    mbs_multiBSDegElevf ( ncurves, spdimen, indegree, inlastknot, inknots,
                          inpitch, inctlpoints, outdegree-indegree,
                          &outdegree, &lknt, knt, ptch, ctlp, false );
  }
  else {
        /* just copy data */
    lknt = inlastknot;
    ptch = spdimen*(lknt-indegree);
    knt = pkv_GetScratchMemf ( lknt+1 );
    ctlp = pkv_GetScratchMemf ( ncurves*ptch );
    if ( !knt || !ctlp )
      goto failure;
    memcpy ( knt, inknots, (lknt+1)*sizeof(float) );
    pkv_Selectf ( ncurves, ptch, inpitch, ptch, inctlpoints, ctlp );
  }
        /* set the end knots as in the final sequence */
  mbs_multiBSChangeLeftKnotsf ( ncurves, spdimen, outdegree, knt,
                                ptch, ctlp, outknots );
  mbs_multiBSChangeRightKnotsf ( ncurves, spdimen, outdegree, lknt, knt,
                                 ptch, ctlp, &outknots[outlastknot-outdegree] );
        /* it is assumed that now the knot sequence of the curve */
        /* is a subsequence of the final knot sequence */
  mbs_multiOsloInsertKnotsf ( ncurves, spdimen, outdegree,
                              lknt, knt, ptch, ctlp,
                              outlastknot, outknots, outpitch, outctlpoints );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiAdjustBSCRepf*/

