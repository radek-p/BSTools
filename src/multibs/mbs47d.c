
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


int mbs_SetKnotd ( int lastknot, double *knots,
                   int knotnum, int mult, double t )
{
  int i;

  if ( knotnum < 0 || knotnum > lastknot ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    return -1;
  }
  for ( i = 0; i < mult && knotnum > i; i++ )
    knots[knotnum-i] = t;
  if ( pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double),
                      0, lastknot-1, &knots[1] ) != SORT_OK )
    return -1;
  if ( knots[0] > knots[1] )
    knots[0] = knots[1];
  if ( knots[lastknot] < knots[lastknot-1] )
    knots[lastknot] = knots[lastknot-1];
  knotnum = mbs_FindKnotIntervald ( 0, lastknot, knots, t, NULL );
  return knotnum;
} /*mbs_SetKnotd*/

static void _mbs_AuxSetKCd ( int K, int lastknot, double *knots, double T,
                             int knotnum, int mult, double t )
{
  int i, j;

  for ( i = 1-mult; i <= 0; i++ )
    knots[knotnum+i] = t;
  for ( j = 1; knotnum+j*K-mult <= lastknot; j++ )
    for ( i = 1-mult; i <= 0 && knotnum+j*K+i <= lastknot; i++ )
      knots[knotnum+j*K+i] = t + (double)j*T;
  for ( j = 1; knotnum-j*K >= 0; j++ )
    for ( i = 0; i < mult && knotnum-j*K-i >= 0; i++ )
      knots[knotnum-j*K-i] = t - (double)j*T;
} /*_mbs_AuxSetKCd*/

int mbs_SetKnotClosedd ( int degree, int lastknot, double *knots, double T,
                         int knotnum, int mult, double t )
{
  int    K, kn, i, m;
  double tT;

  K = lastknot-2*degree;
  if ( lastknot <= 3*degree || mult > K ||
       knotnum < 0 || knotnum >= lastknot )
    return -1;

  if ( knotnum-mult+1 < degree )
    { kn = knotnum+K;  tT = t+T; }
  else if ( knotnum > lastknot-degree )
    { kn = knotnum-K;  tT = t-T; }
  else
    { kn = knotnum;  tT = t; }

  while ( kn > 1 && tT < knots[kn-mult] ) {
    m = 1;
    while ( kn-1-m > 0 && knots[kn-1-m] == knots[kn-1] )
      m++;
    _mbs_AuxSetKCd ( K, lastknot, knots, T, kn, m, knots[kn-1] );
    kn -= m;  knotnum -= m;
    if ( kn <= m )
      { kn += K;  tT += T; }
  }
  while ( kn < lastknot-1 && tT > knots[kn+1] ) {
    m = 1;
    while ( kn+1+m < lastknot && knots[kn+1+m] == knots[kn+1] )
      m++;
    _mbs_AuxSetKCd ( K, lastknot, knots, T, kn, m, knots[kn+m] );
    kn += m;  knotnum += m;
    if ( kn >= lastknot )
      { kn -= K;  tT -= T; }
  }

  _mbs_AuxSetKCd ( K, lastknot, knots, T, kn, mult, tT );

  if ( kn != knotnum ) { /* remove rounding errors */
    for ( i = 0; i < mult && knotnum-i >= 0; i++ )
      knots[knotnum-i] = t;
  }
    /* rounding errors may still make the knot sequence non monotone, */
    /* the following statement is supposed to extinguish that */
  for ( i = knotnum-1; i >= 0; i-- )
    if ( knots[i] > knots[i+1] )
      knots[i] = knots[i+1];
  for ( i = knotnum+1; i <= lastknot; i++ )
    if ( knots[i] < knots[i-1] )
      knots[i] = knots[i-1];

  while ( knotnum > mult && knots[knotnum] > t )
    knotnum--;
  while ( knotnum < lastknot-1 && knots[knotnum+1] <= t )
    knotnum++;
  return knotnum;
} /*mbs_SetKnotClosedd*/

