
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <pthread.h>

#define CONST_

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

#include "raybezprivatef.h"

/* ////////////////////////////////////////////////////////////////////////// */
boolean _rbez_ConvexHullTest1f ( int n, float *ff )
{
  int i;

  if ( ff[0] < 0.0 ) {
    for ( i = 1; i <= n; i++ )
      if ( ff[i] >= 0.0 )
        return true;
  }
  else if (  ff[0] > 0.0 ) {
    for ( i = 1; i <= n; i++ )
      if ( ff[i] <= 0.0 )
        return true;
  }
  else
    return true;
  return false;
} /*_rbez_ConvexHullTest1f*/

boolean _rbez_UniquenessTest1f ( int n, float *ff )
{
  int i, nps, nns, lastsgn, nextsgn;

  nps = nns = 0;
  lastsgn = -2;
  for ( i = 0; i <= n; i++ ) {
    if ( ff[i] < 0.0 )
      nextsgn = -1;
    else
      nextsgn = ff[i] > 0.0;
    if ( nextsgn != lastsgn ) {
      if ( nextsgn >= 0 )
        nps ++;
      if ( nextsgn <= 0 )
        nns ++;
      if ( nps > 1 || nns > 1 )
        return false;
    }
    lastsgn = nextsgn;
  }
  return true;
} /*_rbez_UniquenessTest1f*/

boolean _rbez_NewtonMethod1f ( int n, float *ff, float *z )
{
#define MAXITER 7
#define EPS     1.0e-6
#define DELTA   1.0e-6
  float t, d, f, df;
  int   i;

  t = 0.5;
  for ( i = 0; i < MAXITER; i++ ) {
    mbs_BCHornerDerC1f ( n, ff, t, &f, &df );
    if ( f*f < EPS*EPS ) {
      *z = t;
      return true;
    }
    t -= d = f/df;
    if ( d*d < DELTA*DELTA ) {
      *z = t;
      return true;
    }
    if ( t < 0.0 || t > 1.0 )
      return false;
  }
  return false;
} /*_rbez_NewtonMethod1f*/

