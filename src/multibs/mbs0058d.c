
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

static int FindMaxInt ( int a, int b, int c, int d, int e )
{
  if ( b > a ) a = b;  if ( c > a ) a = c;  if ( d > a ) a = d;
  return e > a ? e : a;
} /*FindMaxInt*/

boolean mbs_BezC1CoonsToBezd ( int spdimen,
                               int degc00, const double *c00,
                               int degc01, const double *c01,
                               int degc10, const double *c10,
                               int degc11, const double *c11,
                               int degd00, const double *d00,
                               int degd01, const double *d01,
                               int degd10, const double *d10,
                               int degd11, const double *d11,
                               int *n, int *m, double *p )
{
  void    *sp;
  double  *workspace;
  int     degu, degv, size;
  boolean result;

  degu = FindMaxInt ( degc00, degc01, degc10, degc11, 3 );
  degv = FindMaxInt ( degd00, degd01, degd10, degd11, 3 );
  size = (18 + 3*(degu+1)*(degv+1) + 4*(max(degu,degv)+1))*spdimen;
  sp = workspace = pkv_GetScratchMemd ( size );
  if ( !workspace )
    return false;
  result = _mbs_BezC1CoonsToBezd ( spdimen,
                   degc00, c00, degc01, c01, degc10, c10, degc11, c11,
                   degd00, d00, degd01, d01, degd10, d10, degd11, d11,
                   n, m, p, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BezC1CoonsToBezd*/

