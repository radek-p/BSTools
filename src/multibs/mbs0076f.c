
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

static int FindMaxInt ( int a, int b, int c, int d, int e, int f, int g )
{
  if ( b > a ) a = b;  if ( c > a ) a = c;  if ( d > a ) a = d;
  if ( e > a ) a = e;  if ( f > a ) a = f;
  return g > a ? g : a;
} /*FindMaxInt*/

boolean mbs_BezC2CoonsToBezf ( int spdimen,
                               int degc00, const float *c00,
                               int degc01, const float *c01,
                               int degc02, const float *c02,
                               int degc10, const float *c10,
                               int degc11, const float *c11,
                               int degc12, const float *c12,
                               int degd00, const float *d00,
                               int degd01, const float *d01,
                               int degd02, const float *d02,
                               int degd10, const float *d10,
                               int degd11, const float *d11,
                               int degd12, const float *d12,
                               int *n, int *m, float *p )
{
  void    *sp;
  int     degu, degv;
  float   *workspace;
  boolean result;


  degu = FindMaxInt ( degc00, degc01, degc02, degc10, degc11, degc12, 5 );
  degv = FindMaxInt ( degd00, degd01, degd02, degd10, degd11, degd12, 5 );
  sp = workspace = pkv_GetScratchMemf ( (39 + 3*(degu+1)*(degv+1) +
                                         6*(max(degu,degv)+1))*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_BezC2CoonsToBezf ( spdimen,
                   degc00, c00, degc01, c01, degc02, c02,
                   degc10, c10, degc11, c11, degc12, c12,
                   degd00, d00, degd01, d01, degd02, d02,
                   degd10, d10, degd11, d11, degd12, d12,
                   n, m, p, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BezC2CoonsToBezf*/

