
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
#include "multibs.h"


/* /////////////////////////////////////////// */
/* Division of Bezier curves with the de Casteljau algorithm */

void mbs_multiBisectBezCurvesf ( int degree, int ncurves,
                                 int spdimen, int pitch,
                                 float *ctlp, float *ctlq )
{
  int i, j, k, l, m;
  float *p, *q;

  p = ctlp;
  q = ctlq;
  for ( k = 0; k < ncurves; k++ ) {
    memcpy ( q, p, spdimen*sizeof(float) );
    for ( j = 1; j <= degree; j++ ) {
      for ( i = 0, m = 0; i <= degree-j; i++, m += spdimen )
        for ( l = 0; l < spdimen; l++ )
          p[m+l] = (float)(0.5*(p[m+l]+p[m+spdimen+l]));
      memcpy ( &q[j*spdimen], p, spdimen*sizeof(float) );
    }
    p += pitch;
    q += pitch;
  }
} /*mbs_multiBisectBezCurvesf*/

void mbs_multiDivideBezCurvesf ( int degree, int ncurves,
                                 int spdimen, int pitch,
                                 float t,
                                 float *ctlp, float *ctlq )
{
  int i, j, k, l, m;
  float *p, *q;
  float s;

  s = (float)(1.0-t);
  p = ctlp;
  q = ctlq;
  for ( k = 0; k < ncurves; k++ ) {
    memcpy ( q, p, spdimen*sizeof(float) );
    for ( j = 1; j <= degree; j++ ) {
      for ( i = 0, m = 0; i <= degree-j; i++, m += spdimen )
        for ( l = 0; l < spdimen; l++ )
          p[m+l] = s*p[m+l] + t*p[m+spdimen+l];
      memcpy ( &q[j*spdimen], p, spdimen*sizeof(float) );
    }
    p += pitch;
    q += pitch;
  }
} /*mbs_multiDivideBezCurvesf*/

