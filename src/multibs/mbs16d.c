
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

void mbs_multiBisectBezCurvesd ( int degree, int ncurves,
                                 int spdimen, int pitch,
                                 double *ctlp, double *ctlq )
{
  int i, j, k, l, m;
  double *p, *q;

  p = ctlp;
  q = ctlq;
  for ( k = 0; k < ncurves; k++ ) {
    memcpy ( q, p, spdimen*sizeof(double) );
    for ( j = 1; j <= degree; j++ ) {
      for ( i = 0, m = 0; i <= degree-j; i++, m += spdimen )
        for ( l = 0; l < spdimen; l++ )
          p[m+l] = 0.5*(p[m+l]+p[m+spdimen+l]);
      memcpy ( &q[j*spdimen], p, spdimen*sizeof(double) );
    }
    p += pitch;
    q += pitch;
  }
} /*mbs_multiBisectBezCurvesd*/

void mbs_multiDivideBezCurvesd ( int degree, int ncurves,
                                 int spdimen, int pitch,
                                 double t,
                                 double *ctlp, double *ctlq )
{
  int i, j, k, l, m;
  double *p, *q;
  double s;

  s = 1.0-t;
  p = ctlp;
  q = ctlq;
  for ( k = 0; k < ncurves; k++ ) {
    memcpy ( q, p, spdimen*sizeof(double) );
    for ( j = 1; j <= degree; j++ ) {
      for ( i = 0, m = 0; i <= degree-j; i++, m += spdimen )
        for ( l = 0; l < spdimen; l++ )
          p[m+l] = s*p[m+l] + t*p[m+spdimen+l];
      memcpy ( &q[j*spdimen], p, spdimen*sizeof(double) );
    }
    p += pitch;
    q += pitch;
  }
} /*mbs_multiDivideBezCurvesd*/

