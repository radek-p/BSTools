
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
/* finding B-spline representations of derivatives of B-spline curves */

void mbs_multiFindBSDerivativef ( int degree, int lastknot, const float *knots,
                                  int ncurves, int spdimen,
                                  int pitch, const float *ctlpoints,
                                  int *lastdknot, float *dknots,
                                  int dpitch, float *dctlpoints )
{
  int   i;
  float a;

  if ( lastdknot )
    *lastdknot = lastknot-2;
  if ( dknots )
    memcpy ( dknots, &knots[1], (lastknot-1)*sizeof(float) );

  for ( i = 0; i <= lastknot-degree-2; i++ ) {
    a = knots[i+degree+1]-knots[i+1];
    if ( a ) {
      a = (float)degree/a;
      pkn_MatrixMDifferencef ( ncurves, spdimen,
                               pitch, &ctlpoints[(i+1)*spdimen],
                               pitch, &ctlpoints[i*spdimen], a,
                               dpitch, &dctlpoints[i*spdimen] );
    }
    else
      pkv_ZeroMatf ( ncurves, spdimen, dpitch, &dctlpoints[i*spdimen] );
  }
} /*mbs_multiFindBSDerivativef*/

