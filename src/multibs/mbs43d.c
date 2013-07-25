
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
/* finding Bezier representations of derivatives of Bezier curves */

void mbs_multiFindBezDerivatived ( int degree, int ncurves, int spdimen,
                                   int pitch, const double *ctlpoints,
                                   int dpitch, double *dctlpoints )
{
  int i;

  for ( i = 0; i <= degree-1; i++ )
    pkn_MatrixMDifferenced ( ncurves, spdimen,
                             pitch, &ctlpoints[(i+1)*spdimen],
                             pitch, &ctlpoints[i*spdimen], (double)degree,
                             dpitch, &dctlpoints[i*spdimen] );
} /*mbs_multiFindBezDerivatived*/

