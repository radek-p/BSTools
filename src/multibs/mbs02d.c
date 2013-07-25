
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
void mbs_multiReverseBSCurved ( int degree, int lastknot, double *knots,
                                int ncurves, int spdimen, int pitch,
                                double *ctlpoints )
{
  int i;

  if ( knots ) {
    for ( i = 0; i <= lastknot; i++ )
      knots[i] = -knots[i];
    pkv_ReverseMatd ( lastknot+1, 1, 1, knots );
  }
  else
    lastknot = 2*degree+1;

  for ( i = 0;  i < ncurves;  i++, ctlpoints += pitch )
    pkv_ReverseMatd ( lastknot-degree, spdimen, spdimen, ctlpoints );
} /*mbs_multiReverseBSCurved*/

