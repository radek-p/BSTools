
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

double mbs_GrevilleAbscissad ( int degree, double *knots, int i )
{
  int j;
  double xi;

  xi = knots[i+1];
  for ( j = 2; j <= degree; j++ )
    xi += knots[i+j];
  return xi/(double)degree;
} /*mbs_GrevilleAbscissad*/

