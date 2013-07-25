
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"


double pkn_SecondNormd ( int spdimen, const double *b )
{
  int    i;
  double s;

  s = 0.0;
  for ( i = 0; i < spdimen; i++ )
    s += b[i]*b[i];
  return sqrt ( s );
} /*pkn_SecondNormd*/

