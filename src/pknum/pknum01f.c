
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


double pkn_ScalarProductf ( int spdimen, const float *a, const float *b )
{
  double s;
  int    i;

  s = a[0]*b[0];
  for ( i = 1; i < spdimen; i++ )
    s += a[i]*b[i];
  return s;
} /*pkn_ScalarProductf*/

