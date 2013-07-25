
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "pkvaria.h"

#undef CONST_
#define CONST_

#include "pknum.h"

/* ///////////////////////////////////////////////////////////////////////// */
void pkn_MVectorSumf ( int m, int n, float *sum, ... )
{
  va_list ap;
  int     i, j;
  float   *x;

  va_start ( ap, sum );
  x = va_arg ( ap, float* );
  for ( j = 0; j < n; j++ )
    sum[j] = x[j];
  for ( i = 1; i < m; i++ ) {
    x = va_arg ( ap, float* );
    for ( j = 0; j < n; j++ )
      sum[j] += x[j];
  }
  va_end ( ap );
} /*pkn_MVectorSumf*/

void pkn_MVectorLinCombf ( int m, int n, float *sum, ... )
{
  va_list ap;
  int     i, j;
  float   *x;
  double  c;

  va_start ( ap, sum );
  x = va_arg ( ap, float* );
  c = va_arg ( ap, double );
  for ( j = 0; j < n; j++ )
    sum[j] = (float)(c*x[j]);
  for ( i = 1; i < m; i++ ) {
    x = va_arg ( ap, float* );
    c = va_arg ( ap, double );
    for ( j = 0; j < n; j++ )
      sum[j] += (float)(c*x[j]);
  }
  va_end ( ap );
} /*pkn_MVectorLinCombf*/

