
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
void pkn_MVectorSumd ( int m, int n, double *sum, ... )
{
  va_list ap;
  int     i, j;
  double  *x;

  va_start ( ap, sum );
  x = va_arg ( ap, double* );
  for ( j = 0; j < n; j++ )
    sum[j] = x[j];
  for ( i = 1; i < m; i++ ) {
    x = va_arg ( ap, double* );
    for ( j = 0; j < n; j++ )
      sum[j] += x[j];
  }
  va_end ( ap );
} /*pkn_MVectorSumd*/

void pkn_MVectorLinCombd ( int m, int n, double *sum, ... )
{
  va_list ap;
  int     i, j;
  double  *x;
  double  c;

  va_start ( ap, sum );
  x = va_arg ( ap, double* );
  c = va_arg ( ap, double );
  for ( j = 0; j < n; j++ )
    sum[j] = c*x[j];
  for ( i = 1; i < m; i++ ) {
    x = va_arg ( ap, double* );
    c = va_arg ( ap, double );
    for ( j = 0; j < n; j++ )
      sum[j] += c*x[j];
  }
  va_end ( ap );
} /*pkn_MVectorLinCombd*/

