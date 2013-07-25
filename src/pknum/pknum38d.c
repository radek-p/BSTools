
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
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

boolean pkn_NRBSymFindEigenvalueIntervald ( int n, const int *prof,
                                            double *a, double **row,
                                            double *amin, double *amax )
{
  void   *sp;
  double **r;
  int    i, j;
  double bmin, bmax, gr, b;

  sp = pkv_GetScratchMemTop ();
  if ( !row ) {
    r = pkv_GetScratchMem ( n*sizeof(double*) );
    if ( !r )
      goto failure;
    if ( !pkn_NRBFindRowsd ( n, prof, a, r ) )
      goto failure;
  }
  else
    r = row;

  bmin = bmax = r[0][0];
  for ( i = 0; i < n; i++ ) {
    b = r[i][i];
        /* find the radius of a Gershgorin circle */
    gr = 0.0;
    for ( j = prof[i]; j < i; j++ )
      gr += fabs ( r[i][j] );
    for ( j = i+1; j < n; j++ )
      if ( i >= prof[j] )
        gr += fabs ( r[j][i] );
    bmin = min ( bmin, b-gr );
    bmax = max ( bmax, b+gr );
  }

  *amin = bmin;
  *amax = bmax;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_NRBSymFindEigenvalueIntervald*/

