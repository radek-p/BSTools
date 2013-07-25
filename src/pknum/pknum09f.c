
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


/* ///////////////////////////////////////////////////////////////////////// */
void pkn_PrintBandmf ( int ncols, const bandm_profile *aprof, const float *a )
{
  int   i, j, k, nrows;

  nrows = 0;
  for ( j = 0; j < ncols; j++ ) {
    k = aprof[j].firstnz + (aprof[j+1].ind-aprof[j].ind);
    nrows = max ( nrows, k );
  }
  for ( i = 0; i < nrows; i++ ) {
    for ( j = 0; j < ncols; j++ ) {
      k = aprof[j].firstnz + (aprof[j+1].ind-aprof[j].ind);
      if ( i >= aprof[j].firstnz && i < k )
        printf ( "%6.3f ", a[aprof[j].ind+i-aprof[j].firstnz] );
      else
        printf ( " 0.0   " );
    }
    printf ( "\n" );
  }
  printf ( "\n" );
} /*pkn_PrintBandmf*/

void pkn_PrintBandmRowSumf ( int ncols, const bandm_profile *aprof, const float *a )
{
  int   i, j, k, nrows;
  float b, s;

  nrows = 0;
  for ( j = 0; j < ncols; j++ ) {
    k = aprof[j].firstnz + (aprof[j+1].ind-aprof[j].ind);
    nrows = max ( nrows, k );
  }
  printf ( "nrows=%d, ncols=%d\n", nrows, ncols );
  for ( i = 0; i < nrows; i++ ) {
    s = 0.0;
    for ( j = 0; j < ncols; j++ ) {
      k = aprof[j].firstnz + (aprof[j+1].ind-aprof[j].ind);
      if ( i >= aprof[j].firstnz && i < k ) {
        b = a[aprof[j].ind+i-aprof[j].firstnz];
        s += b;
        printf ( "%6.3f ", b );
      }
      else
        printf ( " 0.0   " );
    }
    printf ( ";%6.3f\n", s );
  }
  printf ( "\n" );
} /*pkn_PrintBandmRowSumf*/

/* ///////////////////////////////////////////////////////////////////////// */
void pkn_PrintMatf ( int nrows, int ncols, const float *a )
{
  int i, j, k;

  for ( i = 0, k = 0;  i < nrows;  i++ ) {
    for ( j = 0;  j < ncols;  j++, k++ )
      printf ( "%6.3f ", a[k] );
    printf ( "\n" );
  }
} /*pkn_PrintMatf*/

