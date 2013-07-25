
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005,2007                             */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"


/* ///////////////////////////////////////////////// */
/* rectangular matrix multiplication                 */
/* C = A^T*B */

void pkn_MultTMatrixd ( int nrows_a, int rowlen_a, int pitch_a, const double *a,
                        int rowlen_b, int pitch_b, const double *b,
                        int pitch_c, double *c )
{
  int    i, j, k, ka, kb, ic;
  double s;

  for ( i = ic = 0;  i < rowlen_a;  i++, ic += pitch_c )
    for ( j = 0; j < rowlen_b; j++ ) {
      s = a[i]*b[j];
      for ( k = 1, ka = pitch_a, kb = pitch_b;
            k < nrows_a;
            k++, ka += pitch_a, kb += pitch_b )
        s += a[ka+i]*b[kb+j];
      c[ic+j] = (double)s;
    }
} /*pkn_MultTMatrixd*/

/* ///////////////////////////////////////////////// */
/* rectangular matrix multiplication and addition    */
/* D = C + A^T*B */

void pkn_MultTMatrixAddd ( int nrows_a, int rowlen_a, int pitch_a,
                           const double *a,
                           int rowlen_b, int pitch_b, const double *b,
                           int pitch_c, double *c )
{
  int    i, j, k, ka, kb, ic;
  double s;

  for ( i = ic = 0;  i < rowlen_a;  i++, ic += pitch_c )
    for ( j = 0; j < rowlen_b; j++ ) {
      s = a[i]*b[j];
      for ( k = 1, ka = pitch_a, kb = pitch_b;
            k < nrows_a;
            k++, ka += pitch_a, kb += pitch_b )
        s += a[ka+i]*b[kb+j];
      c[ic+j] += (double)s;
    }
} /*pkn_MultTMatrixAddd*/

/* ///////////////////////////////////////////////// */
/* rectangular matrix multiplication and subtraction */
/* D = C - A^T*B */

void pkn_MultTMatrixSubd ( int nrows_a, int rowlen_a, int pitch_a,
                           const double *a,
                           int rowlen_b, int pitch_b, const double *b,
                           int pitch_c, double *c )
{
  int    i, j, k, ka, kb, ic;
  double s;

  for ( i = ic = 0;  i < rowlen_a;  i++, ic += pitch_c )
    for ( j = 0; j < rowlen_b; j++ ) {
      s = a[i]*b[j];
      for ( k = 1, ka = pitch_a, kb = pitch_b;
            k < nrows_a;
            k++, ka += pitch_a, kb += pitch_b )
        s += a[ka+i]*b[kb+j];
      c[ic+j] -= (double)s;
    }
} /*pkn_MultTMatrixSubd*/

/* ///////////////////////////////////////////////// */
/* rectangular matrix multiplication                 */
/* C = A*B^T */

void pkn_MultMatrixTd ( int nrows_a, int rowlen_a, int pitch_a, const double *a,
                        int nrows_b, int pitch_b, const double *b,
                        int pitch_c, double *c )
{
  int    i, j, k, ia, jb, ic;
  double s;

  for ( i = ia = ic = 0;  i < nrows_a;  i++, ia += pitch_a, ic += pitch_c )
    for ( j = jb = 0;  j < nrows_b;  j++, jb += pitch_b ) {
      s = a[ia]*b[jb];
      for ( k = 1; k < rowlen_a; k++ )
        s += a[ia+k]*b[jb+k];
      c[ic+j] = (double)s;
    }
} /*pkn_MultMatrixTd*/

/* ///////////////////////////////////////////////// */
/* rectangular matrix multiplication and addition    */
/* D = C + A*B^T */

void pkn_MultMatrixTAddd ( int nrows_a, int rowlen_a, int pitch_a, const double *a,
                           int nrows_b, int pitch_b, const double *b,
                           int pitch_c, double *c )
{
  int    i, j, k, ia, jb, ic;
  double s;

  for ( i = ia = ic = 0;
        i < nrows_a;
        i++, ia += pitch_a, ic += pitch_c )
    for ( j = jb = 0;  j < nrows_b;  j++, jb += pitch_b ) {
      s = a[ia]*b[jb];
      for ( k = 1; k < rowlen_a; k++ )
        s += a[ia+k]*b[jb+k];
      c[ic+j] += (double)s;
    }
} /*pkn_MultMatrixTAddd*/

/* ///////////////////////////////////////////////// */
/* rectangular matrix multiplication and subtraction */
/* D = C - A*B^T */

void pkn_MultMatrixTSubd ( int nrows_a, int rowlen_a, int pitch_a, const double *a,
                           int nrows_b, int pitch_b, const double *b,
                           int pitch_c, double *c )
{
  int    i, j, k, ia, jb, ic;
  double s;

  for ( i = ia = ic = 0;
        i < nrows_a;
        i++, ia += pitch_a, ic += pitch_c )
    for ( j = jb = 0;  j < nrows_b;  j++, jb += pitch_b ) {
      s = a[ia]*b[jb];
      for ( k = 1; k < rowlen_a; k++ )
        s += a[ia+k]*b[jb+k];
      c[ic+j] -= (double)s;
    }
} /*pkn_MultMatrixTSubd*/

