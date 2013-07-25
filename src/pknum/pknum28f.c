
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pkvaria.h"
#include "pknum.h"

/* //////////////////////////////////////////////// */
/* Fast Fourier Transform, for n being a power of 2 */

int pkn_FFTPermutef ( int n, int rowlen, int pitch, complexf *a )
{
  int i, j, k, l, m;

  if ( n < 1 )
    return -1;
  l = 0;  m = n;
  while ( !(m & 0x01) )
    l++,  m /= 2;
  if ( m != 1 )
    return -1;
  for ( j = m = n / 2, i = 1;  i < n-1;  i++ ) {
    if ( i < j )
      pkv_Exchange ( &a[i*pitch], &a[j*pitch], rowlen*sizeof(complexf) );
    k = m;
    while ( k <= j ) { j = j-k;  k /= 2; }
      j += k;
  }
  return l;
} /*pkn_FFTPermutef*/

boolean pkn_FFTf ( int n, int rowlen, int pitch, complexf *a )
{
  complexf t, u, w;
  int      i, j, k, l, m, p, q, ipq, ippq;

  if ( n == 1 )
    return true;
  if ( ( l = pkn_FFTPermutef ( n, rowlen, pitch, a ) ) < 0 )
    return false;
  for ( m = k = 1; k <= l; k++ ) {
    p = m;  m += m;
    u.x = 1.0;  u.y = 0.0;
    w.x = (float)cos ( PI/(float)p );  w.y = (float)(-sin ( PI/(float)p ));
    for ( j = 0; j < p; j++ ) {
      for ( i = j; i < n; i += m ) {
        ipq = i*pitch;  ippq = (i+p)*pitch;
        for ( q = 0;  q < rowlen;  q++, ipq++, ippq++ ) {
          t.x = a[ippq].x*u.x-a[ippq].y*u.y;
          t.y = a[ippq].x*u.y+a[ippq].y*u.x;
          a[ippq].x = a[ipq].x-t.x;  a[ippq].y = a[ipq].y-t.y;
          a[ipq].x += t.x;           a[ipq].y += t.y;
        }
      }
      t.x = u.x*w.x-u.y*w.y;  u.y = u.x*w.y+u.y*w.x;  u.x = t.x;
    }
  }
  return true;
} /*pkn_FFTf*/

boolean pkn_InvFFTf ( int n, int rowlen, int pitch, complexf *a )
{
  complexf t, u, w;
  int      i, j, k, l, m, p, q, ipq, ippq;

  if ( n == 1 )
    return true;
  if ( ( l = pkn_FFTPermutef ( n, rowlen, pitch, a ) ) < 0 )
    return false;
  for ( m = k = 1; k <= l; k++ ) {
    p = m;  m += m;
    u.x = 1.0;  u.y = 0.0;
    w.x = (float)cos ( PI/(float)p );  w.y = (float)sin ( PI/(float)p );
    for ( j = 0; j < p; j++ ) {
      for ( i = j; i < n; i += m ) {
        ipq = i*pitch;  ippq = (i+p)*pitch;
        for ( q = 0;  q < rowlen;  q++, ipq++, ippq++ ) {
          t.x = a[ippq].x*u.x-a[ippq].y*u.y;
          t.y = a[ippq].x*u.y+a[ippq].y*u.x;
          a[ippq].x = a[ipq].x-t.x;  a[ippq].y = a[ipq].y-t.y;
          a[ipq].x += t.x;           a[ipq].y += t.y;
        }
      }
      t.x = u.x*w.x-u.y*w.y;  u.y = u.x*w.y+u.y*w.x;  u.x = t.x;
    }
  }
  pkn_MultMatrixNumf ( n, 2*rowlen, 2*pitch, &a[0].x, 1.0/(float)n,
                       2*pitch, &a[0].x );
  return true;
} /*pkn_InvFFTf*/

