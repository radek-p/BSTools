
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

void mbs_multiBCDegRedf ( int ncurves, int spdimen,
                          int inpitch, int indegree, const float *inctlpoints,
                          int deltadeg,
                          int outpitch, int *outdegree, float *outctlpoints )
{
  void  *sp;
  float *a, *b, e;
  int   i, j, n, k, c, d;

  if ( deltadeg > indegree )
    return;
  else if ( !deltadeg ) {
    pkv_Selectf ( ncurves, spdimen*(indegree+1), inpitch, outpitch,
                  inctlpoints, outctlpoints );
    *outdegree = indegree;
    return;
  }
  sp = pkv_GetScratchMemTop ();
  n = indegree-deltadeg;
  a = pkv_GetScratchMemf ( ((n+1)*(n+2))/2 );
  b = pkv_GetScratchMemf ( (n+1)*(indegree+1) );
  if ( a && b ) {
        /* setup the matrices of scalar products */
    c = pkv_Binom ( 2*n, n );
    e = 1.0/(float)((2*n+1)*c);
    for ( i = k = 0;  i <= n;  i++, k++ ) {
      d = c;
      for ( j = 0;  j < i;  j++, k++ ) {
        a[k] = (float)d*e;
        d = ((d*(i+j+1))/(j+1)*(n-j))/(2*n-i-j);
      }
      b[k] = (float)d*e;
      c = (c*(n-i))/(2*n-i);
    }
    c = pkv_Binom ( indegree+n, n );
    e = 1.0/(float)((n+indegree+1)*c);
    for ( i = k = 0;  i <= n;  i++, k++ ) {
      d = c;
      for ( j = 0;  j < indegree;  j++, k++ ) {
        a[k] = (float)d*e;
        d = ((d*(i+j+1))/(j+1)*(indegree-j))/(n-i+indegree-j);
      }
      b[k] = (float)d*e;
      c = c*(n-i)/(n+indegree-i);
    }
    if ( !pkn_CholeskyDecompf ( n+1, a ) )
      goto out;
        /* degree reduction for each curve */
    for ( i = 0; i < ncurves; i++ ) {
      pkn_MultMatrixf ( n+1, indegree+1, indegree+1, b,
                        spdimen, spdimen, &inctlpoints[i*inpitch],
                        spdimen, &outctlpoints[i*outpitch] );
      pkn_LowerTrMatrixSolvef ( n+1, a, spdimen, spdimen,
                                &outctlpoints[i*outpitch],
                                spdimen, &outctlpoints[i*outpitch] );
      pkn_UpperTrMatrixSolvef ( n+1, a, spdimen, spdimen,
                                &outctlpoints[i*outpitch],
                                spdimen, &outctlpoints[i*outpitch] );
    }
    *outdegree = n;
  }
out:
  pkv_SetScratchMemTop ( sp );
} /*mbs_multiBCDegRedf*/

