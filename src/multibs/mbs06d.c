
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* Knot insertion                              */

int mbs_multiKnotInsd ( int degree, int *lastknot,
                        double *knots,
                        int ncurves, int spdimen, int inpitch, int outpitch,
                        double *ctlpoints, double t )
{
  int i, j, k, l, r, mp, nkn, is, ll;
  long double alpha, beta;

  if ( outpitch != inpitch )
    pkv_Rearranged ( ncurves, (*lastknot-degree)*spdimen,
                     inpitch, outpitch, ctlpoints );

  /* find the proper interval and current multiplicity the knot t */
  nkn = *lastknot;
  k = nkn-1;
  while ( t < knots[k] )
    k--;
  r = 0;  i = k;
  while ( i >= 1 && t == knots[i] ) {
    i--;
    r++;
  }

  /* actual knot insertion */
  mp = nkn-degree-(k-r);    /* number of points to move for each curve */
  if ( mp > 0 )
    pkv_Moved ( ncurves, mp*spdimen, outpitch, spdimen,
                &ctlpoints[(k-r)*spdimen] );
  for ( i = k-r, is = i*spdimen;
        i >= k-degree+1;
        i--, is -= spdimen ) {
    alpha = (t-knots[i])/(knots[i+degree]-knots[i]);
    beta = 1.0-alpha;
    for ( l = 0, ll = 0; l < ncurves; l++, ll += outpitch )
      for ( j = 0; j < spdimen; j++ )
        ctlpoints[ll+is+j] =
          (double)(beta*ctlpoints[ll+is-spdimen+j] + alpha*ctlpoints[ll+is+j]);
  }
  if ( nkn > k+1 )
    memmove ( &knots[k+2], &knots[k+1], (nkn-k)*sizeof(double) );
  knots[k+1] = t;
  (*lastknot)++;
  return k;
} /*mbs_multiKnotInsd*/


