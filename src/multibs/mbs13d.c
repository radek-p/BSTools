
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
/* computing values of B-spline functions at t */
/* fnz - index of the first nonzero function,  */
/* whose value is stored in bfv[0].            */

void mbs_deBoorBasisd ( int degree, int lastknot, const double *knots,  
                        double t, int *fnz, int *nnz, double *bfv )
{
  int i, j, k, l, r;
  double alpha, beta;
  double *auxkn, *auxbfv;
  int   auxsize;

  k = mbs_FindKnotIntervald ( -1, lastknot, knots, t, &r );

                     /* first, there are two special cases to deal with */
                     /* copying data, extending the knot sequence       */
                     /* and a recursive call to this procedure          */
  if ( k < degree ) {
    if ( k < 0 ) {
      *fnz = *nnz = 0;
      return;
    }
    auxkn = pkv_GetScratchMemd ( auxsize = 3*(degree+1) );
    auxbfv = &auxkn[2*(degree+1)];
    l = degree-k;
    memcpy ( &auxkn[l], knots, (2*(degree+1)-l)*sizeof(double) );
    for ( i = 0; i < l; i++ )
      auxkn[i] = knots[0];
    mbs_deBoorBasisd ( degree, 2*degree+1, auxkn, t, fnz, nnz, auxbfv );
    if ( *nnz > l ) {
      *nnz -= l;
      memcpy ( bfv, &auxbfv[l], (*nnz)*sizeof(double) );
    }
    else
      *nnz = 0;
    *fnz = 0;
    pkv_FreeScratchMemd ( auxsize );
  }
  else if ( k >= lastknot-degree ) {
    if ( r == degree+1 ) {
      *fnz = lastknot-degree-1;  *nnz = 1;
      bfv[0] = 1.0;
      return;
    }
    if ( k >= lastknot ) {
      *fnz = *nnz = 0;
      return;
    }
    auxkn = pkv_GetScratchMemd ( auxsize = 3*(degree+1) );
    auxbfv = &auxkn[2*(degree+1)];
    l = k+1+degree-lastknot;
    memcpy ( auxkn, &knots[k-degree], (2*(degree+1)-l)*sizeof(double) );
    for ( i = 2*(degree+1)-l; i < 2*(degree+1); i++ )
      auxkn[i] = knots[lastknot];
    mbs_deBoorBasisd ( degree, 2*degree+1, auxkn, t, fnz, nnz, auxbfv );
    *fnz = k-degree;
    if ( (*fnz)+(*nnz) > lastknot-degree )
      *nnz = lastknot-degree-(*fnz);
    if ( *nnz > 0 )
      memcpy ( bfv, auxbfv, (*nnz)*sizeof(double) );
    pkv_FreeScratchMemd ( auxsize );
  }
  else {
                     /* now the regular case, */
                     /* with knots[degree] <= t < knots[lastknot-degree] */
    if ( r > degree )
      r = degree;
    *fnz = l = k-degree;
    *nnz = degree+1-r;

    bfv[degree] = 1.0;
    for ( j = 1; j <= degree; j++ ) {
      beta = (knots[k+1]-t)/(knots[k+1]-knots[k-j+1]);
      bfv[k-j-l] = beta*bfv[k-j+1-l];
      for ( i = k-j+1; i < k; i++ ) {
        alpha = 1.0-beta;
        beta = (knots[i+j+1]-t)/(knots[i+j+1]-knots[i+1]);
        bfv[i-l] = alpha*bfv[i-l] + beta*bfv[i+1-l];
      }
      bfv[k-l] = (1.0-beta)*bfv[k-l];
    }
  }
} /*mbs_deBoorBasisd*/

