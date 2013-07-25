
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"


/* ////////////////////////////////////////// */
/* Oslo algorithm of multiple knot insertion. */
/* The procedures here serve for the purpose  */
/* of building a B-spline basis change, from  */
/* the one with an initial knot sequence to   */
/* one with additional knots.                 */

boolean mbs_OsloKnotsCorrectf ( int lastuknot, const float *uknots,
                                int lastvknot, const float *vknots )
{
  int i, k;

  /* two conditions to be checked: both knot sequences must be nondecreasing */
  for ( i = 0; i < lastuknot; i++ )
    if ( uknots[i] > uknots[i+1] )
      return false;
  for ( k = 0; k < lastvknot; k++ )
    if ( vknots[k] > vknots[k+1] )
      return false;

  /* and the first sequence is a subsequence of the second one.              */
  for ( i = 0, k = 0; i <= lastuknot; i++ ) {
    while ( vknots[k] < uknots[i] ) {
      k++;
      if ( k > lastvknot )
        return false;
    }
    if ( vknots[k] > uknots[i] )
      return false;
    k++;
  }
  return true;
} /*mbs_OsloKnotsCorrectf*/

int mbs_BuildOsloMatrixProfilef ( int degree,
                                  int lastuknot, const float *uknots,
                                  int lastvknot, const float *vknots,
                                  bandm_profile *prof )
{
  int i, r, nc, ra, rb, ai, bi, asize;

  asize = 0;
  nc = lastuknot-degree;
  for ( i = 0; i < nc; i++ ) {
    ra = mbs_KnotMultiplicityf ( degree+1, &uknots[i], uknots[i] );
    ai = mbs_FindKnotIntervalf ( -1, lastvknot, vknots, uknots[i], NULL );
    ai += 1-ra;
    rb = mbs_KnotMultiplicityf ( degree+1, &uknots[i], uknots[i+degree+1] );
    bi = mbs_FindKnotIntervalf ( -1, lastvknot, vknots, uknots[i+degree+1], &r );
    if ( r > rb )
      { bi -= r-rb;  r = rb; }
    bi -= degree;
    prof[i].firstnz = ai;
    prof[i].ind = asize;
    asize += bi-ai;
  }
  prof[nc].firstnz = 0;
  prof[nc].ind = asize;
  return asize;
} /*mbs_BuildOsloMatrixProfilef*/

static void _mbs_FindOsloCoefficientsf ( int degree, const float *uknots,
                                         const float *vknots,
                                         int *nnz, double *oc )
{
  int    i, j, l;
  float  tk;
  double alpha, beta;

  l = mbs_FindKnotIntervalf ( -1, degree+1, uknots, vknots[0], NULL );
/*  memset ( oc, 0, (degree+1)*sizeof(double) ); */
  oc[l] = 1.0;
  for ( j = 1; j <= degree; j++ ) {
    tk = vknots[j];
    if ( j <= l ) {
      beta = (uknots[l+1]-tk)/(uknots[l+1]-uknots[l-j+1]);
      oc[l-j] = beta*oc[l-j+1];
    }
    else
      beta = 0.0;
    for ( i = l-j+1; i < l; i++ ) {
      if ( i > -1 ) {
        alpha = 1.0-beta;
        beta = (uknots[i+j+1]-tk)/(uknots[i+j+1]-uknots[i+1]);
        if ( i >= 0 )
          oc[i] = alpha*oc[i] + beta*oc[i+1];
      }
      else
        beta = 0.0;
    }
    oc[l] *= 1.0-beta;
  }
  *nnz = l+1;
} /*_mbs_FindOsloCoefficientsf*/

void mbs_BuildOsloMatrixf ( int degree, int lastuknot, const float *uknots,
                            const float *vknots,
                            const bandm_profile *prof, float *a )
{
  int    i, j, k, l, nnz;
  int    nc, size_oc;
  double *oc;

/* It is assumed that the prof array already contains the data computed  */
/* by the mbs_BuildOsloMatrixProfilef procedure for the knot sequences   */
/* specified by the parameters. Therefore here we just follow this data. */

  oc = pkv_GetScratchMemd ( size_oc = degree+1 );
  nc = lastuknot-degree;
  for ( i = j = 0; i < nc; i++ ) {
    for ( ; j < prof[i].firstnz+prof[i+1].ind-prof[i].ind; j++ ) {
      _mbs_FindOsloCoefficientsf ( degree, &uknots[i], &vknots[j], &nnz, oc );
      if ( i+nnz > nc )
        nnz = nc-i;
      for ( k = 0; k < nnz; k++ ) {
        if ( j >= prof[i+k].firstnz &&
             j < prof[i+k].firstnz + prof[i+k+1].ind - prof[i+k].ind ) {
          l = prof[i+k].ind + j - prof[i+k].firstnz;
          a[l] = (float)oc[k];
        }
      }
    }
  }
  pkv_FreeScratchMemd ( size_oc );
} /*mbs_BuildOsloMatrixf*/

