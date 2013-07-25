
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

#include "multibs.h"


boolean mbs_ApproxBSKnotsValidf ( int degree, int lastknot, const float *knots,
                                  int lastiknot, const float *iknots )
{
  int i, k;
                       /* check that the knot sequences are correct */
  for ( k = 0; k < lastknot; k++ )
    if ( knots[k] > knots[k+1] )         /* should be nondecreasing */
      return false;
  for ( k = 1; k < lastknot-degree-1; k++ )
    if ( knots[k] == knots[k+degree] )   /* multiplicities at most degree */
      return false;
  for ( i = 0; i < lastiknot; i++ )
    if ( iknots[k] >= iknots[k+1] )      /* this sequence must be strictly */
      return false;                      /* increasing */
  if ( iknots[0] < knots[degree] ||
       iknots[lastiknot] > knots[lastknot-degree] )
    return false;          /* knots of interpolation must be in the domain */

  return true;
} /*mbs_ApproxBSKnotsValidf*/

int mbs_ApproxBSBandmSizef ( int degree, const float *knots,
                             int lastiknot, const float *iknots )
{
  int i, k, r, msize;

  msize = 0;
  for ( i = 0, k = degree; i <= lastiknot; i++ ) {
                         /* find the proper interval */
    while ( iknots[i] >= knots[k+1] )
      k++;
                         /* find the multiplicity */
    r = 0;
    while ( iknots[i] == knots[k-r] && r < degree )
      r++;
                         /* add the number of nonzero coefficients */
    msize += degree-r+1;
  }

  return msize;
} /*mbs_ApproxBSBandmSizef*/

boolean mbs_ConstructApproxBSProfilef ( int degree, int lastknot,
                                        const float *knots,
                                        int lastiknot, const float *iknots,
                                        bandm_profile *prof )
{
  int i, j, k, r;

                       /* now the proper algorithm */
                         /* initialize the profile empty */
  for ( k = 0; k <= lastknot-degree; k++ ) {
    prof[k].firstnz = lastiknot+1;
    prof[k].ind = 0;
  }
  for ( i = 0, k = degree; i <= lastiknot; i++ ) {
                         /* find the proper interval */
    while ( iknots[i] >= knots[k+1] )
      k++;
                         /* find the multiplicity */
    r = 0;
    while ( iknots[i] == knots[k-r] && r < degree )
      r++;
                         /* update the matrix profile */
    for ( j = k-degree; j <= k-r; j++ ) {
      if ( prof[j].firstnz > i ) {
        prof[j].firstnz = i;
        prof[j+1].ind = 1;
      }
      else
        prof[j+1].ind = i-prof[j].firstnz+1;
    }
  }
                         /* now compute the indexes in the profile */
  for ( k = 0, i = 0; k < lastknot-degree; k++ )
    prof[k+1].ind += prof[k].ind;

  return true;
} /*mbs_ConstructApproxBSProfilef*/

boolean mbs_ConstructApproxBSMatrixf ( int degree, int lastknot,
                                       const float *knots,
                                       int lastiknot, const float *iknots,
                                       int *nrows, int *ncols,
                                       bandm_profile *prof,
                                       float *a )
{
  int   i, j, nr, scr_size;
  float *bfv;
  int   fnz, nnz;

  if ( mbs_ConstructApproxBSProfilef (degree, lastknot, knots,
                                      lastiknot, iknots, prof) ) {
    if ( (bfv = pkv_GetScratchMemf ( scr_size = degree+1 )) ) {
      *ncols = lastknot-degree;
      *nrows = nr = lastiknot+1;
      for ( i = 0; i < nr; i++ ) {
        mbs_deBoorBasisf ( degree, lastknot, knots, iknots[i], &fnz, &nnz, bfv );
        for ( j = fnz; j < fnz+nnz; j++ )
          a[prof[j].ind+i-prof[j].firstnz] = bfv[j-fnz];
      }
      pkv_FreeScratchMemf ( scr_size );
      return true;
    }
    else
      return false;
  }
  else
    return false;
} /*mbs_ConstructApproxBSMatrixf*/

boolean mbs_multiConstructApproxBSCf ( int degree, int lastknot,
                                       const float *knots,
                                       int lastpknot, const float *pknots,
                                       int ncurves, int spdimen,
                                       int ppitch, const float *ppoints,
                                       int bcpitch, float *ctlpoints )
{
  int           ncols, nrows, msize, qsize, rsize;
  float         *amat, *qmat, *rmat, *pp;
  bandm_profile *aprof, *qprof, *rprof;
  int           i;
  boolean       result;
  void          *scratchsp;

  result = false;
  scratchsp = pkv_GetScratchMemTop ();
  if ( mbs_ApproxBSKnotsValidf(degree,lastknot,knots,lastpknot,pknots) ) {
    ncols = lastknot-degree;
    nrows = lastpknot+1;
    aprof = (bandm_profile*)pkv_GetScratchMem((ncols+1)*sizeof(bandm_profile));
    qprof = (bandm_profile*)pkv_GetScratchMem((ncols+1)*sizeof(bandm_profile));
    rprof = (bandm_profile*)pkv_GetScratchMem((ncols+1)*sizeof(bandm_profile));
    msize = mbs_ApproxBSBandmSizef ( degree, knots, lastpknot, pknots );
    amat = pkv_GetScratchMemf ( msize );
    if ( !aprof || !qprof || !rprof || !amat )
      goto cleanup;
    if ( !mbs_ConstructApproxBSMatrixf(degree, lastknot, knots,
               lastpknot, pknots, &nrows, &ncols, aprof, amat) )
      goto cleanup;
    pkn_BandmFindQRMSizes ( ncols, aprof, &qsize, &rsize );
    qmat = pkv_GetScratchMemf ( qsize );
    rmat = pkv_GetScratchMemf ( rsize );
    pp = pkv_GetScratchMemf ( nrows*spdimen );
    if ( !qmat || !rmat || !pp )
      goto cleanup;
    pkn_BandmQRDecomposeMatrixf ( nrows, ncols, aprof, amat,
                                  qprof, qmat, rprof, rmat );
    for ( i = 0;
          i < ncurves;
          i++, ppoints += ppitch, ctlpoints += bcpitch ) {
      memcpy ( pp, ppoints, nrows*spdimen*sizeof(float) );
      pkn_multiBandmReflectVectorf ( ncols, qprof, qmat, spdimen, pp );
      pkn_multiBandmMultInvUTMVectorf ( ncols, rprof, rmat, spdimen, pp, pp );
      memcpy ( ctlpoints, pp, ncols*spdimen*sizeof(float) );
    }

    result = true;
  }

cleanup:
  pkv_SetScratchMemTop ( scratchsp );
  return result;
} /*mbs_multiConstructApproxBSCf*/

