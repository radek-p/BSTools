
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"


/* solve a regular linear least squares problem with a band matrix */

boolean pkn_multiBandmSolveRLSQf ( int nrows, int ncols,
                                   const bandm_profile *aprof, const float *a,
                                   int nrsides, int spdimen,
                                   int bpitch, const float *b,
                                   int xpitch, float *x )
{
  void          *sp;
  int           qsize, rsize, i;
  bandm_profile *qprof, *rprof;
  float         *qa, *ra, *y;

  sp = pkv_GetScratchMemTop ();

  pkn_BandmFindQRMSizes ( ncols, aprof, &qsize, &rsize );
  qprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  rprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  qa = pkv_GetScratchMemf ( qsize );
  ra = pkv_GetScratchMemf ( rsize );
  y = pkv_GetScratchMemf ( nrows*spdimen );
  if ( !qprof || !rprof || !qa || !ra || !y ) {
    PKV_SIGNALERROR ( LIB_PKNUM, 2, ERRMSG_2 );
    goto failure;
  }

  pkn_BandmQRDecomposeMatrixf ( nrows, ncols, aprof, a, qprof, qa, rprof, ra );
  for ( i = 0; i < nrsides; i++ ) {
    memcpy ( y, &b[i*bpitch], nrows*spdimen*sizeof(float) );
    pkn_multiBandmReflectVectorf ( ncols, qprof, qa, spdimen, y );
    pkn_multiBandmMultInvUTMVectorf ( ncols, rprof, ra, spdimen, y,
                                      &x[i*xpitch] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_multiBandmSolveRLSQf*/

