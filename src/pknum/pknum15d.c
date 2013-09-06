
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


/* solve a single dual linear least squares problem with a band matrix */

void pkn_multiBandmSolveDLSQd ( int nrows, int ncols,
                                const bandm_profile *atprof, const double *at,
                                int nrsides, int spdimen,
                                int bpitch, const double *b,
                                int x0pitch, const double *x0,
                                int xpitch, double *x )
{
  void          *sp;
  int           qsize, rsize, i;
  bandm_profile *qprof, *rprof;
  double        *qa, *ra;

  sp = pkv_GetScratchMemTop ();

  pkn_BandmFindQRMSizes ( ncols, atprof, &qsize, &rsize );
  qprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  rprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  qa = pkv_GetScratchMemd ( qsize );
  ra = pkv_GetScratchMemd ( rsize );
  if ( !qprof || !rprof || !qa || !ra ) {
    PKV_SIGNALERROR ( LIB_PKNUM, 2, ERRMSG_2 );
    exit ( 1 );
  }

  pkn_BandmQRDecomposeMatrixd ( nrows, ncols, atprof, at, qprof, qa, rprof, ra );
  for ( i = 0; i < nrsides; i++ ) {
    if ( x0 ) {
      memcpy ( &x[i*xpitch], &x0[i*x0pitch], spdimen*nrows*sizeof(double) );
      pkn_multiBandmReflectVectord ( ncols, qprof, qa, spdimen, x );
    }
    else
      memset ( &x[i*xpitch], 0, spdimen*nrows*sizeof(double) );
    if ( b )
      pkn_multiBandmMultInvTrUTMVectord ( ncols, rprof, ra, spdimen,
                                          &b[i*bpitch], &x[i*xpitch] );
    else
      memset ( &x[i*xpitch], 0, spdimen*ncols*sizeof(double) );
    pkn_multiBandmInvReflectVectord ( ncols, qprof, qa, spdimen, &x[i*xpitch] );
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_multiBandmSolveDLSQd*/

