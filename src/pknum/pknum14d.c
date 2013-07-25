
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

#include "msgpool.h"


/* solve a regular linear least squares problem with a band matrix */

void pkn_multiBandmSolveRLSQd ( int nrows, int ncols,
                                const bandm_profile *aprof, const double *a,
                                int nrsides, int spdimen,
                                int bpitch, const double *b,
                                int xpitch, double *x )
{
  void          *sp;
  int           qsize, rsize, i;
  bandm_profile *qprof, *rprof;
  double        *qa, *ra, *y;

  sp = pkv_GetScratchMemTop ();

  pkn_BandmFindQRMSizes ( ncols, aprof, &qsize, &rsize );
  qprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  rprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  qa = pkv_GetScratchMemd ( qsize );
  ra = pkv_GetScratchMemd ( rsize );
  y = pkv_GetScratchMemd ( nrows*spdimen );
  if ( !qprof || !rprof || !qa || !ra || !y ) {
    pkv_SignalError ( LIB_PKNUM, 0, ERRMSG_0 );
    exit ( 1 );
  }

  pkn_BandmQRDecomposeMatrixd ( nrows, ncols, aprof, a, qprof, qa, rprof, ra );
  for ( i = 0; i < nrsides; i++ ) {
    memcpy ( y, &b[i*bpitch], nrows*spdimen*sizeof(double) );
    pkn_multiBandmReflectVectord ( ncols, qprof, qa, spdimen, y );
    pkn_multiBandmMultInvUTMVectord ( ncols, rprof, ra, spdimen, y,
                                      &x[i*xpitch] );
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_multiBandmSolveRLSQd*/

