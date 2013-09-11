
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* /////////////////////////////////////////// */
/* constructing Bezier curves of interpolation, solving the Hermite */
/* interpolation problem with knots of interpolation 0 and 1, */
/* with arbitrary multiplicities. */

boolean mbs_multiInterp2knHermiteBezf ( int ncurves, int spdimen, int degree,
                                        int nlbc, int lbcpitch, const float *lbc,
                                        int nrbc, int rbcpitch, const float *rbc,
                                        int pitch, float *ctlpoints )
{
  int i, j;

  if ( nlbc+nrbc != degree+1 ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    return false;
  }

        /* get the interpolation conditions */
  pkv_Selectf ( ncurves, spdimen*nlbc, lbcpitch, pitch, lbc, ctlpoints );
          /* at the right end reverse */
  for ( i = 0; i < ncurves; i++ )
    pkv_Selectf ( nrbc, spdimen, spdimen, -spdimen, &rbc[i*rbcpitch],
                  &ctlpoints[i*pitch+degree*spdimen] );

        /* now the proper computation */
          /* at the left end */
  for ( j = nlbc-2; j >= 0; j-- )
    for ( i = 0; i <= nlbc-2-j; i++ )
      pkn_AddMatrixMf ( ncurves, spdimen, pitch, &ctlpoints[(j+i)*spdimen],
                   pitch, &ctlpoints[(j+i+1)*spdimen], 1.0/(float)(degree-j),
                   pitch, &ctlpoints[(j+i+1)*spdimen] );
          /* at the right end */
  for ( j = nrbc-2; j >= 0; j-- )
    for ( i = degree-j; i >= degree-nrbc+2; i-- )
      pkn_AddMatrixMf ( ncurves, spdimen, pitch, &ctlpoints[i*spdimen],
                   pitch, &ctlpoints[(i-1)*spdimen], -1.0/(float)(degree-j),
                   pitch, &ctlpoints[(i-1)*spdimen] );
  return true;
} /*mbs_multiInterp2knHermiteBezf*/

