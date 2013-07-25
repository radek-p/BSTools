
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

#include "msgpool.h"

/* /////////////////////////////////////////// */
/* constructing spline curves of interpolation, solving the Hermite */
/* interpolation problem with knots of interpolation at the end points */
/* of the domain, with arbitrary multiplicities. */

void mbs_multiInterp2knHermiteBSf ( int ncurves, int spdimen, int degree,
                                    int lastknot, const float *knots,
                                    int nlbc, int lbcpitch, const float *lbc,
                                    int nrbc, int rbcpitch, const float *rbc,
                                    int pitch, float *ctlpoints )
{
  int i, j;

  if ( nlbc+nrbc != lastknot-degree ) {
    pkv_SignalError ( LIB_MULTIBS, 33, ERRMSG_0 );
    exit ( 1 );
  }

        /* get the interpolation conditions */
  pkv_Selectf ( ncurves, spdimen*nlbc, lbcpitch, pitch, lbc, ctlpoints );
          /* at the right end reverse */
  for ( i = 0; i < ncurves; i++ )
    pkv_Selectf ( nrbc, spdimen, spdimen, -spdimen, &rbc[i*rbcpitch],
                  &ctlpoints[i*pitch+(lastknot-degree-1)*spdimen] );

        /* now the proper computation */
          /* at the left end */
  for ( j = nlbc-2; j >= 0; j-- )
    for ( i = 0; i <= nlbc-2-j; i++ )
      pkn_AddMatrixMf ( ncurves, spdimen, pitch, &ctlpoints[(j+i)*spdimen],
                        pitch, &ctlpoints[(j+i+1)*spdimen],
                        (knots[i+degree+1]-knots[i+j+1])/(float)(degree-j),
                        pitch, &ctlpoints[(j+i+1)*spdimen] );
          /* at the right end */
  for ( j = nrbc-2; j >= 0; j-- )
    for ( i = lastknot-degree-j-1; i >= lastknot-degree-nrbc+1; i-- )
      pkn_AddMatrixMf ( ncurves, spdimen, pitch, &ctlpoints[i*spdimen],
                        pitch, &ctlpoints[(i-1)*spdimen],
                        -(knots[i+degree]-knots[i+j])/(float)(degree-j),
                        pitch, &ctlpoints[(i-1)*spdimen] );
} /*mbs_multiInterp2knHermiteBSf*/

