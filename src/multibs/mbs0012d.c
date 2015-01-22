
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
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

boolean _mbs_BCHornerDerPd ( int degreeu, int degreev, int spdimen,
                             const double *ctlpoints,
                             double u, double v,
                             double *p, double *du, double *dv,
                             double *workspace )
{
  double *w1, *w2, *w4, *w5, *w6;

  if ( degreeu == 0 ) {
    if ( degreev == 0 ) {
      memcpy ( p, ctlpoints, spdimen*sizeof(double) );
      memset ( dv, 0, spdimen*sizeof(double) );
    }
    else {
      if ( !_mbs_multiBCHornerDerd ( degreev, 1, spdimen, 0, ctlpoints, v, p, dv,
                                     workspace ) )
        return false;
    }
    memset ( du, 0, spdimen*sizeof(double) );
  }
  else if ( degreev == 0 ) {
    if ( !_mbs_multiBCHornerDerd ( degreeu, 1, spdimen, 0, ctlpoints, u, p, du,
                                   workspace ) )
      return false;
    memset ( dv, 0, spdimen*sizeof(double) );
  }
  else {
    w1 = &workspace[spdimen];
    w2 = &workspace[2*spdimen];
    w4 = &workspace[4*spdimen];
    w5 = &workspace[5*spdimen];
    w6 = &workspace[6*spdimen];
    if ( !mbs_multiBCHornerd ( degreeu-1, 2, spdimen*(degreev+1),
                               spdimen*(degreev+1), ctlpoints, u, w6 ) )
      return false;
    if ( !mbs_multiBCHornerd ( degreev-1, 2, spdimen, spdimen,
                               w6, v, workspace ) )
      return false;
    if ( !mbs_multiBCHornerd ( degreev-1, 2, spdimen, spdimen,
                               &w6[spdimen*(degreev+1)], v, w2 ) )
      return false;
    pkn_SubtractMatrixd ( 2, spdimen, spdimen, w2, spdimen, workspace,
                          spdimen, w4 );
    pkn_MatrixLinCombd ( 1, spdimen, 0, w4, (1.0-v)*(double)degreeu,
                         0, w5, v*(double)degreeu, 0, du );
    pkn_MatrixLinCombd ( 2, spdimen, spdimen, workspace, 1.0-u,
                         spdimen, w2, u, spdimen, workspace );
    pkn_MatrixMDifferenced ( 1, spdimen, 0, w1, 0, workspace,
                             (double)degreev, 0, dv );
    pkn_MatrixLinCombd ( 1, spdimen, 0, workspace, 1.0-v, 0, w1, v, 0, p );
  }
  return true;
} /*_mbs_BCHornerDerPd*/

