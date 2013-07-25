
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "egholed.h"

#include "eghprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
/* The domain boundary curves are cubic, therefore a 6th order quadrature    */
/* is equal (up to rounding errors) to the actual integral we need. Here the */
/* 6th order Gauss-Legendre quadrature (with 3 knots) is used.               */

double gh_HoleDomainAread ( GHoleDomaind *domain, boolean symmetric )
{
  point2d          *surrpc, p;
  vector2d         d;
  double           qknots[3], qcoeff[3], f, area;
  int              hole_k, i, j;

  hole_k = domain->hole_k;
  pkn_QuadGaussLegendre6d ( 0.0, 1.0, 3, qknots, qcoeff );
  area = 0.0;
  if ( symmetric ) {
    surrpc = _gh_GetDomSurrndPatchd ( domain, 0, 1 );
    for ( j = 0; j < 3; j++ ) {
      mbs_BCHornerDerC2d ( 3, surrpc, qknots[j], &p, &d );
      f = p.x*d.y-p.y*d.x;
      area -= f*qcoeff[j];
    }
    area *= (double)hole_k;
  }
  else {
    for ( i = 0; i < hole_k; i++ ) {
      surrpc = _gh_GetDomSurrndPatchd ( domain, i, 1 );
      for ( j = 0; j < 3; j++ ) {
        mbs_BCHornerDerC2d ( 3, surrpc, qknots[j], &p, &d );
        f = p.x*d.y-p.y*d.x;
        area -= f*qcoeff[j];
      }
      surrpc = _gh_GetDomSurrndPatchd ( domain, i, 2 );
      for ( j = 0; j < 3; j++ ) {
        mbs_BCHornerDerC2d ( 3, surrpc, qknots[j], &p, &d );
        f = p.x*d.y-p.y*d.x;
        area += f*qcoeff[j];
      }
    }
    area *= 0.5;
  }
  return area;
} /*gh_HoleDomainAread*/

