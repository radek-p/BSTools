
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_BCHornerC2Rd ( int degree, const point3d *ctlpoints, double t,
                           point2d *cpoint )
{
  point3d hcpoint;

  if ( mbs_multiBCHornerd ( degree, 1, 3, 0, (double*)ctlpoints, t,
                           (double*)&hcpoint ) ) {
    Point3to2d ( &hcpoint, cpoint );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerC2Rd*/

