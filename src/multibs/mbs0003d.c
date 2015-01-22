
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

boolean mbs_BCHornerC3Rd ( int degree, const point4d *ctlpoints, double t,
                           point3d *cpoint )
{
  point4d hcpoint;

  if ( mbs_multiBCHornerd ( degree, 1, 4, 0, (double*)ctlpoints, t,
                            (double*)&hcpoint ) ) {
    Point4to3d ( &hcpoint, cpoint );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerC3Rd*/

