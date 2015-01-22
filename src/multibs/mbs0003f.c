
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

boolean mbs_BCHornerC3Rf ( int degree, const point4f *ctlpoints, float t,
                           point3f *cpoint )
{
  point4f hcpoint;

  if ( mbs_multiBCHornerf ( degree, 1, 4, 0, (float*)ctlpoints, t,
                            (float*)&hcpoint ) ) {
    Point4to3f ( &hcpoint, cpoint );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerC3Rf*/

