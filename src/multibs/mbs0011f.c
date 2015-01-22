
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

boolean mbs_BCHornerDerC3Rf ( int degree, const point4f *ctlpoints, float t,
                              point3f *p, vector3f *d )
{
  point4f hp, hd;

  if ( !mbs_BCHornerDerC4f ( degree, ctlpoints, t, &hp, &hd ) )
    return false;
  Point4to3f ( &hp, p );
  memcpy ( d, &hd, sizeof(vector3f) );
  AddVector3Mf ( d, p, -hd.w, d );
  MultVector3f ( 1.0/hp.w, d, d );
  return true;
} /*mbs_BCHornerDerC3Rf*/

