
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

boolean mbs_BCHornerDerP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                              float u, float v,
                              point3f *p, vector3f *du, vector3f *dv )
{
  float   *workspace;
  boolean result;

  workspace = pkv_GetScratchMemf ( (6+2*(degreev+1))*3 );
  if ( workspace ) {
    result = _mbs_BCHornerDerP3Rf ( degreeu, degreev, ctlpoints, u, v,
                                    p, du, dv, workspace );
    pkv_SetScratchMemTop ( workspace );
    return result;
  }
  else
    return false;
} /*mbs_BCHornerDerP3Rf*/

