
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

boolean mbs_BCHornerDer2P3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                               double u, double v,
                               point3d *p, vector3d *du, vector3d *dv,
                               vector3d *duu, vector3d *duv, vector3d *dvv )
{
  void    *sp;
  double  *aux;
  int     size;
  boolean result;

  size = (36+3*(degreev+1))*4;
  sp = aux = pkv_GetScratchMemd ( size );
  if ( aux ) {
    result = _mbs_BCHornerDer2P3Rd ( degreeu, degreev, ctlpoints, u, v,
                                     p, du, dv, duu, duv, dvv, aux );
    pkv_SetScratchMemTop ( sp );
    return result;
  }
  else
    return false;
} /*mbs_BCHornerDer2P3Rd*/

