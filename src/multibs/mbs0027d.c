
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

boolean mbs_BCHornerDer2Pd ( int degreeu, int degreev, int spdimen,   
                             const double *ctlpoints,   
                             double u, double v,   
                             double *p, double *du, double *dv,
                             double *duu, double *duv, double *dvv )
{
  void    *sp;
  double  *aux;
  int     size;
  boolean result;

  size = (36+3*(degreev+1))*spdimen;
  sp = aux = pkv_GetScratchMemd ( size );
  if ( aux ) {
    result = _mbs_BCHornerDer2Pd ( degreeu, degreev, spdimen, ctlpoints,
                                   u, v, p, du, dv, duu, duv, dvv, aux );
    pkv_SetScratchMemTop ( sp );
    return result;
  }
  else
    return false;
} /*mbs_BCHornerDer2Pd*/

