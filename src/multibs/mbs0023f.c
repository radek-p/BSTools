
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

boolean mbs_BCHornerDer2C2Rf ( int degree, const point3f *ctlpoints, float t,
                               point2f *p, vector2f *d1, vector2f *d2 )
{
  void    *sp;
  float   *aux;
  boolean result;

  sp = aux = pkv_GetScratchMemf ( 3*3 );
  if ( aux ) {
    result = _mbs_BCHornerDer2C2Rf ( degree, ctlpoints, t, p, d1, d2, aux );
    pkv_SetScratchMemTop ( sp );
    return result;
  }
  else
    return false;
} /*mbs_BCHornerDer2C2Rf*/

