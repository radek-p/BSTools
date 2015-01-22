
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

boolean mbs_BCHornerDerPd ( int degreeu, int degreev, int spdimen,
                            const double *ctlpoints,
                            double u, double v,
                            double *p, double *du, double *dv )
{
  void    *sp;
  double  *aux;
  int     scr_size;
  boolean result;

  if ( degreeu == 0 ) {
    if ( degreev == 0 )
      return _mbs_BCHornerDerPd ( degreeu, degreev, spdimen, ctlpoints,
                                  u, v, p, du, dv, NULL );
    else {
      scr_size = 2*spdimen;
      goto call_it;
    }
  }
  else if ( degreev == 0 ) {
    scr_size = 2*spdimen;
    goto call_it;
  }
  else {
    scr_size = (6+2*(degreev+1))*spdimen;
call_it:
    sp = pkv_GetScratchMemTop ();
    aux = pkv_GetScratchMemd ( scr_size );
    if ( aux )
      result = _mbs_BCHornerDerPd ( degreeu, degreev, spdimen, ctlpoints,
                                    u, v, p, du, dv, aux );
    else
      result = false;
    pkv_SetScratchMemTop ( sp );
    return result;
  }
} /*mbs_BCHornerDerPd*/

