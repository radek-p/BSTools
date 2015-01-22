
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

boolean _mbs_BCHornerPf ( int degreeu, int degreev, int spdimen,
                          const float *ctlpoints,
                          float u, float v, float *ppoint,
                          float *workspace )
{
  if ( !mbs_multiBCHornerf ( degreeu, 1, spdimen*(degreev+1), 0, ctlpoints,
                             u, workspace ) )
    return false;
  return mbs_multiBCHornerf ( degreev, 1, spdimen, 0, workspace, v, ppoint );
} /*_mbs_BCHornerPf*/

