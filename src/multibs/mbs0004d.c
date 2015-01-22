
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

boolean _mbs_BCHornerPd ( int degreeu, int degreev, int spdimen,
                          const double *ctlpoints,
                          double u, double v, double *ppoint,
                          double *workspace )
{
  if ( !mbs_multiBCHornerd ( degreeu, 1, spdimen*(degreev+1), 0, ctlpoints,
                             u, workspace ) )
    return false;
  return mbs_multiBCHornerd ( degreev, 1, spdimen, 0, workspace, v, ppoint );
} /*_mbs_BCHornerPf*/

