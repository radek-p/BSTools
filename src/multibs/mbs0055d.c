
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

boolean _mbs_BezC1CoonsFindCornersd ( int spdimen,
                                     int degc00, const double *c00,
                                     int degc01, const double *c01,
                                     int degc10, const double *c10,
                                     int degc11, const double *c11,
                                     double *pcorners,
                                     double *workspace )
{
  if ( !_mbs_multiBCHornerDerd ( degc00, 1, spdimen, 0, c00, 0.0,
                  &pcorners[0], &pcorners[spdimen*8], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerd ( degc00, 1, spdimen, 0, c00, 1.0,
                  &pcorners[spdimen*4], &pcorners[spdimen*12], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerd ( degc10, 1, spdimen, 0, c10, 0.0,
                  &pcorners[spdimen*1], &pcorners[spdimen*9], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerd ( degc10, 1, spdimen, 0, c10, 1.0,
                  &pcorners[spdimen*5], &pcorners[spdimen*13], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerd ( degc01, 1, spdimen, 0, c01, 0.0,
                  &pcorners[spdimen*2], &pcorners[spdimen*10], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerd ( degc01, 1, spdimen, 0, c01, 1.0,
                  &pcorners[spdimen*6], &pcorners[spdimen*14], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerd ( degc11, 1, spdimen, 0, c11, 0.0,
                  &pcorners[spdimen*3], &pcorners[spdimen*11], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerd ( degc11, 1, spdimen, 0, c11, 1.0,
                  &pcorners[spdimen*7], &pcorners[spdimen*15], workspace ) )
    return false;
  return true;
} /*_mbs_BezC1CoonsFindCornersd*/

