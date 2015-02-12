
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

boolean _mbs_BezC2CoonsFindCornersd ( int spdimen,
                                      int degc00, const double *c00,
                                      int degc01, const double *c01,
                                      int degc02, const double *c02,
                                      int degc10, const double *c10,
                                      int degc11, const double *c11,
                                      int degc12, const double *c12,
                                      double *pcorners,
                                      double *workspace )
{
  if ( !_mbs_multiBCHornerDer2d ( degc00, 1, spdimen, 0, c00, 0.0,
        &pcorners[0], &pcorners[spdimen*12], &pcorners[spdimen*24], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc00, 1, spdimen, 0, c00, 1.0,
        &pcorners[spdimen*6], &pcorners[spdimen*18], &pcorners[spdimen*30], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc10, 1, spdimen, 0, c10, 0.0,
        &pcorners[spdimen*1], &pcorners[spdimen*13], &pcorners[spdimen*25], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc10, 1, spdimen, 0, c10, 1.0,
        &pcorners[spdimen*7], &pcorners[spdimen*19], &pcorners[spdimen*31], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc01, 1, spdimen, 0, c01, 0.0,
        &pcorners[spdimen*2], &pcorners[spdimen*14], &pcorners[spdimen*26], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc01, 1, spdimen, 0, c01, 1.0,
        &pcorners[spdimen*8], &pcorners[spdimen*20], &pcorners[spdimen*32], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc11, 1, spdimen, 0, c11, 0.0,
        &pcorners[spdimen*3], &pcorners[spdimen*15], &pcorners[spdimen*27], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc11, 1, spdimen, 0, c11, 1.0,
        &pcorners[spdimen*9], &pcorners[spdimen*21], &pcorners[spdimen*33], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc02, 1, spdimen, 0, c02, 0.0,
        &pcorners[spdimen*4], &pcorners[spdimen*16], &pcorners[spdimen*28], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc02, 1, spdimen, 0, c02, 1.0,
        &pcorners[spdimen*10], &pcorners[spdimen*22], &pcorners[spdimen*34], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc12, 1, spdimen, 0, c12, 0.0,
        &pcorners[spdimen*5], &pcorners[spdimen*17], &pcorners[spdimen*29], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2d ( degc12, 1, spdimen, 0, c12, 1.0,
        &pcorners[spdimen*11], &pcorners[spdimen*23], &pcorners[spdimen*35], workspace ) )
    return false;
  return true;
} /*_mbs_BezC2CoonsFindCornersd*/

