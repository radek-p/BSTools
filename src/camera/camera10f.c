
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "camera.h"

#include "msgpool.h"


void CameraSetDepthRangef ( CameraRecf *CPos, float zmin, float zmax )
{
  if ( zmax > zmin ) {
    CPos->zmin = zmin;
    CPos->zmax = zmax;
    CPos->ncplanes = 6;
  }
  else
    CPos->ncplanes = 4;
} /*CameraSetDepthRangef*/

