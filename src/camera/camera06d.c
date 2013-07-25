
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "camera.h"

#include "msgpool.h"


void CameraSetFd ( CameraRecd *CPos, double f )
{
  CPos->vd.persp.f = f;
  CameraSetMappingd ( CPos );
} /*CameraSetFd*/

void CameraZoomd ( CameraRecd *CPos, double fchange )
{
  CPos->vd.persp.f = fchange * CPos->vd.persp.f;
  CameraSetMappingd ( CPos );
} /*CameraZoomd*/

