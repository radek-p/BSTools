
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

#include "cprivate.h"
#include "msgpool.h"

void CameraMoveToGf ( CameraRecf *CPos, point3f *pos )
{
  /* This procedure sets the camera position to pos, given in global */
  /* coordinates. */
  CPos->position = *pos;
  CameraSetMappingf ( CPos );
  _CameraUpdateRotCentref ( CPos );
} /*CameraMoveToGf*/

void CameraTurnGf ( CameraRecf *CPos, float _psi, float _theta, float _phi )
{
  /* This procedure sets the camera Euler angles, turning it around */
  /* the its current rotations centre */
  /* Find the camera position */
  CameraRotGf ( CPos, -CPos->phi, -CPos->theta, -CPos->psi );
  /* The assignment below anihilates rounding errors of angles */
  CPos->psi = CPos->theta = CPos->phi = 0.0;
  CameraRotGf ( CPos, _psi, _theta, _phi );
} /*CameraTurnGf*/

void CameraMoveGf ( CameraRecf *CPos, vector3f *v )
{
  AddVector3f ( &CPos->position, v, &CPos->position );
  CameraSetMappingf ( CPos );
  _CameraUpdateRotCentref ( CPos );
}  /*CameraMoveGf*/

