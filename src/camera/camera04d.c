
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

void CameraMoveToGd ( CameraRecd *CPos, point3d *pos )
{
  /* This procedure sets the camera position to pos, given in global */
  /* coordinates. */
  CPos->position = *pos;
  CameraSetMappingd ( CPos );
  _CameraUpdateRotCentred ( CPos );
} /*CameraMoveToGd*/

void CameraTurnGd ( CameraRecd *CPos, double _psi, double _theta, double _phi )
{
  /* This procedure sets the camera Euler angles, turning it around */
  /* the its current rotations centre */
  /* Find the camera position */
  CameraRotGd ( CPos, -CPos->phi, -CPos->theta, -CPos->psi );
  /* The assignment below anihilates rounding errors of angles */
  CPos->psi = CPos->theta = CPos->phi = 0.0;
  CameraRotGd ( CPos, _psi, _theta, _phi );
} /*CameraTurnGd*/

void CameraMoveGd ( CameraRecd *CPos, vector3d *v )
{
  AddVector3d ( &CPos->position, v, &CPos->position );
  CameraSetMappingd ( CPos );
  _CameraUpdateRotCentred ( CPos );
}  /*CameraMoveGd*/

