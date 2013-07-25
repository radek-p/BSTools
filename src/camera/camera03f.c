
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


void CameraInitPosf ( CameraRecf *CPos )
{
  /* Setting initial (arbitrary) position of Camera, in which the camera */
  /* coordinate system is equal to the global system.                    */
  /* Call CameraInitFrame before this function.                          */
  if ( CPos->parallel ) {
    CPos->vd.para.diag = 1.0;
    CPos->vd.para.dim_case = 0;
  }
  else {
    CPos->vd.persp.f = 1.0;
  }
  memset((char *)(&CPos->position), 0, sizeof(point3f));
  CPos->theta = 0.0;
  CPos->psi = 0.0;
  CPos->phi = 0.0;
  CPos->g_centre = CPos->position;
  CPos->c_centre = CPos->position;
  CPos->c_fixed = false;
  CameraSetMappingf ( CPos );
} /*CameraInitPosf*/

void CameraSetRotCentref ( CameraRecf *CPos, point3f *centre,
                           boolean global_coord, boolean global_fixed )
{
  /* This procedure sets the centre point of camera rotations.    */
  /* If global_coord, then it is specified in global coordinates, */
  /* else in camera. if global_fixed, then it is fixed in global  */
  /* coordinate system, else moves together with the camera.      */
  if (global_coord) {
    CPos->g_centre = *centre;
    TransPoint3f ( &CPos->CTr, &CPos->g_centre, &CPos->c_centre );
  } else {
    SetPoint3f ( &CPos->c_centre, centre->x * CPos->xscale,
	         centre->y * CPos->yscale, centre->z );
    TransPoint3f ( &CPos->CTrInv, &CPos->c_centre, &CPos->g_centre );
  }
  CPos->c_fixed = (boolean)(!global_fixed);
} /*CameraSetRotCentref*/

void _CameraUpdateRotCentref ( CameraRecf *CPos )
{
  if (CPos->c_fixed)
    TransPoint3f ( &CPos->CTrInv, &CPos->c_centre, &CPos->g_centre );
  else
    TransPoint3f ( &CPos->CTr, &CPos->g_centre, &CPos->c_centre );
} /*_CameraUpdateRotCentref*/

