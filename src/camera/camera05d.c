
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

void CameraMoveCd ( CameraRecd *CPos, vector3d *v )
{
  trans3d t;
  vector3d vv;

  IdentTrans3d ( &t );
  EulerRotTrans3d ( &t, CPos->phi, CPos->theta, CPos->psi );
  TransVector3d ( &t, v, &vv );
  AddVector3d ( &CPos->position, &vv, &CPos->position );
  CameraSetMappingd ( CPos );
  _CameraUpdateRotCentred ( CPos );
}  /*CameraMoveCd*/

void CameraRotGd ( CameraRecd *CPos, double _psi, double _theta, double _phi )
{
  trans3d t;

  IdentTrans3d ( &t );
  EulerRotTrans3d ( &t, -_psi, -_theta, -_phi );
  ShiftTrans3d ( &t, -CPos->g_centre.x, -CPos->g_centre.y,
	         -CPos->g_centre.z );
  ShiftTrans3d ( &t, CPos->position.x, CPos->position.y,
                     CPos->position.z );
  EulerRotTrans3d ( &t, _phi, _theta, _psi );
  ShiftTrans3d ( &t, CPos->g_centre.x, CPos->g_centre.y,
                     CPos->g_centre.z );
  SetPoint3d ( &CPos->position, t.U0.a14, t.U0.a24, t.U0.a34 );
  CompEulerRotd ( CPos->phi, CPos->theta, CPos->psi, _phi, _theta, _psi,
	          &CPos->phi, &CPos->theta, &CPos->psi );
  CameraSetMappingd ( CPos );
} /*CameraRotGd*/

void CameraRotVGd ( CameraRecd *CPos, vector3d *v, double angle )
{
  /* turns the camera around an arbitrary axis, parallel to the vector v,    */
  /* which is specified in the global coordinate system and must be nonzero. */
  double psi, theta, phi;

  FindRotVEulerd ( v, angle, &psi, &theta, &phi );
  CameraRotGd ( CPos, psi, theta, phi );
} /*CameraRotVGd*/

void CameraRotCd ( CameraRecd *CPos, double _psi, double _theta, double _phi )
{
  trans3d t;
  point3d c;

  /* c - centre in camera unscaled coordinates */
  IdentTrans3d ( &t );
  ShiftTrans3d ( &t, -CPos->position.x, -CPos->position.y,
	             -CPos->position.z );
  EulerRotTrans3d ( &t, -CPos->psi, -CPos->theta, -CPos->phi );
  TransPoint3d ( &t, &CPos->g_centre, &c );
  IdentTrans3d ( &t );
  EulerRotTrans3d ( &t, -CPos->psi, -CPos->theta, -CPos->phi );
  EulerRotTrans3d ( &t, -_psi, -_theta, -_phi );
  ShiftTrans3d ( &t, -c.x, -c.y, -c.z );
  EulerRotTrans3d ( &t, _phi, _theta, _psi );
  ShiftTrans3d ( &t, c.x, c.y, c.z );
  EulerRotTrans3d ( &t, CPos->phi, CPos->theta, CPos->psi );
  ShiftTrans3d ( &t, CPos->position.x, CPos->position.y,
                     CPos->position.z );
  SetPoint3d ( &CPos->position, t.U0.a14, t.U0.a24, t.U0.a34 );
  CompEulerRotd ( _phi, _theta, _psi, CPos->phi, CPos->theta, CPos->psi,
                  &CPos->phi, &CPos->theta, &CPos->psi );
  CameraSetMappingd ( CPos );
} /*CameraRotCd*/

void CameraRotVCd ( CameraRecd *CPos, vector3d *v, double angle )
{
  /* turns the camera around an arbitrary axis, parallel to the vector v,    */
  /* which is specified in the camera coordinate system and must be nonzero. */
  double psi, theta, phi;

  FindRotVEulerd ( v, angle, &psi, &theta, &phi );
  CameraRotCd ( CPos, psi, theta, phi );
} /*CameraRotVCd*/

