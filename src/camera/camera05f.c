
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

void CameraMoveCf ( CameraRecf *CPos, vector3f *v )
{
  trans3f t;
  vector3f vv;

  IdentTrans3f ( &t );
  EulerRotTrans3f ( &t, CPos->phi, CPos->theta, CPos->psi );
  TransVector3f ( &t, v, &vv );
  AddVector3f ( &CPos->position, &vv, &CPos->position );
  CameraSetMappingf ( CPos );
  _CameraUpdateRotCentref ( CPos );
}  /*CameraMoveCf*/

void CameraRotGf ( CameraRecf *CPos, float _psi, float _theta, float _phi )
{
  trans3f t;

  IdentTrans3f ( &t );
  EulerRotTrans3f ( &t, -_psi, -_theta, -_phi );
  ShiftTrans3f ( &t, -CPos->g_centre.x, -CPos->g_centre.y,
	         -CPos->g_centre.z );
  ShiftTrans3f ( &t, CPos->position.x, CPos->position.y,
                     CPos->position.z );
  EulerRotTrans3f ( &t, _phi, _theta, _psi );
  ShiftTrans3f ( &t, CPos->g_centre.x, CPos->g_centre.y,
                     CPos->g_centre.z );
  SetPoint3f ( &CPos->position, t.U0.a14, t.U0.a24, t.U0.a34 );
  CompEulerRotf ( CPos->phi, CPos->theta, CPos->psi, _phi, _theta, _psi,
	          &CPos->phi, &CPos->theta, &CPos->psi );
  CameraSetMappingf ( CPos );
} /*CameraRotGf*/

void CameraRotVGf ( CameraRecf *CPos, vector3f *v, float angle )
{
  /* turns the camera around an arbitrary axis, parallel to the vector v,    */
  /* which is specified in the global coordinate system and must be nonzero. */
  float psi, theta, phi;

  FindRotVEulerf ( v, angle, &psi, &theta, &phi );
  CameraRotGf ( CPos, psi, theta, phi );
} /*CameraRotVGf*/

void CameraRotCf ( CameraRecf *CPos, float _psi, float _theta, float _phi )
{
  trans3f t;
  point3f c;

  /* c - centre in camera unscaled coordinates */
  IdentTrans3f ( &t );
  ShiftTrans3f ( &t, -CPos->position.x, -CPos->position.y,
	             -CPos->position.z );
  EulerRotTrans3f ( &t, -CPos->psi, -CPos->theta, -CPos->phi );
  TransPoint3f ( &t, &CPos->g_centre, &c );
  IdentTrans3f ( &t );
  EulerRotTrans3f ( &t, -CPos->psi, -CPos->theta, -CPos->phi );
  EulerRotTrans3f ( &t, -_psi, -_theta, -_phi );
  ShiftTrans3f ( &t, -c.x, -c.y, -c.z );
  EulerRotTrans3f ( &t, _phi, _theta, _psi );
  ShiftTrans3f ( &t, c.x, c.y, c.z );
  EulerRotTrans3f ( &t, CPos->phi, CPos->theta, CPos->psi );
  ShiftTrans3f ( &t, CPos->position.x, CPos->position.y,
                     CPos->position.z );
  SetPoint3f ( &CPos->position, t.U0.a14, t.U0.a24, t.U0.a34 );
  CompEulerRotf ( _phi, _theta, _psi, CPos->phi, CPos->theta, CPos->psi,
                  &CPos->phi, &CPos->theta, &CPos->psi );
  CameraSetMappingf ( CPos );
} /*CameraRotCf*/

void CameraRotVCf ( CameraRecf *CPos, vector3f *v, float angle )
{
  /* turns the camera around an arbitrary axis, parallel to the vector v,    */
  /* which is specified in the camera coordinate system and must be nonzero. */
  float psi, theta, phi;

  FindRotVEulerf ( v, angle, &psi, &theta, &phi );
  CameraRotCf ( CPos, psi, theta, phi );
} /*CameraRotVCf*/

