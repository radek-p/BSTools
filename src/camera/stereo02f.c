
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////// */
/* Defines two cameras for stereo viewing 3D scenes.                         */
/* Allows to manipulate with both cameras together.                          */
/* After setting the position of cameras, use procedures of the Camera unit  */
/* in order to project points to images corresponding to the left and right  */
/* eye.                                                                      */
/* ///////////////////////////////////////////////////////////////////////// */


#include <string.h>
#include <math.h>

#include "stereo.h"

#include "cprivate.h"
#include "msgpool.h"

void StereoMoveGf ( StereoRecf *Stereo, vector3f *v )
{
  AddVector3f ( &Stereo->position, v, &Stereo->position );
  StereoSetMappingf ( Stereo );
  _StereoUpdateRotCentref ( Stereo );
} /*StereoMoveGf*/

void StereoMoveCf ( StereoRecf *Stereo, vector3f *v )
{
  trans3f  t;
  vector3f vv;

  IdentTrans3f ( &t );
  EulerRotTrans3f ( &t, Stereo->left.phi, Stereo->left.theta, Stereo->left.psi );
  TransVector3f ( &t, v, &vv );
  AddVector3f ( &Stereo->position, &vv, &Stereo->position );
  StereoSetMappingf ( Stereo );
  _StereoUpdateRotCentref ( Stereo );
} /*StereoMoveCf*/

void StereoRotGf ( StereoRecf *Stereo, float _psi, float _theta, float _phi )
{
  trans3f t;

  IdentTrans3f ( &t );
  EulerRotTrans3f ( &t, -_psi, -_theta, -_phi );
  ShiftTrans3f ( &t, -Stereo->left.g_centre.x, -Stereo->left.g_centre.y,
                 -Stereo->left.g_centre.z );
  ShiftTrans3f ( &t, Stereo->position.x, Stereo->position.y,
                 Stereo->position.z );
  EulerRotTrans3f ( &t, _phi, _theta, _psi );
  ShiftTrans3f ( &t, Stereo->left.g_centre.x, Stereo->left.g_centre.y,
                 Stereo->left.g_centre.z );
  SetPoint3f ( &Stereo->position, t.U0.a14, t.U0.a24, t.U0.a34 );
  CompEulerRotf ( Stereo->left.phi, Stereo->left.theta, Stereo->left.psi,
                  _phi, _theta, _psi,
                  &Stereo->left.phi, &Stereo->left.theta, &Stereo->left.psi );
  StereoSetMappingf ( Stereo );
} /*StereoRotGf*/

void StereoRotVGf ( StereoRecf *Stereo, vector3f *v, float angle )
{
    /* turns the stereo around an arbitrary axis, parallel to the vector v,    */
    /* which is specified in the global coordinate system and must be nonzero. */

  float psi, theta, phi;

  FindRotVEulerf ( v, angle, &psi, &theta, &phi );
  StereoRotGf ( Stereo, psi, theta, phi );
} /*StereoRotVGf*/

void StereoRotCf ( StereoRecf *Stereo, float _psi, float _theta, float _phi )
{
  trans3f t;
  point3f c;

                           /* c - centre in stereo coordinates */
  IdentTrans3f ( &t );
  ShiftTrans3f ( &t, -Stereo->position.x, -Stereo->position.y, -Stereo->position.z );
  EulerRotTrans3f ( &t, -Stereo->left.psi, -Stereo->left.theta, -Stereo->left.phi );
  TransPoint3f ( &t, &Stereo->left.g_centre, &c );

  IdentTrans3f ( &t );
  EulerRotTrans3f ( &t, -Stereo->left.psi, -Stereo->left.theta, -Stereo->left.phi );
  EulerRotTrans3f ( &t, -_psi, -_theta, -_phi );
  ShiftTrans3f    ( &t, -c.x, -c.y, -c.z );
  EulerRotTrans3f ( &t, _phi, _theta, _psi );
  ShiftTrans3f    ( &t, c.x, c.y, c.z );
  EulerRotTrans3f ( &t, Stereo->left.phi, Stereo->left.theta, Stereo->left.psi );
  ShiftTrans3f    ( &t, Stereo->position.x, Stereo->position.y, Stereo->position.z );
  SetPoint3f ( &Stereo->position, t.U0.a14, t.U0.a24, t.U0.a34 );
  CompEulerRotf ( _phi, _theta, _psi,
                 Stereo->left.phi, Stereo->left.theta, Stereo->left.psi,
                 &Stereo->left.phi, &Stereo->left.theta, &Stereo->left.psi );
  StereoSetMappingf ( Stereo );
} /*StereoRotCf*/

void StereoRotVCf ( StereoRecf *Stereo, vector3f *v, float angle )
{
    /* turns the stereo around an arbitrary axis, parallel to the vector v,    */
    /* which is specified in the stereo coordinate system and must be nonzero. */
  float psi, theta, phi;

  FindRotVEulerf ( v, angle, &psi, &theta, &phi );
  StereoRotCf ( Stereo, psi, theta, phi );
} /*StereoRotVCf*/

void StereoZoomf ( StereoRecf *Stereo, float fchange )
{
  StereoSetDimf ( Stereo, fchange*Stereo->left.vd.persp.f,
                  Stereo->d, Stereo->l );
  StereoSetMappingf ( Stereo );
} /*StereoZoomf*/

