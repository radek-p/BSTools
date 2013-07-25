
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

void StereoMoveGd ( StereoRecd *Stereo, vector3d *v )
{
  AddVector3d ( &Stereo->position, v, &Stereo->position );
  StereoSetMappingd ( Stereo );
  _StereoUpdateRotCentred ( Stereo );
} /*StereoMoveGd*/

void StereoMoveCd ( StereoRecd *Stereo, vector3d *v )
{
  trans3d  t;
  vector3d vv;

  IdentTrans3d ( &t );
  EulerRotTrans3d ( &t, Stereo->left.phi, Stereo->left.theta, Stereo->left.psi );
  TransVector3d ( &t, v, &vv );
  AddVector3d ( &Stereo->position, &vv, &Stereo->position );
  StereoSetMappingd ( Stereo );
  _StereoUpdateRotCentred ( Stereo );
} /*StereoMoveCd*/

void StereoRotGd ( StereoRecd *Stereo, double _psi, double _theta, double _phi )
{
  trans3d t;

  IdentTrans3d ( &t );
  EulerRotTrans3d ( &t, -_psi, -_theta, -_phi );
  ShiftTrans3d ( &t, -Stereo->left.g_centre.x, -Stereo->left.g_centre.y,
                 -Stereo->left.g_centre.z );
  ShiftTrans3d ( &t, Stereo->position.x, Stereo->position.y,
                 Stereo->position.z );
  EulerRotTrans3d ( &t, _phi, _theta, _psi );
  ShiftTrans3d ( &t, Stereo->left.g_centre.x, Stereo->left.g_centre.y,
                 Stereo->left.g_centre.z );
  SetPoint3d ( &Stereo->position, t.U0.a14, t.U0.a24, t.U0.a34 );
  CompEulerRotd ( Stereo->left.phi, Stereo->left.theta, Stereo->left.psi,
                  _phi, _theta, _psi,
                  &Stereo->left.phi, &Stereo->left.theta, &Stereo->left.psi );
  StereoSetMappingd ( Stereo );
} /*StereoRotGd*/

void StereoRotVGd ( StereoRecd *Stereo, vector3d *v, double angle )
{
    /* turns the stereo around an arbitrary axis, parallel to the vector v,    */
    /* which is specified in the global coordinate system and must be nonzero. */

  double psi, theta, phi;

  FindRotVEulerd ( v, angle, &psi, &theta, &phi );
  StereoRotGd ( Stereo, psi, theta, phi );
} /*StereoRotVGd*/

void StereoRotCd ( StereoRecd *Stereo, double _psi, double _theta, double _phi )
{
  trans3d t;
  point3d c;

                           /* c - centre in stereo coordinates */
  IdentTrans3d ( &t );
  ShiftTrans3d ( &t, -Stereo->position.x, -Stereo->position.y, -Stereo->position.z );
  EulerRotTrans3d ( &t, -Stereo->left.psi, -Stereo->left.theta, -Stereo->left.phi );
  TransPoint3d ( &t, &Stereo->left.g_centre, &c );

  IdentTrans3d ( &t );
  EulerRotTrans3d ( &t, -Stereo->left.psi, -Stereo->left.theta, -Stereo->left.phi );
  EulerRotTrans3d ( &t, -_psi, -_theta, -_phi );
  ShiftTrans3d    ( &t, -c.x, -c.y, -c.z );
  EulerRotTrans3d ( &t, _phi, _theta, _psi );
  ShiftTrans3d    ( &t, c.x, c.y, c.z );
  EulerRotTrans3d ( &t, Stereo->left.phi, Stereo->left.theta, Stereo->left.psi );
  ShiftTrans3d    ( &t, Stereo->position.x, Stereo->position.y, Stereo->position.z );
  SetPoint3d ( &Stereo->position, t.U0.a14, t.U0.a24, t.U0.a34 );
  CompEulerRotd ( _phi, _theta, _psi,
                 Stereo->left.phi, Stereo->left.theta, Stereo->left.psi,
                 &Stereo->left.phi, &Stereo->left.theta, &Stereo->left.psi );
  StereoSetMappingd ( Stereo );
} /*StereoRotCd*/

void StereoRotVCd ( StereoRecd *Stereo, vector3d *v, double angle )
{
    /* turns the stereo around an arbitrary axis, parallel to the vector v,    */
    /* which is specified in the stereo coordinate system and must be nonzero. */
  double psi, theta, phi;

  FindRotVEulerd ( v, angle, &psi, &theta, &phi );
  StereoRotCd ( Stereo, psi, theta, phi );
} /*StereoRotVCd*/

void StereoZoomd ( StereoRecd *Stereo, double fchange )
{
  StereoSetDimd ( Stereo, fchange*Stereo->left.vd.persp.f,
                  Stereo->d, Stereo->l );
  StereoSetMappingd ( Stereo );
} /*StereoZoomd*/

