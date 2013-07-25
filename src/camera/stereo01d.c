
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
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

void StereoInitFramed ( StereoRecd *Stereo, boolean upside,
                        short width, short height, short xmin, short ymin,
                        double aspect, int ncplanes )
{
  CameraInitFramed ( &Stereo->left,  false, upside,
                     width, height, xmin, ymin, aspect, ncplanes );
  CameraInitFramed ( &Stereo->right, false, upside,
                     width, height, xmin, ymin, aspect, ncplanes );
} /*StereoInitFramed*/

void StereoSetDimd ( StereoRecd *Stereo, double f, double d, double l )
{
  double dx, w, h;

  Stereo->d = d;
  Stereo->l = l;
  Stereo->left.vd.persp.f  = f;
  Stereo->right.vd.persp.f = f;
  w = (double)Stereo->left.width;
  h = ((double)Stereo->left.height)/(Stereo->left.aspect);
  dx = 0.5*d*f/l*sqrt(w*w+h*h);
  Stereo->left.vd.persp.dxi0  = -dx;
  Stereo->right.vd.persp.dxi0 = dx;
} /*StereoSetDimd*/

void StereoSetMagd ( StereoRecd *Stereo, char mag )
{
  CameraSetMagd ( &Stereo->left,  mag );
  CameraSetMagd ( &Stereo->right, mag );
} /*StereoSetMagd*/

void StereoSetDepthRanged ( StereoRecd *Stereo, double zmin, double zmax )
{
  CameraSetDepthRanged ( &Stereo->left,  zmin, zmax );
  CameraSetDepthRanged ( &Stereo->right, zmin, zmax );
} /*StereoSetDepthRanged*/

void StereoSetMappingd ( StereoRecd *Stereo )
{
  point3d p;

        /* find the transformation STr from global coordinates to stereo */
  IdentTrans3d ( &Stereo->STr );
  ShiftTrans3d ( &Stereo->STr,
                 -Stereo->position.x, -Stereo->position.y, -Stereo->position.z );
  EulerRotTrans3d ( &Stereo->STr,
                    -Stereo->left.psi, -Stereo->left.theta, -Stereo->left.phi );

         /* calculate the position of each eye in the global coordinate system */
         /* find the transformation STrInv from stereo coordinates to global   */
  IdentTrans3d ( &Stereo->STrInv );
  EulerRotTrans3d ( &Stereo->STrInv,
                    Stereo->left.phi, Stereo->left.theta, Stereo->left.psi );
  ShiftTrans3d ( &Stereo->STrInv,
                 Stereo->position.x, Stereo->position.y, Stereo->position.z );
             /* find the position of the left eye */
  SetPoint3d ( &p, -0.5*Stereo->d, 0.0, 0.0 );
  TransPoint3d ( &Stereo->STrInv, &p, &Stereo->left.position );
             /* then the right eye */
  SetPoint3d ( &p, +0.5*Stereo->d, 0.0, 0.0 );
  TransPoint3d ( &Stereo->STrInv, &p, &Stereo->right.position );

         /* now, copy the other parameters */
  Stereo->right.vd.persp.f = Stereo->left.vd.persp.f;
  Stereo->right.psi   = Stereo->left.psi;
  Stereo->right.theta = Stereo->left.theta;
  Stereo->right.phi   = Stereo->left.phi;
         /* finally set the cameras of both eyes */
  CameraSetMappingd ( &Stereo->left );
  CameraSetMappingd ( &Stereo->right );
} /*StereoSetMappingd*/

void StereoInitPosd ( StereoRecd *Stereo )
{
  Stereo->left.vd.persp.f = 1.0; /* this is to prevent from the damage, caused by  */
  Stereo->d = 0.0;               /* uninitialised parameters. The values here make */
  Stereo->l = 1.0;               /* no sense, and should be given via SetStereoDim */
  memset ( &Stereo->position, 0, sizeof(point3d) );
  Stereo->left.theta = 0.0;
  Stereo->left.psi   = 0.0;
  Stereo->left.phi   = 0.0;
  Stereo->left.g_centre = Stereo->position;
  Stereo->left.c_centre = Stereo->position;
  Stereo->left.c_fixed  = false;
  StereoSetMappingd ( Stereo );
} /*StereoInitPosd*/

void StereoSetRotCentred ( StereoRecd *Stereo,
                           point3d *centre,
                           boolean global_coord, boolean global_fixed )
{
    /* This procedure sets the centre point of cameras rotations.   */
    /* If global_coord, then it is specified in global coordinates, */
    /* else in stereo. if global_fixed, then it is fixed in global  */
    /* coordinate system, else moves together with the cameras.     */

  if ( global_coord ) {
    Stereo->left.g_centre = *centre;
    TransPoint3d ( &Stereo->STr,
                   &Stereo->left.g_centre, &Stereo->left.c_centre );
  }
  else {
    Stereo->left.c_centre = *centre;
    TransPoint3d ( &Stereo->STrInv,
                   &Stereo->left.c_centre, &Stereo->left.g_centre );
  }
  Stereo->left.c_fixed = (boolean)(!global_fixed);
} /*StereoSetRotCentred*/

void _StereoUpdateRotCentred ( StereoRecd *Stereo )
{
  if ( Stereo->left.c_fixed )
    TransPoint3d ( &Stereo->STrInv,
                   &Stereo->left.c_centre, &Stereo->left.g_centre );
  else
    TransPoint3d ( &Stereo->STr,
                   &Stereo->left.g_centre, &Stereo->left.c_centre );
} /*_StereoUpdateRotCentred*/

