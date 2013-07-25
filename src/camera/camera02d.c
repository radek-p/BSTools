
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "camera.h"

#include "msgpool.h"


void CameraProjectPoint3d ( CameraRecd *CPos, const point3d *p, point3d *q )
{
  /* for a given point p (in global coordinates) this procedure finds its */
  /* projection (q.x,q.y) (in image coordinates). q.z is the depth of the */
  /* point p in the camera/image coordinate system and may be put to a    */
  /* z-buffer */
  if ( CPos->parallel ) {
    TransPoint3d ( &CPos->CTr, p, q );
  }
  else {
    TransPoint3d ( &CPos->CTr, p, q );
    if ( q->z > 0.0 ) {
      q->x = q->x/q->z + CPos->vd.persp.xi0;
      q->y = q->y/q->z + CPos->vd.persp.eta0;
    }
  }
  if ( CPos->upside )
    q->y = 2.0*CPos->ymin + CPos->height - q->y;
} /*CameraProjectPoint3d*/

void CameraUnProjectPoint3d ( CameraRecd *CPos, const point3d *p, point3d *q )
{
  point3d r;

  r = *p;
  if ( CPos->upside )
    r.y = 2.0*CPos->ymin + CPos->height - r.y;
  if ( CPos->parallel ) {
    TransPoint3d ( &CPos->CTrInv, &r, q );
  }
  else {
    r.x = (r.x-CPos->vd.persp.xi0)*r.z;
    r.y = (r.y-CPos->vd.persp.eta0)*r.z;
    TransPoint3d ( &CPos->CTrInv, &r, q );
  }
} /*CameraUnProjectPoint3d*/

void CameraProjectPoint2d ( CameraRecd *CPos, const point2d *p, point2d *q )
{
  point3d r;

  SetPoint3d ( &r, p->x, p->y, 0.0 );
  CameraProjectPoint3d ( CPos, &r, &r );
  SetPoint2d ( q, r.x, r.y );
} /*CameraProjectPoint2d*/

void CameraUnProjectPoint2d ( CameraRecd *CPos, const point2d *p, point2d *q )
{
  point3d r;

  SetPoint3d ( &r, p->x, p->y, 0.0 );
  CameraUnProjectPoint3d ( CPos, &r, &r );
  SetPoint2d ( q, r.x, r.y );
} /*CameraUnProjectPoint2d*/

void CameraProjectPoint3Rd ( CameraRecd *CPos, const point4d *p, point3d *q )
{
  point3d p0;

  Point4to3d ( p, &p0 );
  CameraProjectPoint3d ( CPos, &p0, q );
} /*CameraProjectPoint3Rd*/

void CameraUnProjectPoint3Rd ( CameraRecd *CPos, const point3d *p, double w,
                               point4d *q )
{
  point3d q0;

  CameraUnProjectPoint3d ( CPos, p, &q0 );
  Point3to4d ( &q0, w, q );
} /*CameraUnProjectPoint3Rd*/

void CameraProjectPoint2Rd ( CameraRecd *CPos, const point3d *p, point2d *q )
{
  point2d p0;

  Point3to2d ( p, &p0 );
  CameraProjectPoint2d ( CPos, &p0, q );
} /*CameraProjectPoint2Rd*/

void CameraUnProjectPoint2Rd ( CameraRecd *CPos, const point2d *p, double w,
                               point3d *q )
{
  point2d q0;

  CameraUnProjectPoint2d ( CPos, p, &q0 );
  Point2to3d ( &q0, w, q );
} /*CameraUnProjectPoint2Rd*/

void CameraProjectPoint3d2s ( CameraRecd *CPos, const point3d *p, point2s *q )
{
  point3d r;

  CameraProjectPoint3d ( CPos, p, &r );
  q->x = (short)(r.x+0.5);
  q->y = (short)(r.y+0.5);
} /*CameraProjectPoint3d2s*/

void CameraProjectPoint3Rd2s ( CameraRecd *CPos, const point4d *p, point2s *q )
{
  point3d p0;

  Point4to3d ( p, &p0 );
  CameraProjectPoint3d2s ( CPos, &p0, q );
} /*CameraProjectPoint3Rd2s*/

void CameraRayOfPixeld ( CameraRecd *CPos, double xi, double eta, ray3d *ray )
{
  /* for a given point (xi, eta) on the image, this procedure */
  /* finds the ray specified in the global coordinates        */
  point3d q;

  if ( CPos->upside )
    eta = 2.0*CPos->ymin + CPos->height - eta;
  if ( CPos->parallel ) {
    SetPoint3d ( &q, xi, eta, 0.0 );
    TransPoint3d ( &CPos->CTrInv, &q, &ray->p );
    SetPoint3d ( &q, xi, eta, 1.0 );
    TransPoint3d ( &CPos->CTrInv, &q, &q );
    SubtractPoints3d ( &q, &ray->p, &ray->v );
    NormalizeVector3d ( &ray->v );
  }
  else {
    ray->p = CPos->position;
    SetPoint3d ( &q, xi - CPos->vd.persp.xi0, eta - CPos->vd.persp.eta0, 1.0 );
    TransPoint3d ( &CPos->CTrInv, &q, &q );
    SubtractPoints3d ( &q, &ray->p, &ray->v );
    NormalizeVector3d ( &ray->v );
  }
} /*CameraRayOfPixeld*/

