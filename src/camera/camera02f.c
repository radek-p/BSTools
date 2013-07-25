
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


void CameraProjectPoint3f ( CameraRecf *CPos, const point3f *p, point3f *q )
{
  /* for a given point p (in global coordinates) this procedure finds its */
  /* projection (q.x,q.y) (in image coordinates). q.z is the depth of the */
  /* point p in the camera/image coordinate system and may be put to a    */
  /* z-buffer */
  if ( CPos->parallel ) {
    TransPoint3f ( &CPos->CTr, p, q );
  }
  else {
    TransPoint3f ( &CPos->CTr, p, q );
    if ( q->z > 0.0 ) {
      q->x = q->x/q->z + CPos->vd.persp.xi0;
      q->y = q->y/q->z + CPos->vd.persp.eta0;
    }
  }
  if ( CPos->upside )
    q->y = (float)(2.0*CPos->ymin + CPos->height - q->y);
} /*CameraProjectPoint3f*/

void CameraUnProjectPoint3f ( CameraRecf *CPos, const point3f *p, point3f *q )
{
  point3f r;

  r = *p;
  if ( CPos->upside )
    r.y = (float)(2.0*CPos->ymin + CPos->height - r.y);
  if ( CPos->parallel ) {
    TransPoint3f ( &CPos->CTrInv, &r, q );
  }
  else {
    r.x = (r.x-CPos->vd.persp.xi0)*r.z;
    r.y = (r.y-CPos->vd.persp.eta0)*r.z;
    TransPoint3f ( &CPos->CTrInv, &r, q );
  }
} /*CameraUnProjectPoint3f*/

void CameraProjectPoint2f ( CameraRecf *CPos, const point2f *p, point2f *q )
{
  point3f r;

  SetPoint3f ( &r, p->x, p->y, 0.0 );
  CameraProjectPoint3f ( CPos, &r, &r );
  SetPoint2f ( q, r.x, r.y );
} /*CameraProjectPoint2f*/

void CameraUnProjectPoint2f ( CameraRecf *CPos, const point2f *p, point2f *q )
{
  point3f r;

  SetPoint3f ( &r, p->x, p->y, 0.0 );
  CameraUnProjectPoint3f ( CPos, &r, &r );
  SetPoint2f ( q, r.x, r.y );
} /*CameraUnProjectPoint2f*/

void CameraProjectPoint3Rf ( CameraRecf *CPos, const point4f *p, point3f *q )
{
  point3f p0;

  Point4to3f ( p, &p0 );
  CameraProjectPoint3f ( CPos, &p0, q );
} /*CameraProjectPoint3Rd*/

void CameraUnProjectPoint3Rf ( CameraRecf *CPos, const point3f *p, float w,
                               point4f *q )
{
  point3f q0;

  CameraUnProjectPoint3f ( CPos, p, &q0 );
  Point3to4f ( &q0, w, q );
} /*CameraUnProjectPoint3Rd*/

void CameraProjectPoint2Rf ( CameraRecf *CPos, const point3f *p, point2f *q )
{
  point2f p0;

  Point3to2f ( p, &p0 );
  CameraProjectPoint2f ( CPos, &p0, q );
} /*CameraProjectPoint2Rd*/

void CameraUnProjectPoint2Rf ( CameraRecf *CPos, const point2f *p, float w,
                               point3f *q )
{
  point2f q0;

  CameraUnProjectPoint2f ( CPos, p, &q0 );
  Point2to3f ( &q0, w, q );
} /*CameraUnProjectPoint2Rd*/

void CameraProjectPoint3f2s ( CameraRecf *CPos, const point3f *p, point2s *q )
{
  point3f r;

  CameraProjectPoint3f ( CPos, p, &r );
  q->x = (short)(r.x+0.5);
  q->y = (short)(r.y+0.5);
} /*CameraProjectPoint3f2s*/

void CameraProjectPoint3Rf2s ( CameraRecf *CPos, const point4f *p, point2s *q )
{
  point3f p0;

  Point4to3f ( p, &p0 );
  CameraProjectPoint3f2s ( CPos, &p0, q );
} /*CameraProjectPoint3Rf2s*/

void CameraRayOfPixelf ( CameraRecf *CPos, float xi, float eta, ray3f *ray )
{
  /* for a given point (xi, eta) on the image, this procedure */
  /* finds the ray specified in the global coordinates        */
  point3f q;

  if ( CPos->upside )
    eta = (float)(2.0*CPos->ymin + CPos->height - eta);
  if ( CPos->parallel ) {
    SetPoint3f ( &q, xi, eta, 0.0 );
    TransPoint3f ( &CPos->CTrInv, &q, &ray->p );
    SetPoint3f ( &q, xi, eta, 1.0 );
    TransPoint3f ( &CPos->CTrInv, &q, &q );
    SubtractPoints3f ( &q, &ray->p, &ray->v );
    NormalizeVector3f ( &ray->v );
  }
  else {
    ray->p = CPos->position;
    SetPoint3f ( &q, xi - CPos->vd.persp.xi0, eta - CPos->vd.persp.eta0, 1.0 );
    TransPoint3f ( &CPos->CTrInv, &q, &q );
    SubtractPoints3f ( &q, &ray->p, &ray->v );
    NormalizeVector3f ( &ray->v );
  }
} /*CameraRayOfPixelf*/

