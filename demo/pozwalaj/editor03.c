
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>   
#include <stdio.h>
#include <math.h>  
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>  
#include <GL/glu.h> 
#include <GL/glx.h>

#include "pkvaria.h" 
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

/*#include "widgets.h"*/
#include "editor.h"  
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"


void GeomObjectSetupIniPoints ( char spdimen, boolean rational,
                                char *cpdimen,
                                int np, const double *pt4, double *pt )
{
  int  i, j, k;
  char cpdim;

  if ( rational ) {
    switch ( spdimen ) {
  case 2:
      *cpdimen = cpdim = 3;
      for ( i = j = k = 0;  i < np;  i++, j += 4, k += 3 )
        pt[k] = pt4[j],  pt[k+1] = pt4[j+1],  pt[k+2] = pt4[j+3];
      break;
  case 3:
      *cpdimen = cpdim = 4;
      memcpy ( pt, pt4, np*4*sizeof(double) );
      break;
  default:
      return;
    }
  }
  else {
    *cpdimen = cpdim = spdimen;
    switch ( spdimen ) {
  case 2:
      pkv_Selectd ( np, 2, 4, 2, pt4, pt );
      break;
  case 3:
      pkv_Selectd ( np, 3, 4, 3, pt4, pt );
      break;
  default:
      return;
    }
  }
} /*GeomObjectSetupIniPoints*/

void GeomObjectSetupWeightPoints ( char cpdimen, int np, double *cpoints,
                                   double *weightpoints )
{
  pkn_AddMatrixd ( 1, cpdimen*(np-1), 0, cpoints, 0, &cpoints[(int)cpdimen],
                   0, weightpoints );
/*
  pkn_MultMatrixNumd ( 1, cpdimen*(np-1), 0, weightpoints, 0.5,
                       0, weightpoints );
*/
} /*GeomObjectSetupWeightPoints*/

void GeomObjectUnmarkPoints ( int np, byte *mkcp )
{
  int i;

  for ( i = 0; i < np; i++ )
    mkcp[i] &= ~(char)(MASK_CP_MOVEABLE);
} /*GeomObjectUnmarkPoints*/

static void CVP ( char cpdimen, boolean rational, double *cp, point3d *p )
{
  if ( rational ) {
    switch ( cpdimen ) {
  case 2:
      SetPoint3d ( p, cp[0]/cp[1], 0.0, 0.0 );
      break;
  case 3:
      SetPoint3d ( p, cp[0]/cp[2], cp[1]/cp[2], 0.0 );
      break;
  case 4:
      SetPoint3d ( p, cp[0]/cp[3], cp[1]/cp[3], cp[2]/cp[3] );
      break;
  default:
      return;
    }
  }
  else {
    switch ( cpdimen ) {
  case 2:
      SetPoint3d ( p, cp[0], cp[1], 0.0 );
      break;
  case 3:
  case 4:
      memcpy ( p, cp, sizeof(point3d) );
      break;
  default:
      return;
    }
  }
} /*CVP*/

void GeomObjectFindBBox ( char cpdimen, boolean rational, int np, double *cp,
                          Box3d *box )
{
  int     i, k;
  point3d p;

  if ( np > 0 ) {
    CVP ( cpdimen, rational, cp, &p );
    box->x0 = box->x1 = p.x;
    box->y0 = box->y1 = p.y;
    box->z0 = box->z1 = p.z;
    for ( i = 1, k = cpdimen;  i < np;  i++, k += cpdimen ) {
      CVP ( cpdimen, rational, &cp[k], &p );
      if ( p.x < box->x0 )      box->x0 = p.x;
      else if ( p.x > box->x1 ) box->x1 = p.x;
      if ( p.y < box->y0 )      box->y0 = p.y;
      else if ( p.y > box->y1 ) box->y1 = p.y;
      if ( p.z < box->z0 )      box->z0 = p.z;
      else if ( p.z > box->z1 ) box->z1 = p.z;
    }
  }
} /*GeomObjectFindBBox*/

void GeomObjectSetCurrentPointPtr ( double *currpt, char dim, boolean rational )
{
  current_point = currpt;
  current_point_dim = dim;
  current_point_rational = rational;
} /*GeomObjectSetCurrentPointPtr*/

boolean GeomObjectFindNearestPoint ( char cpdimen, char spdimen,
                                     int np, double *cp, byte *mkcp,
                                     byte mask, CameraRecd *CPos,
                                     short x, short y, int *dist )
{
  int     i, k;
  int     d0, d1;
  point3d p, q;
  char    m;
  boolean ok;

  d0 = *dist;
  for ( i = k = 0;  i < np;  i++, k += cpdimen ) {
    m = mkcp != NULL ? mkcp[i] & mask : true;
    if ( m ) {
      switch ( cpdimen ) {
    case 2:
        SetPoint3d ( &p, cp[k], cp[k+1], 0.0 );
        ok = CameraClipPoint3d ( CPos, &p, &q );
        break;
    case 3:
        if ( spdimen == 2 ) {
          SetPoint3d ( &p, cp[k]/cp[k+2], cp[k+1]/cp[k+2], 0.0 );
          ok = CameraClipPoint3d ( CPos, &p, &q );
        }
        else
          ok = CameraClipPoint3d ( CPos, (point3d*)&cp[k], &q );
        break;
    case 4:
        Point4to3d ( (point4d*)&cp[k], &p );
        ok = CameraClipPoint3d ( CPos, &p, &q );
        break;
    default:
        return false;
      }
      if ( ok ) {
        d1 = (int)(fabs(x-q.x) + fabs(y-q.y));
        if ( d1 < d0 ) {
          d0 = d1;
          current_point_ind = i;
        }
      }
    }
  }
  if ( d0 < *dist ) {
    GeomObjectSetCurrentPointPtr ( &cp[current_point_ind*cpdimen],
                                   cpdimen, cpdimen > spdimen );
    *dist = d0;
    return true;
  }
  else
    return false;
} /*GeomObjectFindNearestPoint*/

static boolean _GeomObjectFindObjCPoint ( char spdimen, CameraRecd *CPos,
                              short x, short y, int *dist, geom_object *obj )
{
  switch ( obj->obj_type ) {
case GO_BEZIER_CURVE:
    if ( GeomObjectBezierCurveFindCPoint ( (GO_BezierCurve*)obj,
                                           CPos, x, y, dist ) ) {
      currentp_go = obj;
      return true;
    }
    break;
case GO_BSPLINE_CURVE:
    if ( GeomObjectBSplineCurveFindCPoint ( (GO_BSplineCurve*)obj,
                                            CPos, x, y, dist ) ) {
      currentp_go = obj;
      return true;
    }
    break;
case GO_BEZIER_PATCH:
    if ( GeomObjectBezierPatchFindCPoint ( (GO_BezierPatch*)obj,
                                           CPos, x, y, dist ) ) {
      currentp_go = obj;
      return true;
    }
    break;
case GO_BSPLINE_PATCH:
    if ( GeomObjectBSplinePatchFindCPoint ( (GO_BSplinePatch*)obj,
                                            CPos, x, y, dist ) ) {
      currentp_go = obj;
      return true;
    }
    break;
case GO_BSPLINE_MESH:
    if ( GeomObjectBSplineMeshFindCPoint ( (GO_BSplineMesh*)obj,
                                           CPos, x, y, dist ) ) {
      currentp_go = obj;
      return true;
    }
    break;
case GO_BSPLINE_HOLE:
    if ( GeomObjectBSplineHoleFindCPoint ( (GO_BSplineHole*)obj,
                                           CPos, x, y, dist ) ) {
      currentp_go = obj;
      return true;
    }
    break;
default:
    break;
  }
  return false;
} /*_GeomObjectFindObjCPoint*/

boolean GeomObjectFindObjCPoint ( char spdimen, CameraRecd *CPos, short x, short y,   
                                  geom_object *obj )
{
  int dist;

  currentp_go = NULL;
  dist = MAXPIXDIST+1;
  return _GeomObjectFindObjCPoint ( spdimen, CPos, x, y, &dist, obj );
} /*GeomObjectFindObjCPoint*/

boolean GeomObjectFindCPoint ( char spdimen, CameraRecd *CPos, short x, short y )
{
  geom_object *go;
  int         dist;

  currentp_go = NULL;
  dist = MAXPIXDIST+1;
  for ( go = first_go; go; go = go->next )
    if ( (go == current_go || go->active) && go->spdimen == spdimen )
      _GeomObjectFindObjCPoint ( spdimen, CPos, x, y, &dist, go );
  return currentp_go != NULL;
} /*GeomObjectFindCPoint*/

void GeomObjectSetPoint ( CameraRecd *CPos, short x, short y,
                          char cpdimen, char spdimen, double *pt )
{
  point3d q3;
  point2d q2;

  switch ( cpdimen ) {
case 2:
    SetPoint2d ( &q2, (double)x, (double)y );
    CameraUnProjectPoint2d ( CPos, &q2, (point2d*)pt );
    break;
case 3:
    switch ( spdimen ) {
  case 2:
      SetPoint2d ( &q2, (double)x, (double)y );
      CameraUnProjectPoint2Rd ( CPos, &q2, pt[2], (point3d*)pt );
      break;
  case 3:
      CameraProjectPoint3d ( CPos, (point3d*)pt, &q3 );
      q3.x = (double)x;  q3.y = (double)y;
      CameraUnProjectPoint3d ( CPos, &q3, (point3d*)pt );
      break;
  default:
      break;
    }
    break;
case 4:
    CameraProjectPoint3Rd ( CPos, (point4d*)pt, &q3 );
    q3.x = (double)x;  q3.y = (double)y;
    CameraUnProjectPoint3Rd ( CPos, &q3, pt[3], (point4d*)pt );
    break;
default:
    break;
  }
} /*GeomObjectSetPoint*/

static double _GeomObjectWeightRatio ( double d, double e )
{
#define MMIN   0.01
#define MMAX 100.0
  double m;

  if ( fabs(d) < 2.0 ) return 1.0;
  else if ( e <= 0.0 ) return MMIN;
  else if ( e >= d )   return MMAX;
  else {
    m = e/(d-e);
    m = max ( m, MMIN );
    return min ( m, MMAX );
  }
#undef MMIN
#undef MMAX
} /*_GeomObjectWeightRatio*/

void GeomObjectSetWeightPoint ( CameraRecd *CPos, short x, short y,
                                char cpdimen, double *pt, double *wpt )
{
  vector2d v, w;
  point2d  q;
  double   m, nw, r;

  switch ( cpdimen ) {
case 3: {
      point2d p0, p1;

      CameraProjectPoint2Rd ( CPos, (point3d*)pt, &p0 );
      CameraProjectPoint2Rd ( CPos, (point3d*)&pt[3], &p1 );
      SetPoint2d ( &q, (double)x, (double)y );
      SubtractPoints2d ( &p1, &p0, &v );
      SubtractPoints2d ( &q, &p0, &w );
      m = _GeomObjectWeightRatio ( DotProduct2d ( &v, &v ),
                                   DotProduct2d ( &v, &w ) );
      nw = m*pt[2];
      r = nw/pt[5];
      MultVector2d ( r, (vector2d*)&pt[3], (vector2d*)&pt[3] );
      pt[5] = nw;
      GeomObjectSetupWeightPoints ( 3, 2, pt, wpt );
    }
    break;
case 4: {
      point3d p0, p1;

      CameraProjectPoint3Rd ( CPos, (point4d*)pt, &p0 );
      CameraProjectPoint3Rd ( CPos, (point4d*)&pt[4], &p1 );
      SetPoint2d ( &q, (double)x, (double)y );
      SubtractPoints2d ( (point2d*)&p1, (point2d*)&p0, &v );
      SubtractPoints2d ( &q, (point2d*)&p0, &w );
      m = _GeomObjectWeightRatio ( DotProduct2d ( &v, &v ),
                                   DotProduct2d ( &v, &w ) );
      nw = m*pt[3];
      r = nw/pt[7];
      MultVector3d ( r, (vector3d*)&pt[4], (vector3d*)&pt[4] );
      pt[7] = nw;
      GeomObjectSetupWeightPoints ( 4, 2, pt, wpt );
    }
    break;
default:
    break;
  }
} /*GeomObjectSetWeightPoint*/

void GeomObjectSetCPoint ( CameraRecd *CPos, short x, short y )
{
  if ( currentp_go ) {
    switch ( currentp_go->obj_type ) {
  case GO_BEZIER_CURVE:
      GeomObjectBezierCurveSetCPoint ( (GO_BezierCurve*)currentp_go, CPos, x, y );
      break;
  case GO_BSPLINE_CURVE:
      GeomObjectBSplineCurveSetCPoint ( (GO_BSplineCurve*)currentp_go, CPos, x, y );
      break;
  case GO_BEZIER_PATCH:
      GeomObjectBezierPatchSetCPoint ( (GO_BezierPatch*)currentp_go, CPos, x, y );
      break;
  case GO_BSPLINE_PATCH:
      GeomObjectBSplinePatchSetCPoint ( (GO_BSplinePatch*)currentp_go, CPos, x, y );
      break;
  case GO_BSPLINE_MESH:
      GeomObjectBSplineMeshSetCPoint ( (GO_BSplineMesh*)currentp_go, CPos, x, y );
      break;
  case GO_BSPLINE_HOLE:
      GeomObjectBSplineHoleSetCPoint ( (GO_BSplineHole*)currentp_go, CPos, x, y );
      break;
  default:
      break;
    }
  }
} /*GeomObjectSetCPoint*/

void GeomObjectTransformPoint ( char cpdimen, char spdimen,
                                double *scp, double *cp, trans3d *tr )
{
  point3d p, q;

  switch ( cpdimen ) {
case 2:
    SetPoint3d ( &p, scp[0], scp[1], 0.0 );
    TransPoint3d ( tr, &p, &q );
    cp[0] = q.x;  cp[1] = q.y;
    break;
case 3:
    if ( spdimen == 3 )
      TransPoint3d ( tr, (point3d*)scp, (point3d*)cp );
    else {
      SetPoint3d ( &p, scp[0]/scp[2], scp[1]/scp[2], 0.0 );
      TransPoint3d ( tr, &p, &q );
      cp[0] = q.x*scp[2];  cp[1] = q.y*scp[2];  cp[2] = scp[2];
    }
    break;
case 4:
    Trans3Point4d ( tr, (point4d*)scp, (point4d*)cp );
    break;
default:
    break;
  }
} /*GeomObjectTransformPoint*/

boolean GeomObjectPointInBox ( char cpdimen, char spdimen, double *pt,
                               CameraRecd *CPos, Box2s *box )
{
  point3d p, q;

  switch ( cpdimen ) {
case 2:
    SetPoint3d ( &p, pt[0], pt[1], 0.0 );
    break;
case 3:
    if ( spdimen == 3 )
      memcpy ( &p, pt, sizeof(point3d) );
    else
      SetPoint3d ( &p, pt[0]/pt[2], pt[1]/pt[2], 0.0 );
    break;
case 4:
    Point4to3d ( (point4d*)pt, &p );
    break;
default:
    return false;
  }
  if ( CameraClipPoint3d ( CPos, &p, &q ) )
    return q.x >= box->x0 && q.x <= box->x1 &&
           q.y >= box->y0 && q.y <= box->y1;
  else
    return false;
} /*GeomObjectPointInBox*/

void GeomObjectMarkPoints ( char cpdimen, char spdimen,
                            int np, byte *mkcp, double *cp,
                            CameraRecd *CPos, Box2s *box,
                            byte mask, boolean clear )
{
  int i, k;

  for ( i = k = 0;  i < np;  i++, k += cpdimen )
    if ( GeomObjectPointInBox ( cpdimen, spdimen, &cp[k], CPos, box ) ) {
      if ( clear )
        mkcp[i] &= ~mask;
      else
        mkcp[i] |= mask;
    }
} /*GeomObjectMarkPoints*/

void GeomObjectMarkPoint ( int np, byte *mkcp, byte mask, boolean clear )
{
  if ( current_point_ind >= 0 && current_point_ind < np ) {
    if ( clear )
      mkcp[current_point_ind] &= ~mask;
    else
      mkcp[current_point_ind] |= mask;
  }
} /*GeomObjectMarkPoint*/

void GeomObjectMarkCPoints ( char spdimen, CameraRecd *CPos, Box2s *box,
                             boolean clear )
{
  geom_object *go;

  if ( box->x1-box->x0 < 3 && box->y1-box->y0 < 3 ) {
        /* the box is too small - assume that one point has been clicked */
    if ( GeomObjectFindCPoint ( spdimen, CPos, box->x0, box->y0 ) )
      switch ( currentp_go->obj_type ) {
    case GO_BEZIER_CURVE:
        GeomObjectBezierCurveMarkCPoint ( (GO_BezierCurve*)currentp_go,
                                          marking_mask, clear );
        break;
    case GO_BSPLINE_CURVE:
        GeomObjectBSplineCurveMarkCPoint ( (GO_BSplineCurve*)currentp_go,
                                           marking_mask, clear );
        break;
    case GO_BEZIER_PATCH:
        GeomObjectBezierPatchMarkCPoint ( (GO_BezierPatch*)currentp_go,
                                          marking_mask, clear );
        break;
    case GO_BSPLINE_PATCH:
        GeomObjectBSplinePatchMarkCPoint ( (GO_BSplinePatch*)currentp_go,
                                           marking_mask, clear );
        break;
    case GO_BSPLINE_MESH:
        GeomObjectBSplineMeshMarkCPoint ( (GO_BSplineMesh*)currentp_go,
                                          marking_mask, clear );
        break;
    case GO_BSPLINE_HOLE:
        GeomObjectBSplineHoleMarkCPoint ( (GO_BSplineHole*)currentp_go,
                                          marking_mask, clear );
        break;
    default:
        break;
      }
  }
  else {
    for ( go = first_go; go; go = go->next )
      if ( (go == current_go || go->active) && go->spdimen == spdimen ) {
        switch ( go->obj_type ) {
      case GO_BEZIER_CURVE:
          GeomObjectBezierCurveMarkCPoints ( (GO_BezierCurve*)go, CPos, box,
                                              marking_mask, clear );
          break;
      case GO_BSPLINE_CURVE:
          GeomObjectBSplineCurveMarkCPoints ( (GO_BSplineCurve*)go, CPos, box,
                                              marking_mask, clear );
          break;
      case GO_BEZIER_PATCH:
          GeomObjectBezierPatchMarkCPoints ( (GO_BezierPatch*)go, CPos, box,
                                              marking_mask, clear );
          break;
      case GO_BSPLINE_PATCH:
          GeomObjectBSplinePatchMarkCPoints ( (GO_BSplinePatch*)go, CPos, box,
                                              marking_mask, clear );
          break;
      case GO_BSPLINE_MESH:
          GeomObjectBSplineMeshMarkCPoints ( (GO_BSplineMesh*)go, CPos, box,
                                             marking_mask, clear );
          break;
      case GO_BSPLINE_HOLE:
          GeomObjectBSplineHoleMarkCPoints ( (GO_BSplineHole*)go, CPos, box,
                                             marking_mask, clear );
          break;
      default:
          break;
        }
      }
  }
} /*GeomObjectMarkCPoints*/

void GeomObjectSetMarkingMask ( boolean bits[5] )
{
  int i;
  geom_object *go;

  marking_mask = 0;
  for ( i = 0; i < 5; i++ )
    if ( bits[i] )
      marking_mask |= MASK_CP_MARKED_0 << i;
  for ( go = first_go; go; go = go->next )
    go->dlistmask &= ~DLISTMASK_CPOINTS;
} /*GeomObjectChangeMarkingMask*/

void GeomObjectTransformPoints ( char cpdimen, char spdimen,
                                 int np, byte *mkcp, byte mask,
                                 double *savedp, double *cp, trans3d *tr )
{
  int i, k;

  for ( i = k = 0;  i < np;  i++, k += cpdimen )
    if ( mkcp[i] & mask )
  GeomObjectTransformPoint ( cpdimen, spdimen, &savedp[k], &cp[k], tr );
} /*GeomObjectTransformPoints*/

boolean GeomObjectSaveCPoints ( char spdimen )
{
  geom_object *go;
  boolean     isok;

  isok = true;
  for ( go = first_go; go; go = go->next )
    if ( (go == current_go || go->active) && go->spdimen == spdimen ) {
      switch ( go->obj_type ) {
    case GO_BEZIER_CURVE:
        isok &= GeomObjectBezierCurveSaveCPoints ( (GO_BezierCurve*)go );
        break;
    case GO_BSPLINE_CURVE:
        isok &= GeomObjectBSplineCurveSaveCPoints ( (GO_BSplineCurve*)go );
        break;
    case GO_BEZIER_PATCH:
        isok &= GeomObjectBezierPatchSaveCPoints ( (GO_BezierPatch*)go );
        break;
    case GO_BSPLINE_PATCH:
        isok &= GeomObjectBSplinePatchSaveCPoints ( (GO_BSplinePatch*)go );
        break;
    case GO_BSPLINE_MESH:
        isok &= GeomObjectBSplineMeshSaveCPoints ( (GO_BSplineMesh*)go );
        break;
    case GO_BSPLINE_HOLE:
        isok &= GeomObjectBSplineHoleSaveCPoints ( (GO_BSplineHole*)go );
        break;
    default:
        break;
      }
    }
  return isok;
} /*GeomObjectSaveCPoints*/

void GeomObjectUndoLastTransformation ( char spdimen )
{
  geom_object *go;

  for ( go = first_go; go; go = go->next )
    if ( (go == current_go || go->active) && go->spdimen == spdimen ) {
      switch ( go->obj_type ) {
    case GO_BEZIER_CURVE:
        GeomObjectBezierCurveUndoLastTransformation ( (GO_BezierCurve*)go );
        break;
    case GO_BEZIER_PATCH:
        GeomObjectBezierPatchUndoLastTransformation ( (GO_BezierPatch*)go );
        break;
    case GO_BSPLINE_CURVE:
        GeomObjectBSplineCurveUndoLastTransformation ( (GO_BSplineCurve*)go );
        break;
    case GO_BSPLINE_PATCH:
        GeomObjectBSplinePatchUndoLastTransformation ( (GO_BSplinePatch*)go );
        break;
    case GO_BSPLINE_MESH:
        GeomObjectBSplineMeshUndoLastTransformation ( (GO_BSplineMesh*)go );
        break;
    case GO_BSPLINE_HOLE:
        GeomObjectBSplineHoleUndoLastTransformation ( (GO_BSplineHole*)go );
        break;
      }
    }
} /*GeomObjectUndoLastTransformation*/

void GeomObjectTransformCPoints2D ( trans2d *tr, byte mask )
{
  trans3d     tr3;
  geom_object *go;

  IdentTrans3d ( &tr3 );
  tr3.U0.a11 = tr->U0.a11;
  tr3.U0.a12 = tr->U0.a12;
  tr3.U0.a14 = tr->U0.a13;
  tr3.U0.a21 = tr->U0.a21;
  tr3.U0.a22 = tr->U0.a22;
  tr3.U0.a24 = tr->U0.a23;
  tr3.U1.detsgn = tr->U1.detsgn;
  for ( go = first_go; go; go = go->next )
    if ( (go == current_go || go->active) && go->spdimen == 2 ) {
      switch ( go->obj_type ) {
    case GO_BEZIER_CURVE:
        GeomObjectBezierCurveTransformCPoints ( (GO_BezierCurve*)go, &tr3, mask );
        break;
    case GO_BSPLINE_CURVE:
        GeomObjectBSplineCurveTransformCPoints ( (GO_BSplineCurve*)go, &tr3, mask );
        break;
    case GO_BEZIER_PATCH:
        GeomObjectBezierPatchTransformCPoints ( (GO_BezierPatch*)go, &tr3, mask );
        break;
    case GO_BSPLINE_PATCH:
        GeomObjectBSplinePatchTransformCPoints ( (GO_BSplinePatch*)go, &tr3, mask );
        break;
    case GO_BSPLINE_MESH:
        GeomObjectBSplineMeshTransformCPoints ( (GO_BSplineMesh*)go, &tr3, mask );
        break;
    case GO_BSPLINE_HOLE:
        GeomObjectBSplineHoleTransformCPoints ( (GO_BSplineHole*)go, &tr3, mask );
        break;
    default:
        break;
      }
    }    
} /*GeomObjectTransformCPoints2D*/

void GeomObjectTransformCPoints3D ( trans3d *tr, byte mask )
{
  geom_object *go;

  for ( go = first_go; go; go = go->next )
    if ( (go == current_go || go->active) && go->spdimen == 3 ) {
      switch ( go->obj_type ) {
    case GO_BEZIER_CURVE:
        GeomObjectBezierCurveTransformCPoints ( (GO_BezierCurve*)go, tr, mask );
        break;
    case GO_BSPLINE_CURVE:
        GeomObjectBSplineCurveTransformCPoints ( (GO_BSplineCurve*)go, tr, mask );
        break;
    case GO_BEZIER_PATCH:
        GeomObjectBezierPatchTransformCPoints ( (GO_BezierPatch*)go, tr, mask );
        break;
    case GO_BSPLINE_PATCH:
        GeomObjectBSplinePatchTransformCPoints ( (GO_BSplinePatch*)go, tr, mask );
        break;
    case GO_BSPLINE_MESH:
        GeomObjectBSplineMeshTransformCPoints ( (GO_BSplineMesh*)go, tr, mask );
        break;
    case GO_BSPLINE_HOLE:
        GeomObjectBSplineHoleTransformCPoints ( (GO_BSplineHole*)go, tr, mask );
        break;
    default:
        break;
      }
    }    
} /*GeomObjectTransformCPoints3D*/

boolean GeomObjectGetPointCoord ( geom_object *go, int p,
                                  int *spdimen, int *cpdimen, double **pc )
{
  switch ( go->obj_type ) {
case GO_BEZIER_CURVE:
    return GeomObjectBezierCurveGetPointCoord ( (GO_BezierCurve*)go, p,
                                                spdimen, cpdimen, pc );
case GO_BSPLINE_CURVE:
    return GeomObjectBSplineCurveGetPointCoord ( (GO_BSplineCurve*)go, p,
                                                 spdimen, cpdimen, pc );
case GO_BEZIER_PATCH:
    return GeomObjectBezierPatchGetPointCoord ( (GO_BezierPatch*)go, p,
                                                 spdimen, cpdimen, pc );
case GO_BSPLINE_PATCH:
    return GeomObjectBSplinePatchGetPointCoord ( (GO_BSplinePatch*)go, p,
                                                 spdimen, cpdimen, pc );
case GO_BSPLINE_MESH:
    return GeomObjectBSplineMeshGetPointCoord ( (GO_BSplineMesh*)go, p,
                                                spdimen, cpdimen, pc );
case GO_BSPLINE_HOLE:
    return GeomObjectBSplineHoleGetPointCoord ( (GO_BSplineHole*)go, p,
                                                spdimen, cpdimen, pc );
default:
    return false;
  }
} /*GeomObjectGetPointCoord*/

void GeomObjectDisplayInfoText ( geom_object *obj )
{
  switch ( obj->obj_type ) {
case GO_BEZIER_CURVE:
    GeomObjectBezierCurveDisplayInfoText ( (GO_BezierCurve*)obj );
    break;
case GO_BSPLINE_CURVE:
    GeomObjectBSplineCurveDisplayInfoText ( (GO_BSplineCurve*)obj );
    break;
case GO_BEZIER_PATCH:
    GeomObjectBezierPatchDisplayInfoText ( (GO_BezierPatch*)obj );
    break;
case GO_BSPLINE_PATCH:
    GeomObjectBSplinePatchDisplayInfoText ( (GO_BSplinePatch*)obj );
    break;
case GO_BSPLINE_MESH:
    GeomObjectBSplineMeshDisplayInfoText ( (GO_BSplineMesh*)obj );
    break;
case GO_BSPLINE_HOLE:
    GeomObjectBSplineHoleDisplayInfoText ( (GO_BSplineHole*)obj );
    break;
default:
    break;
  }
} /*GeomObjectDisplayInfoText*/

