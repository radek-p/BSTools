
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
#include "pkgeomclip.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "editor.h"
#include "editor_bsm.h"


void GeomObjectBSplineMeshFindBBox ( GO_BSplineMesh *obj, Box3d *box )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  GeomObjectFindBBox ( obj->me.cpdimen, obj->rational,  
                       obj->nv, obj->meshvpc, box );
} /*GeomObjectBSplineMeshFindBBox*/

boolean GeomObjectBSplineMeshFindCPoint ( GO_BSplineMesh *obj,
                                          CameraRecd *CPos, short x, short y,
                                          int *dist )
{
  boolean ok;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  ok = GeomObjectFindNearestPoint ( obj->me.cpdimen, obj->me.spdimen,
                                    obj->nv, obj->meshvpc,
                                    obj->mkcp, MASK_CP_MOVEABLE, CPos, x, y, dist );
  return ok;
} /*GeomObjectBSplineMeshFindCPoint*/

void GeomObjectBSplineMeshSetCPoint ( GO_BSplineMesh *obj,
                                     CameraRecd *CPos, short x, short y )
{
  int k;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  if ( current_point_ind >= 0 && current_point_ind < obj->nv ) {
    k = obj->me.cpdimen*current_point_ind;
    GeomObjectSetPoint ( CPos, x, y, obj->me.cpdimen, obj->me.spdimen,
                         &obj->meshvpc[k] );
    obj->special_patches_ok = false;
  }
  obj->me.dlistmask = 0;
} /*GeomObjectBSplineMeshSetCPoint*/

void GeomObjectBSplineMeshMarkCPoints ( GO_BSplineMesh *obj,
                                        CameraRecd *CPos, Box2s *box,
                                        byte mask, int action )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  GeomObjectMarkPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv,
                         obj->mkcp, obj->meshvpc, CPos, box, mask, action );
  obj->me.dlistmask &= ~BSM_DLM_CNET;
} /*GeomObjectBSplineMeshMarkCPoints*/

void GeomObjectBSplineMeshMarkCPoint ( GO_BSplineMesh *obj,
                                       byte mask, int action )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  GeomObjectMarkPoint ( obj->nv, obj->mkcp, mask, action );
  obj->me.dlistmask &= ~BSM_DLM_CNET;
} /*GeomObjectBSplineMeshMarkCPoint*/

static boolean _EdgeBoxCoincidence ( CameraRecd *CPos, Box2s *box,
                                     int cpdimen, int spdimen,
                                     double *v0, double *v1 )
{
  point3d p0, p1, q0, q1;
  Box2d   bb;

  switch ( spdimen ) {
case 2:
    switch ( cpdimen ) {
  case 2:
      SetPoint3d ( &p0, v0[0], v0[1], 0.0 );
      SetPoint3d ( &p1, v1[0], v1[1], 0.0 );
      break;
  case 3:
      SetPoint3d ( &p0, v0[0]/v0[2], v0[1]/v0[2], 0.0 );
      SetPoint3d ( &p1, v1[0]/v1[2], v1[1]/v1[2], 0.0 );
      break;
  default:
      return false;
    }
    break;
case 3:
    switch ( cpdimen ) {
  case 3:
      memcpy ( &p0, v0, sizeof(point3d) );
      memcpy ( &p1, v1, sizeof(point3d) );
      break;
  case 4:
      Point4to3d ( (point4d*)v0, &p0 );
      Point4to3d ( (point4d*)v1, &p1 );
      break;
  default:
      return false;
    }
    break;
default:
    return false;
  }
  if ( !CameraClipLine3d ( CPos, &p0, 0.0, &p1, 1.0, &q0, &q1 ) )
    return false;
  bb.x0 = box->x0;  bb.x1 = box->x1;
  bb.y0 = box->y0;  bb.y1 = box->y1;
  return LiangBarskyClip2d ( (point2d*)&q0, (point2d*)&q1, 0.0, 1.0,
                             &bb, NULL, NULL ) != PKGEOM_CLIP_NONE;
} /*_EdgeBoxCoincidence*/

void GeomObjectBSplineMeshMarkHalfedges ( GO_BSplineMesh *obj,
                                          CameraRecd *CPos, Box2s *box,
                                          byte mask, int action )
{
  int         nhe, cpdimen, spdimen, i;
  BSMhalfedge *mhe;
  double      *cp;
  byte        *mkhe;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  nhe  = obj->nhe;
  mhe  = obj->meshhe;
  cp   = obj->meshvpc;
  mkhe = obj->mkhe;
  cpdimen = obj->me.cpdimen;
  spdimen = obj->me.spdimen;
  for ( i = 0; i < nhe; i++ )
    if ( _EdgeBoxCoincidence ( CPos, box, cpdimen, spdimen,
                               &cp[cpdimen*mhe[i].v0], &cp[cpdimen*mhe[i].v1] ) ) {
      switch ( action ) {
    case MARK_SELECT:
        mkhe[i] |= mask;
        break;
    case MARK_UNSELECT:
        mkhe[i] &= ~mask;
        break;
    case MARK_TGSELECT:
        mkhe[i] ^= mask;
        break;
    default:
        break;
      }
    }

  obj->me.dlistmask &= ~BSM_DLM_CNET;
} /*GeomObjectBSplineMeshMarkHalfedges*/

boolean GeomObjectBSplineMeshSaveCPoints ( GO_BSplineMesh *obj )
{
  int size;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  if ( obj->savedcpoints ) free ( obj->savedcpoints );
  size = obj->nv*obj->me.cpdimen*sizeof(double);
  obj->savedcpoints = malloc ( size );
  if ( !obj->savedcpoints )
    return false;
  memcpy ( obj->savedcpoints, obj->meshvpc, size );
  obj->savedsize = size;
  return true;
} /*GeomObjectBSplineMeshSaveCPoints*/

void GeomObjectBSplineMeshUndoLastTransformation ( GO_BSplineMesh *obj )
{
  int size;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  size = obj->nv*obj->me.cpdimen*sizeof(double);
  if ( obj->savedcpoints && obj->savedsize == size ) {
    memcpy ( obj->meshvpc, obj->savedcpoints, size );
    obj->special_patches_ok = false;
    obj->me.dlistmask = 0;
  }
} /*GeomObjectBSplineMeshUndoLastTransformation*/

void GeomObjectBSplineMeshTransformCPoints ( GO_BSplineMesh *obj,
                                             trans3d *tr, byte mask )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  if ( obj->savedcpoints ) {
    GeomObjectTransformPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv,
                                obj->mkcp, mask,
                                obj->savedcpoints, obj->meshvpc, tr );
    obj->special_patches_ok = false;
    obj->me.dlistmask = 0;
  }
} /*GeomObjectBSplineMeshTransformCPoints*/

boolean GeomObjectBSplineMeshGetPointCoord ( GO_BSplineMesh *obj, int p,          
                                  int *spdimen, int *cpdimen, double **pc )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  if ( p < 0 || p >= obj->nv )
    return false;
  *spdimen = obj->me.spdimen;
  *cpdimen = obj->me.cpdimen;
  *pc = &obj->meshvpc[obj->me.cpdimen*p];
  return true;
} /*GeomObjectBSplineMeshGetPointCoord*/

