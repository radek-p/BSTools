
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

#include "editor.h"  
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"

#define PARENT_SIDE
#include "pozwalajipc.h"
#undef PARENT_SIDE

geom_object *first_go = NULL, *last_go = NULL, *current_go = NULL,
            *currentp_go = NULL;
int         current_point_ind;
double      *current_point = NULL;
char        current_point_dim;
boolean     current_point_rational;
byte        marking_mask = MASK_CP_MARKED_0;

void GeomObjectInitList ( void )
{
  first_go = last_go = current_go = NULL;
} /*InitGeomObjectList*/

void GeomObjectPurgeList ( void )
{
  while ( first_go ) {
    current_go = first_go;
    GeomObjectDeleteCurrent ();
  }
} /*GeomObjectPurgeList*/

void GeomObjectDeleteCurrent ( void )
{
  geom_object *obj;

  obj = current_go;
  if ( !obj )
    return;

  if ( obj->bound_with_a_child )
    BindChildToGeomObject ( NULL );
  GeomObjectProcessDeletedDep ( obj );
  GeomObjectDeleteDependencies ( obj );
        /* extract the object from the object list */
  if ( obj == first_go ) {
    first_go = obj->next;
    if ( first_go ) {
      first_go->prev = NULL;
      current_go = first_go;
    }
    else
      current_go = last_go = NULL;
  }
  else if ( obj == last_go ) {
    current_go = last_go = obj->prev;
    last_go->next = NULL;
  }
  else {
    current_go = obj->next;
    current_go->prev = obj->prev;
    obj->prev->next = current_go;
  }
        /* deallocate the object */
  switch ( obj->obj_type ) {
case GO_BEZIER_CURVE:
    GeomObjectDeleteBezierCurve ( (GO_BezierCurve*)obj );
    break;
case GO_BSPLINE_CURVE:
    GeomObjectDeleteBSplineCurve ( (GO_BSplineCurve*)obj );
    break;
case GO_BEZIER_PATCH:
    GeomObjectDeleteBezierPatch ( (GO_BezierPatch*)obj );
    break;
case GO_BSPLINE_PATCH:
    GeomObjectDeleteBSplinePatch ( (GO_BSplinePatch*)obj );
    break;
case GO_BSPLINE_MESH:
    GeomObjectDeleteBSplineMesh ( (GO_BSplineMesh*)obj );
    break;
case GO_BSPLINE_HOLE:
    GeomObjectDeleteBSplineHole ( (GO_BSplineHole*)obj );
    break;
default:  /* should never happen */
    break;
  }
  GeomObjectSortDependencies ();
} /*GeomObjectDeleteCurrent*/

boolean GeomObjectCopyCurrent ( void )
{
  geom_object *obj;

  if ( !current_go )
    return false;
  switch ( current_go->obj_type ) {
case GO_BEZIER_CURVE:
    obj = GeomObjectCopyBezierCurve ( (GO_BezierCurve*)current_go );
    break;
case GO_BSPLINE_CURVE:
    obj = GeomObjectCopyBSplineCurve ( (GO_BSplineCurve*)current_go );
    break;
case GO_BEZIER_PATCH:
    obj = GeomObjectCopyBezierPatch ( (GO_BezierPatch*)current_go );
    break;
case GO_BSPLINE_PATCH:
    obj = GeomObjectCopyBSplinePatch ( (GO_BSplinePatch*)current_go );
    break;
case GO_BSPLINE_MESH:
    obj = GeomObjectCopyBSplineMesh ( (GO_BSplineMesh*)current_go );
    break;
case GO_BSPLINE_HOLE:
    obj = GeomObjectCopyBSplineHole ( (GO_BSplineHole*)current_go );
    break;
default:
    return false;
  }
  if ( !obj )
    return false;
  memcpy ( obj->colour, current_go->colour, 3*sizeof(double) );
      /* insert the new object to the list */
  obj->prev = current_go;
  obj->next = current_go->next;
  current_go->next = obj;
  if ( current_go == last_go )
    last_go = obj;
  current_go = obj;
  return true;
} /*GeomObjectCopyCurrent*/

int GeomObjectNumber ( void )
{
  int         cnt;
  geom_object *go;

  for ( go = first_go, cnt = 0;  go;  go = go->next )
    cnt ++;
  return cnt;
} /*GeomObjectNumber*/

int GeomObjectCurrentNumber ( void )
{
  int         cnt;
  geom_object *go;

  for ( go = first_go, cnt = 0;  go;  go = go->next, cnt ++ )
    if ( go == current_go )
      return cnt;
  return -1;
} /*GeomObjectCurrentNumber*/

boolean GeomObjectSelect ( int objno )
{
  int         cnt;
  geom_object *go;

  for ( go = first_go, cnt = 0;  go;  go = go->next, cnt ++ )
    if ( cnt == objno ) {
      current_go = go;
      return true;
    }
  return false;
} /*GeomObjectSelect*/

boolean GeomObjectIsInTheList ( geom_object *obj )
{
  geom_object *go;

  if ( !obj )
    return false;
  for ( go = first_go; go; go = go->next )
    if ( obj == go )
      return true;
  return false;
} /*GeomObjectIsInTheList*/

void GeomObjectDisplayActive ( char spdimen )
{
  geom_object *go;

  for ( go = first_go; go; go = go->next )
    if ( (go == current_go || go->active) && go->spdimen == spdimen ) {
      glPushMatrix ();
      if ( go == current_go && go->display_pretrans )
        xgle_MultMatrix3d ( &go->pretrans );
      switch ( go->obj_type ) {
    case GO_BEZIER_CURVE:
        GeomObjectDisplayBezierCurve ( (GO_BezierCurve*)go );
        break;
    case GO_BSPLINE_CURVE:
        GeomObjectDisplayBSplineCurve ( (GO_BSplineCurve*)go );
        break;
    case GO_BEZIER_PATCH:
        GeomObjectDisplayBezierPatch ( (GO_BezierPatch*)go );
        break;
    case GO_BSPLINE_PATCH:
        GeomObjectDisplayBSplinePatch ( (GO_BSplinePatch*)go );
        break;
    case GO_BSPLINE_MESH:
        GeomObjectDisplayBSplineMesh ( (GO_BSplineMesh*)go );
        break;
    case GO_BSPLINE_HOLE:
        GeomObjectDisplayBSplineHole ( (GO_BSplineHole*)go );
        break;
    default:
        break;
      }
      glPopMatrix ();
    }
} /*GeomObjectDisplayActive*/

boolean GeomObjectFindBoundingBox ( char spdimen, Box3d *box )
{
  boolean     itis;
  Box3d       b;
  double      w;
  geom_object *go;

  itis = false;
  for ( go = first_go; go; go = go->next )
    if ( (go == current_go || go->active) && go->spdimen == spdimen ) {
      switch ( go->obj_type ) {
    case GO_BEZIER_CURVE:
        GeomObjectBezierCurveFindBBox ( (GO_BezierCurve*)go, &b );
        break;
    case GO_BSPLINE_CURVE:
        GeomObjectBSplineCurveFindBBox ( (GO_BSplineCurve*)go, &b );
        break;
    case GO_BEZIER_PATCH:
        GeomObjectBezierPatchFindBBox ( (GO_BezierPatch*)go, &b );
        break;
    case GO_BSPLINE_PATCH:
        GeomObjectBSplinePatchFindBBox ( (GO_BSplinePatch*)go, &b );
        break;
    case GO_BSPLINE_MESH:
        GeomObjectBSplineMeshFindBBox ( (GO_BSplineMesh*)go, &b );
        break;
    case GO_BSPLINE_HOLE:
        GeomObjectBSplineHoleFindBBox ( (GO_BSplineHole*)go, &b );
        break;
    default:
        break;
      }
      if ( itis ) {
        if ( b.x0 < box->x0 ) box->x0 = b.x0;
        if ( b.x1 > box->x1 ) box->x1 = b.x1;
        if ( b.y0 < box->y0 ) box->y0 = b.y0;
        if ( b.y1 > box->y1 ) box->y1 = b.y1;
        if ( b.z0 < box->z0 ) box->z0 = b.z0;
        if ( b.z1 > box->z1 ) box->z1 = b.z1;
      }
      else {
        memcpy ( box, &b, sizeof(Box3d) );
        itis = true;
      }
    }
  if ( itis ) {
        /* extend the box by 5% */
    w = 0.025*(box->x1-box->x0);  box->x0 -= w;  box->x1 += w;
    w = 0.025*(box->y1-box->y0);  box->y0 -= w;  box->y1 += w;
    w = 0.025*(box->z1-box->z0);  box->z0 -= w;  box->z1 += w;
    return true;
  }
  else
    return false;
} /*GeomObjectFindBoundingBox*/

geom_object *GeomObjectFindByName ( geom_object *obj, const char *name )
{
  while ( obj ) {
    if ( !strcmp ( name, obj->name ) )
      break;
    obj = obj->next;
  }
  return obj;
} /*GeomObjectFindByName*/

