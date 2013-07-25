
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
#include "editor_bsp.h"


void GeomObjectBSplinePatchFindBBox ( GO_BSplinePatch *obj, Box3d *box )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  GeomObjectFindBBox ( obj->me.cpdimen, obj->rational,
        (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v),
        obj->cpoints, box );
} /*GeomObjectBSplinePatchFindBBox*/

boolean GeomObjectBSplinePatchFindCPoint ( GO_BSplinePatch *obj,
                                           CameraRecd *CPos, short x, short y,
                                           int *dist )
{
  boolean ok;
  int     cn, ncp, i, j, clcKu, clcKv;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  cn = obj->lastknot_v-obj->degree_v;
  ncp = (obj->lastknot_u-obj->degree_u)*cn;
  ok = GeomObjectFindNearestPoint ( obj->me.cpdimen, obj->me.spdimen, ncp,
             obj->cpoints, obj->mkcp, MASK_CP_MOVEABLE, CPos, x, y, dist );
  if ( ok && (obj->closed_u || obj->closed_v) ) {
    clcKu = obj->lastknot_u-2*obj->degree_u;
    clcKv = obj->lastknot_v-2*obj->degree_v;
    i = current_point_ind / cn;
    j = current_point_ind % cn;
    if ( obj->closed_u && i >= clcKu ) i -= clcKu;
    if ( obj->closed_v && j >= clcKv ) j -= clcKv;
    current_point_ind = i*cn+j;
  }
  return ok;
} /*GeomObjectBSplinePatchFindCPoint*/

void GeomObjectBSplinePatchSetCPoint ( GO_BSplinePatch *obj,
                                       CameraRecd *CPos, short x, short y )
{
  int lknu, lknv, dim, degu, degv, ncp, k, l, cn, clcKu, clcKv, i, j;

  if (obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  switch ( obj->bsp_type ) {
case BSP_TYPE_SWEPT:
case BSP_TYPE_SPHERICAL:
case BSP_TYPE_LOFTED:
    return;
default:
    break;
  }
  dim = obj->me.cpdimen;
  lknu = obj->lastknot_u;
  lknv = obj->lastknot_v;
  degu = obj->degree_u;
  degv = obj->degree_v;
  cn = lknv-degv;
  ncp = (lknu-degu)*cn;
  if ( current_point_ind >= 0 && current_point_ind < ncp ) {
    k = dim*current_point_ind;
    GeomObjectSetPoint ( CPos, x, y, dim, obj->me.spdimen,
                         &obj->cpoints[k] );
    i = current_point_ind / cn;
    j = current_point_ind % cn;
    if ( obj->closed_u && i < degu ) {
      clcKu = lknu-2*degu;
      l = k + dim*clcKu*cn;
      memcpy ( &obj->cpoints[l], &obj->cpoints[k], dim*sizeof(double) );
      if ( obj->closed_v && j < degv ) {
        clcKv = lknv-2*degv;
        l = k + dim*clcKv;
        memcpy ( &obj->cpoints[l], &obj->cpoints[k], dim*sizeof(double) );
        l = l + dim*clcKu*cn;
        memcpy ( &obj->cpoints[l], &obj->cpoints[k], dim*sizeof(double) );
      }
    }
    else if ( obj->closed_v && j < degv ) {
      clcKv = lknv-2*degv;
      l = k + dim*clcKv;
      memcpy ( &obj->cpoints[l], &obj->cpoints[k], dim*sizeof(double) );
    }
    switch ( obj->bsp_type ) {
  case BSP_TYPE_BLENDING_G1:
  case BSP_TYPE_BLENDING_G2:
      if ( obj->nharmonic )
        obj->nharmonic = GeomObjectBSplinePatchUpdNHarmonic ( obj );
      break;
  default:
      break;
    }
  }
  obj->me.dlistmask = 0;
} /*GeomObjectBSplinePatchSetCPoint*/

void GeomObjectBSplinePatchMarkCPoints ( GO_BSplinePatch *obj,
                                         CameraRecd *CPos, Box2s *box,
                                         byte mask, boolean clear )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  GeomObjectMarkPoints ( obj->me.cpdimen, obj->me.spdimen,
            (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v),
            obj->mkcp, obj->cpoints, CPos, box, mask, clear );
  obj->me.dlistmask &= ~BSP_DLM_CNET;
} /*GeomObjectBSplinePatchMarkCPoints*/

void GeomObjectBSplinePatchMarkCPoint ( GO_BSplinePatch *obj,
                                        byte mask, boolean clear )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  GeomObjectMarkPoint (
      (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v),
      obj->mkcp, mask, clear );
} /*GeomObjectBSplinePatchMarkCPoint*/

boolean GeomObjectBSplinePatchSaveCPoints ( GO_BSplinePatch *obj )
{
  int ncp, size;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( obj->savedcpoints ) free ( obj->savedcpoints );
  ncp = (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v);
  size = ncp*obj->me.cpdimen*sizeof(double);
  obj->savedcpoints = malloc ( size );
  if ( !obj->savedcpoints )
    return false;
  memcpy ( obj->savedcpoints, obj->cpoints, size );
  obj->savedsize = size;
  return true;
} /*GeomObjectBSplinePatchSaveCPoints*/

void GeomObjectBSplinePatchUndoLastTransformation ( GO_BSplinePatch *obj )
{
  int ncp, size;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  ncp = (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v);
  size = ncp*obj->me.cpdimen*sizeof(double);
  if ( obj->savedcpoints && obj->savedsize == size ) {
    memcpy ( obj->cpoints, obj->savedcpoints, size );
    obj->me.dlistmask = 0;
  }
} /*GeomObjectBSplinePatchUndoLastTransformation*/

void GeomObjectBSplinePatchTransformCPoints ( GO_BSplinePatch *obj,
                                              trans3d *tr, byte mask )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  switch ( obj->bsp_type ) {
case BSP_TYPE_SWEPT:
case BSP_TYPE_SPHERICAL:
case BSP_TYPE_LOFTED:
    return;
default:
    break;
  }
  if ( obj->savedcpoints ) {
    GeomObjectTransformPoints ( obj->me.cpdimen, obj->me.spdimen,
        (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v),
        obj->mkcp, mask, obj->savedcpoints, obj->cpoints, tr );
    switch ( obj->bsp_type ) {
  case BSP_TYPE_BLENDING_G1:
  case BSP_TYPE_BLENDING_G2:
      if ( obj->nharmonic )
        obj->nharmonic = GeomObjectBSplinePatchUpdNHarmonic ( obj );
      break;
  default:
      break;
    }
    obj->me.dlistmask = 0;
  }
} /*GeomObjectBSplinePatchTransformCPoints*/

boolean GeomObjectBSplinePatchGetPointCoord ( GO_BSplinePatch *obj, int p,
                                  int *spdimen, int *cpdimen, double **pc )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( p < 0 || p >= (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v) )
    return false;
  *spdimen = obj->me.spdimen;
  *cpdimen = obj->me.cpdimen;
  *pc = &obj->cpoints[obj->me.cpdimen*p];
  return true;
} /*GeomObjectBSplinePatchGetPointCoord*/

