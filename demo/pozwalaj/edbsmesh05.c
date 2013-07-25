
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
                                        byte mask, boolean clear )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  GeomObjectMarkPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv,
                         obj->mkcp, obj->meshvpc, CPos, box, mask, clear );
  obj->me.dlistmask &= ~BSM_DLM_CNET;
} /*GeomObjectBSplineMeshMarkCPoints*/

void GeomObjectBSplineMeshMarkCPoint ( GO_BSplineMesh *obj,
                                       byte mask, boolean clear )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  GeomObjectMarkPoint ( obj->nv, obj->mkcp, mask, clear );
  obj->me.dlistmask &= ~BSM_DLM_CNET;
} /*GeomObjectBSplineMeshMarkCPoint*/

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

