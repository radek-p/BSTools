
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
#include "editor_bsh.h"


boolean GeomObjectInitBSplineHole ( GO_BSplineHole *obj, boolean rational )
{
  return false;
} /*GeomObjectInitBSplineHole*/

geom_object *GeomObjectAddBSplineHole ( const char *name, boolean rational )
{
  return NULL;
} /*GeomObjectAddBSplineHole*/

geom_object *GeomObjectCopyBSplineHole ( GO_BSplineHole *obj )
{
  return NULL;
} /*GeomObjectCopyBSplineHole*/

void GeomObjectDeleteBSplineHole ( GO_BSplineHole *obj )
{
  free ( obj );
} /*GeomObjectDeleteBSplineHole*/

void GeomObjectDrawBSplineHole ( GO_BSplineHole *obj )
{
} /*GeomObjectDrawBSplineHole*/

void GeomObjectDrawBSplineHNet ( GO_BSplineHole *obj )
{
} /*GeomObjectDrawBSplineHNet*/

void GeomObjectDisplayBSplineHole ( GO_BSplineHole *obj )
{
} /*GeomObjectDisplayBSplineHole*/

void GeomObjectBSplineHoleFindBBox ( GO_BSplineHole *obj, Box3d *box )
{
} /*GeomObjectBSplineHoleFindBBox*/

boolean GeomObjectBSplineHoleFindCPoint ( GO_BSplineHole *obj,
                                          CameraRecd *CPos, short x, short y,
                                          int *dist )
{
  return false;
} /*GeomObjectBSplineHoleFindCPoint*/

void GeomObjectBSplineHoleSetCPoint ( GO_BSplineHole *obj,
                                      CameraRecd *CPos, short x, short y )
{
} /*GeomObjectBSplineHoleSetCPoint*/

void GeomObjectBSplineHoleMarkCPoints ( GO_BSplineHole *obj,
                                        CameraRecd *CPos, Box2s *box,
                                        byte mask, boolean clear )
{
  if ( obj->me.obj_type != GO_BSPLINE_HOLE )
    return;
/* ********* */
} /*GeomObjectBSplineHoleMarkCPoints*/

void GeomObjectBSplineHoleMarkCPoint ( GO_BSplineHole *obj,
                                       byte mask, boolean clear )
{
  if ( obj->me.obj_type != GO_BSPLINE_HOLE )
    return;
  GeomObjectMarkPoint ( 12*obj->hole_k+1, obj->mkcp, mask, clear );
} /*GeomObjectBezierCurveMarkCPoint*/

boolean GeomObjectBSplineHoleSaveCPoints ( GO_BSplineHole *obj )
{
  return false;
} /*GeomObjectBSplineHoleSaveCPoints*/

void GeomObjectBSplineHoleUndoLastTransformation ( GO_BSplineHole *obj )
{
} /*GeomObjectBSplineHoleUndoLastTransformation*/

void GeomObjectBSplineHoleTransformCPoints ( GO_BSplineHole *obj,
                                             trans3d *tr, byte mask )
{
} /*GeomObjectBSplineHoleTransformCPoints*/

boolean GeomObjectBSplineHoleGetPointCoord ( GO_BSplineHole *obj, int p,
                                  int *spdimen, int *cpdimen, double **pc )
{
  return false;
} /*GeomObjectBSplineHoleGetPointCoord*/

boolean GeomObjectWriteBSplineHole ( GO_BSplineHole *obj )
{
  return false;
} /*GeomObjectWriteBSplineHole*/

void GeomObjectReadBSplineHole ( void *usrdata, const char *name, int ident,
                                 int hole_k, const double *knots,
                                 const point2d *domain_cp, const point4d *hole_cp,
                                 int spdimen, boolean rational )
{
} /*GeomObjectReadBSplineHole*/

void GeomObjectBSplineHoleDisplayInfoText ( GO_BSplineHole *obj )
{
  SetStatusText ( "", true );
} /*GeomObjectBSplineHoleDisplayInfoText*/

boolean GeomObjectBSplineHoleProcessDep ( GO_BSplineHole *obj, geom_object *go )
{
  return false;
} /*GeomObjectBSplineHoleProcessDep*/

void GeomObjectBSplineHoleProcessDeletedDep ( GO_BSplineHole *obj, geom_object *go )
{
} /*GeomObjectBSplineHoleProcessDeletedDep*/

