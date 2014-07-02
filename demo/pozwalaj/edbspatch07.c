
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
#include "mengerc.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "render.h"
#include "editor.h"
#include "editor_bsp.h"
#include "editor_bsc.h"

/* ////////////////////////////////////////////////////////////////////////// */
void GeomObjectBSplinePatchAdjustGeneral ( GO_BSplinePatch *obj )
{
  if ( obj->bsp_type != BSP_TYPE_GENERAL ) {
    GeomObjectDeleteDependencies ( &obj->me );
    obj->bsp_type = BSP_TYPE_GENERAL;
  }
} /*GeomObjectBSplinePatchAdjustGeneral*/

/* ////////////////////////////////////////////////////////////////////////// */
GO_BSplineCurve *GeomObjectFindBSCurve2DByName ( char *name )
{
  geom_object *obj;

  for ( obj = first_go;  obj;  obj = obj->next ) {
    obj = GeomObjectFindByName ( obj, name );
    if ( !obj )
      return NULL;
    if ( obj->obj_type == GO_BSPLINE_CURVE && obj->spdimen == 2 )
      return (GO_BSplineCurve*)obj;
  }
  return NULL;
} /*GeomObjectFindBSCurve2DByName*/

boolean GeomObjectBSplinePatchGenSphericalProduct ( GO_BSplinePatch *obj )
{
  void            *sp;
  GO_BSplineCurve *equator, *meridian;
  int             degu, degv, lknu, lknv, ncpe, ncpm, ncp, dim, pitch;
  double          *knotsu, *knotsv, *cpoints, *cpe, *cpm;
  byte            *mkcp;
  int             i;

  sp = pkv_GetScratchMemTop ();
  knotsu = knotsv = cpoints = cpe = cpm = NULL;
  mkcp = NULL;
  if ( obj->me.obj_type != GO_BSPLINE_PATCH ||
       obj->me.maxdn < 2 || !obj->me.dependencies )
    goto failure;
  if ( obj->bsp_type != BSP_TYPE_SPHERICAL ||
       !obj->me.dependencies[0] || !obj->me.dependencies[1] )
    goto failure;
  equator = (GO_BSplineCurve*)obj->me.dependencies[0];
  meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
  if ( equator->me.obj_type != GO_BSPLINE_CURVE ||
       equator->me.spdimen != 2 ||
       meridian->me.obj_type != GO_BSPLINE_CURVE ||
       meridian->me.spdimen != 2 )
    goto failure;
  degu = equator->degree;
  lknu = equator->lastknot;
  degv = meridian->degree;
  lknv = meridian->lastknot;
  ncpe = lknu-degu;
  ncpm = lknv-degv;
  ncp = ncpe*ncpm;
        /* the surface is rational if any of the two curves is rational */
  dim = (equator->me.cpdimen > 2 || meridian->me.cpdimen > 2) ? 4 : 3;
  pitch = ncpm*dim;
        /* reallocate the arrays if necessary */
  if ( lknu != obj->lastknot_u ) {
    knotsu = malloc ( (lknu+1)*sizeof(double) );
    if ( !knotsu ) goto failure;
  }
  if ( lknv != obj->lastknot_v ) {
    knotsv = malloc ( (lknv+1)*sizeof(double) );
    if ( !knotsv ) goto failure;
  }
  if ( lknu != obj->lastknot_u || degu != obj->degree_u ||
       lknv != obj->lastknot_v || degv != obj->degree_v ||
       dim != obj->me.cpdimen ) {
    cpoints = malloc ( ncp*dim*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( !cpoints || !mkcp ) goto failure;
    memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  }
  if ( !equator->rational && dim == 4 ) {
    cpe = pkv_GetScratchMemd ( ncpe*3 );
    if ( !cpe )
      goto failure;
        /* convert to homogeneous form */
    pkv_Selectd ( ncpe, 2, 2, 3, equator->cpoints, cpe );
    for ( i = 0; i < ncpe; i++ )
      cpe[3*i+2] = 1.0;
  }
  else
    cpe = equator->cpoints;
  if ( !meridian->rational && dim == 4 ) {
    cpm = pkv_GetScratchMemd ( ncpm*3 );
    if ( !cpm )
      goto failure;
        /* convert to homogeneous form */
    pkv_Selectd ( ncpm, 2, 2, 3, meridian->cpoints, cpm );
    for ( i = 0; i < ncpm; i++ )
      cpm[3*i+2] = 1.0;
  }
  else
    cpm = meridian->cpoints;
        /* assign the new arrays */
  if ( knotsu ) {
    free ( obj->knots_u );
    obj->knots_u = knotsu;
  }
  if ( knotsv ) {
    free ( obj->knots_v );
    obj->knots_v = knotsv;
  }
  if ( cpoints ) {
    free ( obj->cpoints );
    obj->cpoints = cpoints;
  }
  if ( mkcp ) {
    free ( obj->mkcp );
    obj->mkcp = mkcp;
  }
        /* assign the new numbers of knots, degrees etc., copy the knots */
  obj->lastknot_u = lknu;
  obj->maxknots_u = lknu+1;
  obj->degree_u = degu;
  obj->closed_u = equator->closed;
  memcpy ( obj->knots_u, equator->knots, (lknu+1)*sizeof(double) );
  obj->lastknot_v = lknv;
  obj->maxknots_v = lknv+1;
  obj->degree_v = degv;
  memcpy ( obj->knots_v, meridian->knots, (lknv+1)*sizeof(double) );
  obj->closed_v = meridian->closed;
        /* compute the control points */
  if ( equator->rational || meridian->rational ) {
    obj->rational = true;
    obj->me.cpdimen = 4;
    mbs_SphericalProductRd ( degu, lknu, (point3d*)cpe, degv, lknv, (point3d*)cpm,
                             pitch, (point4d*)obj->cpoints );
  }
  else {
    obj->rational = false;
    obj->me.cpdimen = 3;
    mbs_SphericalProductd ( degu, lknu, (point2d*)cpe, degv, lknv, (point2d*)cpm,
                            pitch, (point3d*)obj->cpoints );
  }
  obj->me.spdimen = 3;
  obj->me.dlistmask = 0;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( knotsu ) free ( knotsu );
  if ( knotsv ) free ( knotsv );
  if ( cpoints ) free ( cpoints );
  if ( mkcp ) free ( mkcp );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplinePatchGenSphericalProduct*/

static void _GeomObjectBSplinePatchModifyMeridian ( geom_object *meridian )
{
  double *mcp;

  mcp = ((GO_BSplineCurve*)meridian)->cpoints;
  mcp[0] = 1.0;
  mcp[1] = 1.0;
  if ( ((GO_BSplineCurve*)meridian)->rational ) {
    mcp[3] = 1.0;
    mcp[4] = -1.0;
    mcp[2] = mcp[5] = 1.0;
  }
  else {
    mcp[2] =  1.0;
    mcp[3] = -1.0;
  }
} /*_GeomObjectBSplinePatchModifyMeridian*/

boolean GeomObjectBSplinePatchAdjustSProduct ( GO_BSplinePatch *obj )
{
  void        *sp;
  geom_object *equator, *meridian, *current;

  sp = pkv_GetScratchMemTop ();
  current = current_go;
  if ( obj->bsp_type != BSP_TYPE_SPHERICAL ) {
        /* curve names ought to be nonempty */
    if ( !obj->eqname[0] ) strcpy ( obj->eqname, "eq" );
    if ( !obj->mername[0] ) strcpy ( obj->mername, "mer" );
    equator = (geom_object*)GeomObjectFindBSCurve2DByName ( obj->eqname );
    if ( !equator ) {
      equator = GeomObjectAddBSplineCurve ( obj->eqname, 2, obj->rational );
      if ( !equator )
        goto failure;
    }
    meridian = (geom_object*)GeomObjectFindBSCurve2DByName ( obj->mername );
    if ( !meridian ) {
      meridian = GeomObjectAddBSplineCurve ( obj->mername, 2, obj->rational );
      if ( !meridian )
        goto failure;
      _GeomObjectBSplinePatchModifyMeridian ( meridian );
    }
    if ( obj->me.maxdn != 2 )
      GeomObjectDeleteDependencies ( (geom_object*)obj );
    if ( !GeomObjectAddDependency ( (geom_object*)obj, 2, 0, equator ) )
      goto failure;
    if ( !GeomObjectAddDependency ( (geom_object*)obj, 2, 1, meridian ) )
      goto failure;
    obj->bsp_type = BSP_TYPE_SPHERICAL;
    if ( !GeomObjectBSplinePatchGenSphericalProduct ( obj ) )
      goto failure;
  }
  current_go = current;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  GeomObjectBSplinePatchAdjustGeneral ( obj );
  current_go = current;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplinePatchAdjustSProduct*/

