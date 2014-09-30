
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


boolean GeomObjectInitBSplinePatch ( GO_BSplinePatch *obj,
                                     char spdimen, boolean rational )
{
  static const double inikn[4] = {0.0,0.0,1.0,1.0};
  static const double inicp[16] =
    {-1.0,-1.0,0.0,1.0, -1.0,1.0,0.0,1.0,
      1.0,-1.0,0.0,1.0,  1.0,1.0,0.0,1.0};

  obj->me.obj_type = GO_BSPLINE_PATCH;
  obj->me.ident = -1;
  obj->me.spdimen = spdimen;
  obj->rational = rational;
  if ( rational )
    obj->me.cpdimen = spdimen+1;
  else
    obj->me.cpdimen = spdimen;
  obj->me.active = false;
  obj->me.name[0] = 0;
  obj->maxknots_u = obj->maxknots_v = 4;
  obj->savedsize = 0;
  obj->knots_u = malloc ( 4*sizeof(double) );
  obj->knots_v = malloc ( 4*sizeof(double) );
  obj->cpoints = malloc ( 4*obj->me.cpdimen*sizeof(double) );
  obj->mkcp = malloc ( 4 );
  if ( !obj->knots_u || !obj->knots_v || !obj->cpoints || !obj->mkcp ) {
    if ( obj->knots_u ) free ( obj->knots_u );
    if ( obj->knots_v ) free ( obj->knots_v );
    if ( obj->cpoints ) free ( obj->cpoints );
    if ( obj->mkcp ) free ( obj->mkcp );
    return false;
  }
  memcpy ( obj->knots_u, inikn, 4*sizeof(double) );
  memcpy ( obj->knots_v, inikn, 4*sizeof(double) );
  GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                             4, inicp, obj->cpoints );
  memset ( obj->mkcp, MASK_CP_MOVEABLE, 4 );
  obj->lastknot_u = obj->lastknot_v = 3;
  obj->degree_u = obj->degree_v = 1;
  obj->closed_u = obj->closed_v = obj->clamped = false;
  obj->dens_u = obj->dens_v = 6;
  obj->view_surf = obj->view_cnet = true;
  obj->me.displaylist = glGenLists ( BSP_NDL );
  obj->bsp_type = BSP_TYPE_GENERAL;
  memset ( obj->blp_range, 0, 4*sizeof(int) );
  obj->blp_C = 0.001;
  obj->nkn1 = 6;
  obj->nkn2 = 8;
  obj->maxit = 20;
  obj->me.colour[0] = obj->me.colour[1] = 1.0;
  obj->me.colour[2] = 0.0;
  obj->me.display_pretrans = false;
  IdentTrans3d ( &obj->me.pretrans );
  return true;
} /*GeomObjectInitBSplinePatch*/

geom_object *GeomObjectAddBSplinePatch ( const char *name,
                                         char spdimen, boolean rational )
{
  GO_BSplinePatch *obj;

  obj = malloc ( sizeof(GO_BSplinePatch) );
  if ( obj ) {
    memset ( obj, 0, sizeof(GO_BSplinePatch) );
    if ( !GeomObjectInitBSplinePatch ( obj, spdimen, rational ) ) {
      free ( obj );
      return NULL;
    }
    strncpy ( obj->me.name, name, MAX_NAME_LENGTH+1 );
    if ( !first_go )
      first_go = last_go = &obj->me;
    else {
      obj->me.prev = last_go;
      last_go->next = &obj->me;
      last_go = &obj->me;
    }
    current_go = &obj->me;
    return &obj->me;
  }
  else
    return NULL;
} /*GeomObjectAddBSplinePatch*/

geom_object *GeomObjectCopyBSplinePatch ( GO_BSplinePatch *obj )
{
  GO_BSplinePatch *copy;
  int             ncp;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return NULL;
  copy = malloc ( sizeof(GO_BSplinePatch) );
  if ( copy ) {
    memset ( copy, 0, sizeof(GO_BSplinePatch) );
    ncp = (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v);
    copy->me.obj_type = GO_BSPLINE_PATCH;
    copy->me.ident = -1;
    copy->me.spdimen = obj->me.spdimen;
    copy->me.cpdimen = obj->me.cpdimen;
    copy->me.active = false;
    copy->me.dependencies = NULL;
    copy->me.maxdn = 0;
    strcpy ( copy->me.name, obj->me.name );
    copy->me.display_pretrans = false;
    copy->me.pretrans = obj->me.pretrans;
    copy->maxknots_u = obj->maxknots_u;
    copy->maxknots_v = obj->maxknots_v;
    copy->knots_u = malloc ( obj->maxknots_u*sizeof(double) );
    copy->knots_v = malloc ( obj->maxknots_v*sizeof(double) );
    copy->cpoints = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
    copy->mkcp = malloc ( ncp );
    if ( !copy->knots_u || !copy->knots_v || !copy->cpoints || !copy->mkcp ) {
      if ( copy->knots_u ) free ( copy->knots_u );
      if ( copy->knots_v ) free ( copy->knots_v );
      if ( copy->cpoints ) free ( copy->cpoints );
      if ( copy->mkcp ) free ( copy->mkcp );
      free ( copy );
      return NULL;
    }
    copy->degree_u = obj->degree_u;
    copy->degree_v = obj->degree_v;
    copy->lastknot_u = obj->lastknot_u;
    copy->lastknot_v = obj->lastknot_v;
    memcpy ( copy->knots_u, obj->knots_u, (obj->lastknot_u+1)*sizeof(double) );
    memcpy ( copy->knots_v, obj->knots_v, (obj->lastknot_v+1)*sizeof(double) );
    memcpy ( copy->cpoints, obj->cpoints, ncp*obj->me.cpdimen*sizeof(double) );
    memset ( copy->mkcp, MASK_CP_MOVEABLE, ncp );
    copy->rational = obj->rational;
    copy->closed_u = obj->closed_u;
    copy->closed_v = obj->closed_v;
    copy->dens_u = obj->dens_u;
    copy->dens_v = obj->dens_v;
    copy->view_surf = copy->view_cnet = true;
    copy->me.displaylist = glGenLists ( BSP_NDL );
    copy->bsp_type = BSP_TYPE_GENERAL;
    copy->blp_C = obj->blp_C;
    memcpy ( copy->blp_range, obj->blp_range, 4*sizeof(int) );
    copy->clamped = obj->clamped;
    copy->nkn1 = obj->nkn1;
    copy->nkn2 = obj->nkn2;
    copy->maxit = obj->maxit;
    GeomObjectBSplinePatchCopyTrimmedDomain ( obj, copy );
    return &copy->me;
  }
  else
    return NULL;
} /*GeomObjectCopyBSplinePatch*/

void GeomObjectDeleteBSplinePatch ( GO_BSplinePatch *obj )
{
  GeomObjectBSplinePatchDeleteTrimmedDomain ( obj );
  GeomObjectBSplinePatchDestroyBlMat ( obj );
  if ( obj->knots_u )      free ( obj->knots_u );
  if ( obj->knots_v )      free ( obj->knots_v );
  if ( obj->cpoints )      free ( obj->cpoints );
  if ( obj->mkcp )         free ( obj->mkcp );
  if ( obj->savedcpoints ) free ( obj->savedcpoints );
  free ( obj );
} /*GeomObjectDeleteBSplinePatch*/

void GeomObjectAssignBSPatch ( GO_BSplinePatch *obj, int spdimen, boolean rational,
                               int degree_u, int lastknot_u, double *knots_u,
                               int degree_v, int lastknot_v, double *knots_v,
                               double *cpoints, byte *mkcp,
                               boolean closed_u, boolean closed_v )
{
  GeomObjectBSplinePatchDestroyBlMat ( obj );
  if ( obj->cpoints ) free ( obj->cpoints );
  if ( obj->knots_u ) free ( obj->knots_u );
  if ( obj->knots_v ) free ( obj->knots_v );
  if ( obj->mkcp ) free ( obj->mkcp );
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  obj->savedsize = 0;
  obj->me.spdimen = spdimen;
  obj->rational = rational;
  if ( rational )
    obj->me.cpdimen = spdimen+1;
  else
    obj->me.cpdimen = spdimen;
  obj->cpoints = cpoints;
  obj->knots_u = knots_u;
  obj->knots_v = knots_v;
  obj->mkcp = mkcp;
  obj->degree_u = degree_u;
  obj->degree_v = degree_v;
  obj->lastknot_u = lastknot_u;
  obj->lastknot_v = lastknot_v;
  obj->maxknots_u = lastknot_u+1;
  obj->maxknots_v = lastknot_v+1;
  obj->closed_u = closed_u;
  obj->closed_v = closed_v;
  obj->me.dlistmask = 0;
} /*GeomObjectAssignBSPatch*/

boolean GeomObjectBSplinePatchSetRational ( GO_BSplinePatch *obj, boolean rational )
{
  int    ncp, dim, cdim, i, j, k, l;
  double *cp, *rcp, w;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  ncp = (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v);
  dim = obj->me.spdimen;
  cdim = dim+1;
  if ( rational ) {
    if ( obj->rational )
      return true;
    cp = malloc ( ncp*cdim*sizeof(double) );
    if ( !cp )
      return false;
    pkv_Selectd ( ncp, dim, dim, cdim, obj->cpoints, cp );
    for ( i = 0, j = dim;  i < ncp;  i++, j += cdim )
      cp[j] = 1.0;
    obj->me.cpdimen = cdim;
    obj->rational = true;
  }
  else {
    if ( !obj->rational )
      return true;
    cp = malloc ( ncp*dim*sizeof(double) );
    if ( !cp )
      return false;
    rcp = obj->cpoints;
    for ( i = j = k = 0;  i < ncp;  i++, j += cdim, k += dim ) {
      w = rcp[j+dim];
      for ( l = 0; l < dim; l++ )
        cp[k+l] = rcp[j+l]/w;
    }
    obj->me.cpdimen = dim;
    obj->rational = false;
  }
  free ( obj->cpoints );
  obj->cpoints = cp;
  if ( obj->savedcpoints )
    { free ( obj->savedcpoints );  obj->savedcpoints = NULL; }
  obj->savedsize = 0;
  obj->me.dlistmask = 0;
  return true;
} /*GeomObjectBSplinePatchSetRational*/

