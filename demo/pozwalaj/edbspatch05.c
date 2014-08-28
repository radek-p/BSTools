
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


boolean GeomObjectWriteBSPAttributes ( GO_BSplinePatch *obj )
{
  void *sp;
  int  *ident;
  int  i, maxdn;

  sp = pkv_GetScratchMemTop ();
  maxdn = obj->me.maxdn;
  if ( maxdn && obj->me.dependencies ) {  /* there are dependencies to write */
    ident = pkv_GetScratchMemi ( maxdn );
    if ( !ident )
      goto failure;
    for ( i = 0; i < maxdn; i++ )
      if ( obj->me.dependencies[i] )
        ident[i] = obj->me.dependencies[i]->ident;
      else
        ident[i] = -1;
    switch ( obj->bsp_type ) {
  case BSP_TYPE_SPHERICAL:
                     /* write the identifiers of the two curves */
      if ( ident[0] >= 0 && ident[1] >= 0 &&
           obj->me.dependencies[0]->obj_type == GO_BSPLINE_CURVE &&
           obj->me.dependencies[1]->obj_type == GO_BSPLINE_CURVE ) {
        bsf_WriteDependencies ( BSF_DEP_SPHERICAL, 2, ident );
      }
      break;
  default:
      break;
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectWriteBSPAttributes*/

boolean GeomObjectWriteBSplinePatch ( GO_BSplinePatch *obj )
{
  int pitch;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  pitch = (obj->lastknot_v-obj->degree_v)*obj->me.cpdimen;
  return bsf_WriteBSplinePatchd ( obj->me.spdimen, obj->me.cpdimen, obj->rational,
               obj->degree_u, obj->lastknot_u, obj->knots_u,
               obj->degree_v, obj->lastknot_v, obj->knots_v,
               obj->closed_u, obj->closed_v,
               pitch, obj->cpoints, obj->me.name, obj->me.ident,
               GeomObjectWriteAttributes,
               (void*)&obj->me );
} /*GeomObjectWriteBSplinePatch*/

void GeomObjectReadBSplinePatch ( void *usrdata,
                 const char *name, int ident,
                 int degreeu, int lastknotu, const double *knotsu, 
                 int degreev, int lastknotv, const double *knotsv,
                 boolean closed_u, boolean closed_v,
                 int pitch, const point4d *cpoints,
                 int spdimen, boolean rational )
{
  GO_BSplinePatch *obj;
  double          *cp, *knu, *knv;
  byte            *mkcp;
  int             ncp, cpdimen;

  obj = (GO_BSplinePatch*)GeomObjectAddBSplinePatch ( name, spdimen, rational );
  if ( obj ) {
    obj->me.ident = ident;
    ncp = (lastknotu-degreeu)*(lastknotv-degreev);
    cpdimen = obj->me.cpdimen;
    cp = malloc ( ncp*cpdimen*sizeof(double) );
    knu = malloc ( (lastknotu+1)*sizeof(double) );
    knv = malloc ( (lastknotv+1)*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( !cp || !knu || !knv || !mkcp ) {
      if ( mkcp ) free ( mkcp );
      if ( knu ) free ( knu );
      if ( knv ) free ( knv );
      if ( cp ) free ( cp );
      GeomObjectDeleteBSplinePatch ( obj );
      return;
    }
    free ( obj->cpoints );
    free ( obj->knots_u );
    free ( obj->knots_v );
    free ( obj->mkcp );
    if ( obj->savedcpoints ) {
      free ( obj->savedcpoints );
      obj->savedcpoints = NULL;
    }
    obj->cpoints = cp;
    obj->knots_u = knu;
    obj->knots_v = knv;
    obj->mkcp = mkcp;
    obj->degree_u = degreeu;
    obj->degree_v = degreev;
    obj->lastknot_u = lastknotu;
    obj->lastknot_v = lastknotv;
    obj->closed_u = closed_u;
    obj->closed_v = closed_v;
    obj->maxknots_u = lastknotu+1;
    obj->maxknots_v = lastknotv+1;
    GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                               ncp, (double*)cpoints, cp );
    memcpy ( knu, knotsu, (lastknotu+1)*sizeof(double) );
    memcpy ( knv, knotsv, (lastknotv+1)*sizeof(double) );
  }
} /*GeomObjectReadBSplinePatch*/

boolean GeomObjectBSPResolveDependencies ( GO_BSplinePatch *obj )
{
  geom_object *go, *eq, *mer;

  switch ( obj->me.filedepname ) {
case BSF_DEP_SPHERICAL:
    if ( obj->me.filedepnum != 2 ||
         obj->me.filedepid[0] < 0 || obj->me.filedepid[1] < 0 )
      goto failure;
        /* try to find the equator and meridian */
    eq = mer = NULL;
    for ( go = first_go; go; go = go->next ) {
      if ( go->ident == obj->me.filedepid[0] ) {
        if ( go->obj_type == GO_BSPLINE_CURVE ) {
          if ( eq )
            goto failure;
          else
            eq = go;
        }
      }
      else if ( go->ident == obj->me.filedepid[1] ) {
        if ( go->obj_type == GO_BSPLINE_CURVE ) {
          if ( mer )
            goto failure;
          else
            mer = go;
        }
      }
    }
        /* both curves found; setup the spherical product */
    if ( eq && mer ) {
      if ( obj->me.maxdn != 2 )
        GeomObjectDeleteDependencies ( (geom_object*)obj );
      if ( !GeomObjectAddDependency ( (geom_object*)obj, 2, 0, eq ) )
        goto failure;
      if ( !GeomObjectAddDependency ( (geom_object*)obj, 2, 1, mer ) )
        goto failure;
          /* names of the two curves must be nonempty */
      if ( !eq->name[0] ) strcpy ( eq->name, "eq" );
      if ( !mer->name[0] ) strcpy ( mer->name, "mer" );
      strncpy ( obj->eqname, eq->name, MAX_NAME_LENGTH );
      strncpy ( obj->mername, mer->name, MAX_NAME_LENGTH );
      obj->bsp_type = BSP_TYPE_SPHERICAL;
      if ( !GeomObjectBSplinePatchGenSphericalProduct ( obj ) )
        goto failure;
    }
    break;

default:
    break;
  }
  return true;

failure:
  GeomObjectDeleteDependencies ( (geom_object*)obj );
  obj->bsp_type = BSP_TYPE_GENERAL;
  return false;
} /*GeomObjectBSPResolveDependencies*/

