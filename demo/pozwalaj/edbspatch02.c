
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


boolean GeomObjectBSplinePatchSetDegreeU ( GO_BSplinePatch *obj, int degu )
{
  void    *sp;
  int     degree_u, lastknot_u, degree_v, lastknot_v;
  int     ku, lknu, ncp, d0, d1, d2, cpdimen;
  double  *knots_u, *cp, *knu, *acp, *aknu;
  byte    *mkcp;
  boolean r;
  
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( obj->bsp_type == BSP_TYPE_BLENDING_G1 ||
       obj->bsp_type == BSP_TYPE_BLENDING_G2 )
    return false;
  if ( degu < 1 || degu > MAX_DEGREE )
    return false;

  sp = pkv_GetScratchMemTop ();
  degree_u = obj->degree_u;
  lastknot_u = obj->lastknot_u;
  knots_u = obj->knots_u;
  degree_v = obj->degree_v;
  lastknot_v = obj->lastknot_v;
  cpdimen = obj->me.cpdimen;

  if ( degu > obj->degree_u ) {
        /* degree elevation */
    ku = mbs_NumKnotIntervalsd ( degree_u, lastknot_u, knots_u );
    for ( d0 = 0; knots_u[degree_u+d0+1] == knots_u[degree_u]; d0++ )
      ;
    for ( d1 = 0; knots_u[lastknot_u-degree_u-d1-1] == knots_u[lastknot_u-degree_u];
          d1++ )
      ;
    lknu = lastknot_u+(ku+1-d0-d1)*(degu-degree_u);
    if ( obj->closed_u ) {
      d2 = mbs_KnotMultiplicityd ( min(lastknot_u-1,degree_u+1), &knots_u[1],
                                   knots_u[degree_u] );
      lknu += d2*(degu-degree_u);
    }
    ncp = (lknu-degu)*(lastknot_v-degree_v);
    knu = malloc ( (lknu+1)*sizeof(double) );
    cp = malloc ( ncp*cpdimen*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( !knu || !cp || !mkcp ) {
      if ( knu ) free ( knu );
      if ( cp ) free ( cp );
      if ( mkcp ) free ( mkcp );
      goto failure;
    }
    if ( obj->closed_u ) {
      if ( !mbs_multiBSDegElevClosedd ( 1, (lastknot_v-degree_v)*cpdimen, degree_u,
                         lastknot_u, knots_u, 0, obj->cpoints, degu-degree_u,
                         &d0, &d1, knu, 0, cp ) )
        goto failure;
    }
    else {
      if ( !mbs_multiBSDegElevd ( 1, (lastknot_v-degree_v)*cpdimen, degree_u, lastknot_u,
                            knots_u, 0, obj->cpoints, degu-degree_u,
                            &d0, &d1, knu, 0, cp, false ) )
        goto failure;
    }

if ( d0 != degu || d1 != lknu ) {
  printf ( "Blad: %d, %d, %d, %d\n", d0, degu, d1, lknu );
  exit ( 1 );
}

  }
  else {
        /* degree reduction - we use auxiliary arrays in the scratch memory */
    aknu = pkv_GetScratchMemd ( lastknot_u+1 );
    acp = pkv_GetScratchMemd ( (lastknot_u-degree_u)*(lastknot_v-degree_v)*cpdimen );
    if ( !aknu || !acp )
      goto failure;
    if ( obj->closed_u ) {
      r = mbs_multiBSDegRedClosedd ( 1, (lastknot_v-degree_v)*cpdimen,
                                 degree_u, lastknot_u, knots_u, 0, obj->cpoints,
                                 degree_u-degu, &d0, &lknu, aknu, 0, acp );
    }
    else {
      r = mbs_multiBSDegRedd ( 1, (lastknot_v-degree_v)*cpdimen,
                           degree_u, lastknot_u, knots_u, 0, obj->cpoints,
                           degree_u-degu, &d0, &lknu, aknu, 0, acp );
    }
    if ( !r )
      goto failure;
    knu = malloc ( (lknu+1)*sizeof(double) );
    ncp = (lknu-degu)*(lastknot_v-degree_v);
    cp = malloc ( ncp*cpdimen*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( !knu || !cp || !mkcp ) {
      if ( knu ) free ( knu );
      if ( cp ) free ( cp );
      if ( mkcp ) free ( mkcp );
      goto failure;
    }
    memcpy ( knu, aknu, (lknu+1)*sizeof(double) );
    memcpy ( cp, acp, ncp*cpdimen*sizeof(double) );
  }
  free ( obj->cpoints );
  free ( obj->knots_u );
  free ( obj->mkcp );
  obj->cpoints = cp;
  obj->knots_u = knu;
  obj->mkcp = mkcp;
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  obj->degree_u = degu;
  obj->lastknot_u = lknu;
  obj->me.dlistmask = 0;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplinePatchSetDegreeU*/

boolean GeomObjectBSplinePatchSetDegreeV ( GO_BSplinePatch *obj, int degv )
{
  void    *sp;
  int     degree_u, lastknot_u, degree_v, lastknot_v;
  int     kv, lknv, ncp, d0, d1, d2, cpdimen;
  double  *knots_v, *cp, *knv, *acp, *aknv;
  byte    *mkcp;
  boolean r;
  
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( obj->bsp_type == BSP_TYPE_BLENDING_G1 ||
       obj->bsp_type == BSP_TYPE_BLENDING_G2 )
    return false;
  if ( degv < 1 || degv > MAX_DEGREE )
    return false;

  sp = pkv_GetScratchMemTop ();
  degree_u = obj->degree_u;
  lastknot_u = obj->lastknot_u;
  degree_v = obj->degree_v;
  lastknot_v = obj->lastknot_v;
  knots_v = obj->knots_v;
  cpdimen = obj->me.cpdimen;

  if ( degv > obj->degree_v ) {
        /* degree elevation */
    kv = mbs_NumKnotIntervalsd ( degree_v, lastknot_v, knots_v );
    for ( d0 = 0; knots_v[degree_v+d0+1] == knots_v[degree_v]; d0++ )
      ;
    for ( d1 = 0; knots_v[lastknot_v-degree_v-d1-1] == knots_v[lastknot_v-degree_v];
          d1++ )
      ;
    lknv = lastknot_v+(kv+1-d0-d1)*(degv-degree_v);
    if ( obj->closed_v ) {
      d2 = mbs_KnotMultiplicityd ( min(lastknot_v-1,degree_v+1), &knots_v[1],
                                   knots_v[degree_v] );
      lknv += d2*(degv-degree_v);
    }
    ncp = (lastknot_u-degree_u)*(lknv-degv);
    knv = malloc ( (lknv+1)*sizeof(double) );
    cp = malloc ( ncp*cpdimen*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( !knv || !cp || !mkcp ) {
      if ( knv ) free ( knv );
      if ( cp ) free ( cp );
      if ( mkcp ) free ( mkcp );
      goto failure;
    }
    if ( obj->closed_v ) {
      if ( !mbs_multiBSDegElevClosedd ( lastknot_u-degree_u, cpdimen, degree_v,
                         lastknot_v, knots_v, (lastknot_v-degree_v)*cpdimen,
                         obj->cpoints, degv-degree_v,
                         &d0, &d1, knv, (lknv-degv)*cpdimen, cp ) )
        goto failure;
    }
    else {
      if ( !mbs_multiBSDegElevd ( lastknot_u-degree_u, cpdimen, degree_v, lastknot_v,
                            knots_v, (lastknot_v-degree_v)*cpdimen,
                            obj->cpoints, degv-degree_v,
                            &d0, &d1, knv, (lknv-degv)*cpdimen, cp, false ) )
        goto failure;
    }

if ( d0 != degv || d1 != lknv ) {
  printf ( "Blad: %d, %d, %d, %d\n", d0, degv, d1, lknv );
  exit ( 1 );
}

  }
  else {
        /* degree reduction - we use auxiliary arrays in the scratch memory */
    aknv = pkv_GetScratchMemd ( lastknot_v+1 );
    acp = pkv_GetScratchMemd ( (lastknot_u-degree_u)*(lastknot_v-degree_v)*cpdimen );
    if ( !aknv || !acp )
      goto failure;
    if ( obj->closed_v ) {
      r = mbs_multiBSDegRedClosedd ( lastknot_u-degree_u, cpdimen,
                                 degree_v, lastknot_v, knots_v,
                                 (lastknot_v-degree_v)*cpdimen, obj->cpoints,
                                 degree_v-degv, &d0, &lknv, aknv,
                                 (lastknot_v-degree_v)*cpdimen, acp );
    }
    else {
      r = mbs_multiBSDegRedd ( lastknot_u-degree_u, cpdimen,
                           degree_v, lastknot_v, knots_v,
                           (lastknot_v-degree_v)*cpdimen, obj->cpoints,
                           degree_v-degv, &d0, &lknv, aknv,
                           (lastknot_v-degree_v)*cpdimen, acp );
    }
    if ( !r )
      goto failure;
    knv = malloc ( (lknv+1)*sizeof(double) );
    ncp = (lastknot_u-degree_u)*(lknv-degv);
    cp = malloc ( ncp*cpdimen*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( !knv || !cp || !mkcp ) {
      if ( knv ) free ( knv );
      if ( cp ) free ( cp );
      if ( mkcp ) free ( mkcp );
      goto failure;
    }
    memcpy ( knv, aknv, (lknv+1)*sizeof(double) );
    pkv_Selectd ( lastknot_u-degree_u, (lknv-degv)*cpdimen,
                  (lastknot_v-degree_v)*cpdimen, (lknv-degv)*cpdimen,
                  acp, cp );
  }
  free ( obj->cpoints );
  free ( obj->knots_v );
  free ( obj->mkcp );
  obj->cpoints = cp;
  obj->knots_v = knv;
  obj->mkcp = mkcp;
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  obj->degree_v = degv;
  obj->lastknot_v = lknv;
  obj->me.dlistmask = 0;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplinePatchSetDegreeV*/

