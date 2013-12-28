
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
               pitch, obj->cpoints, obj->mkcp, obj->me.name,
               (bsf_WriteAttr_fptr)bsf_WriteColour,
               (void*)obj->me.colour );
} /*GeomObjectWriteBSplinePatch*/

void GeomObjectReadBSplinePatch ( void *usrdata,
                 const char *name,
                 int degreeu, int lastknotu, const double *knotsu, 
                 int degreev, int lastknotv, const double *knotsv,
                 boolean closed_u, boolean closed_v,
                 int pitch, const point4d *cpoints,
                 int spdimen, boolean rational, byte *_mkcp )
{
  GO_BSplinePatch *obj;
  double          *cp, *knu, *knv;
  byte            *mkcp;
  int             ncp, cpdimen, i;

  obj = (GO_BSplinePatch*)GeomObjectAddBSplinePatch ( name, spdimen, rational );
  if ( obj ) {
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
    memcpy ( mkcp, _mkcp, ncp );
    for ( i = 0; i < ncp; i++ )
      mkcp[i] |= MASK_CP_MOVEABLE;
  }
} /*GeomObjectReadBSplinePatch*/

