
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
#include "pkrender.h"
#include "xgedit.h"
#include "xgledit.h"

#include "editor.h"
#include "editor_bsp.h"

extern pkRenderer rend;

void GeomObjectBSplinePatchOutputToRenderer3D ( GO_BSplinePatch *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  if ( obj->me.spdimen != 3 )
    return;
  if ( obj->rational )
    RendEnterBSPatch3Rd ( &rend, obj->degree_u, obj->lastknot_u, obj->knots_u,
                          obj->degree_v, obj->lastknot_v, obj->knots_v,
                          (point4d*)obj->cpoints, obj->me.colour );
  else
    RendEnterBSPatch3d ( &rend, obj->degree_u, obj->lastknot_u, obj->knots_u,
                         obj->degree_v, obj->lastknot_v, obj->knots_v,
                         (point3d*)obj->cpoints, obj->me.colour );
} /*GeomObjectBSplinePatchOutputToRenderer3D*/

