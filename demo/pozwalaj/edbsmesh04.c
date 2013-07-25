
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


boolean GeomObjectBSplineMeshSetCurrentVertex ( GO_BSplineMesh *obj, int vn, int cvi )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH ||
       cvi < 0 || cvi >= BSM_NCURRENT_VERTICES )
    return false;
  if ( vn == -2 )
    obj->current_vertex[cvi] = obj->nv-1;
  else if ( vn >= 0 && vn < obj->nv )
    obj->current_vertex[cvi] = vn;
  else
    obj->current_vertex[cvi] = -1;
  obj->me.dlistmask &= ~BSM_DLM_ELEM;
  return true;
} /*GeomObjectBSplineMeshSetCurrentVertex*/

boolean GeomObjectBSplineMeshSetCurrentEdge ( GO_BSplineMesh *obj, int en, int cei )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH ||
       cei < 0 || cei >= BSM_NCURRENT_EDGES )
    return false;
  if ( en == -2 )
    obj->current_edge[cei] = obj->nhe-1;
  else if ( en >= 0 && en < obj->nhe )
    obj->current_edge[cei] = en;
  else
    obj->current_edge[cei] = -1;
  obj->me.dlistmask &= ~BSM_DLM_ELEM;
  return true;
} /*GeomObjectBSplineMeshSetCurrentEdge*/

boolean GeomObjectBSplineMeshSetCurrentFacet ( GO_BSplineMesh *obj, int fn, int cfi )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH ||
       cfi < 0 || cfi >= BSM_NCURRENT_FACETS )
    return false;
  if ( fn == -2 )
    obj->current_facet[cfi] = obj->nfac-1;
  else if ( fn >= 0 && fn < obj->nfac )
    obj->current_facet[cfi] = fn;
  else
    obj->current_facet[cfi] = -1;
  obj->me.dlistmask &= ~BSM_DLM_ELEM;
  return true;
} /*GeomObjectBSplineMeshSetCurrentFacet*/

boolean GeomObjectBSplineMeshSetDegree ( GO_BSplineMesh *obj, int deg )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  if ( deg < 1 || deg > MAX_DEGREE )
    return false;
  obj->degree = deg;
  obj->me.dlistmask &= ~(BSM_DLM_SURF | BSM_DLM_SPECIAL);
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;
} /*GeomObjectBSplineMeshSetDegree*/

boolean GeomObjectBSplineMeshGetCurrentVertices ( GO_BSplineMesh *obj,
                                                  int n, point3d *p )
{
  int i, j;
  double *mvpc;

  if ( obj->me.obj_type != GO_BSPLINE_MESH || n > BSM_NCURRENT_VERTICES )
    return false;
  for ( i = 0; i < n; i++ )
    if ( obj->current_vertex[i] < 0 )
      return false;
  mvpc = obj->meshvpc;
  for ( i = 0; i < n; i++ ) {
    j = obj->current_vertex[i];
    switch ( obj->me.cpdimen ) {
  case 2:
      memcpy ( &p[i], &mvpc[2*j], 2*sizeof(double) );
      p[i].x = 0.0;
      break;
  case 3:
      if ( obj->me.spdimen == 2 )
        SetPoint3d ( &p[i], mvpc[3*j]/mvpc[3*j+2], mvpc[3*j+1]/mvpc[3*j+2], 0.0 );
      else
        memcpy ( &p[i], &mvpc[3*j], 3*sizeof(double) );
      break;
  case 4:
      Point4to3d ( (point4d*)&mvpc[4*j], &p[i] );
      break;
  default:
      return false;
    }
  }
  return true;
} /*GeomObjectBSplineMeshGetCurrentVertices*/

