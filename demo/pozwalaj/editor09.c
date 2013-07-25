
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
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"

#define TAG_WHITE 0
#define TAG_GREY  1
#define TAG_BLACK 2

static geom_object *_GOSortDepDFS ( geom_object *go, geom_object *depl )
{
  int i;

  if ( go->dependencies ) {
    go->tag = TAG_GREY;
    for ( i = 0; i < go->maxdn; i++ )
      if ( go->dependencies[i] )
        depl = _GOSortDepDFS ( go->dependencies[i], depl );
  }
  go->tag = TAG_BLACK;
  go->dep_list = depl;
  return go;
} /*_GOSortDepDFS*/

boolean GeomObjectSortDependencies ( void )
{
  geom_object *go, *depl, *g;

  if ( !first_go )
    return true;
  for ( go = first_go; go; go = go->next ) {
    go->dep_list = NULL;
    go->tag = TAG_WHITE;
  }
  depl = NULL;
  for ( go = first_go; go; go = go->next )
    if ( go->tag == TAG_WHITE ) {
      depl = _GOSortDepDFS ( go, depl );
      if ( !depl )
        return false;
    }
        /* reverse the list */
  while ( depl ) {
    g = depl;          depl = depl->dep_list;
    g->dep_list = go;  go = g;
  }
  return true;
} /*GeomObjectSortDependencies*/

static boolean _GOIsDependent ( geom_object *dgo, geom_object *go )
{
  int i;

  if ( go == dgo )
    return true;
  else if ( dgo->dependencies ) {
    for ( i = 0; i < dgo->maxdn; i++ )
      if ( _GOIsDependent ( dgo->dependencies[i], go ) )
        return true;
  }
  return false;
} /*_GOIsDependent*/

boolean GeomObjectAddDependency ( geom_object *go, int maxdn,
                                  int dn, geom_object *dgo )
{
  geom_object **auxdep;

        /* the object pointed by go will be dependent on the one */
        /* pointed by dgo; changes to *dgo will cause recomputing *go */
        /* but the dependency directed graph must be acyclic */

        /* check, if the new dependency will not introduce a cycle */
  if ( _GOIsDependent ( dgo, go ) )
    return false;
        /* allocate or reallocate the object's dependencies array */
  if ( maxdn <= 0 || dn < 0 || dn >= maxdn )
    return false;
  if ( !go->dependencies ) {
    go->dependencies = malloc ( maxdn*sizeof(geom_object*) );
    if ( !go->dependencies )
      return false;
    memset ( go->dependencies, 0, maxdn*sizeof(geom_object*) );
    go->maxdn = maxdn;
  }
  else if ( go->maxdn < maxdn ) {
    auxdep = malloc ( maxdn*sizeof(geom_object*) );
    if ( !auxdep )
      return false;
    memset ( auxdep, 0, maxdn*sizeof(geom_object*) );
    memcpy ( auxdep, go->dependencies, go->maxdn*sizeof(geom_object*) );
    free ( go->dependencies );
    go->dependencies = auxdep;
    go->maxdn = maxdn;
  }
        /* assign the dependency link */
  go->dependencies[dn] = dgo;
        /* topological sort */
  GeomObjectSortDependencies ();
  return true;
} /*GeomObjectAddDependency*/

void GeomObjectRemoveDependency ( geom_object *go, int dn )
{
  if ( go->dependencies && dn >= 0 && dn < go->maxdn )
    go->dependencies[dn] = NULL;
} /*GeomObjectRemoveDependency*/

void GeomObjectDeleteDependencies ( geom_object *go )
{
  if ( go->dependencies ) {
    free ( go->dependencies );
    go->dependencies = NULL;
  }
  go->maxdn = 0;
} /*GeomObjectDeleteDependencies*/

boolean GeomObjectProcessDependencies ( geom_object *go )
{
  geom_object *dgo;
  int         i;
  boolean     result, changed;

  result = false;
  for ( dgo = first_go; dgo; dgo = dgo->next )
    dgo->tag = TAG_WHITE;
  go->tag = TAG_BLACK;
  for ( dgo = go->dep_list; dgo; dgo = dgo->dep_list ) {
    if ( dgo->dependencies )
      for ( i = 0; i < dgo->maxdn; i++ ) {
        if ( dgo->dependencies[i] &&
             dgo->dependencies[i]->tag != TAG_WHITE ) {
          dgo->tag = TAG_GREY;
          switch ( dgo->obj_type ) {
        case GO_BEZIER_CURVE:
            changed = GeomObjectBezierCurveProcessDep ( (GO_BezierCurve*)dgo, go );
            break;
        case GO_BSPLINE_CURVE:
            changed = GeomObjectBSplineCurveProcessDep ( (GO_BSplineCurve*)dgo, go );
            break;
        case GO_BEZIER_PATCH:
            changed = GeomObjectBezierPatchProcessDep ( (GO_BezierPatch*)dgo, go );
            break;
        case GO_BSPLINE_PATCH:
            changed = GeomObjectBSplinePatchProcessDep ( (GO_BSplinePatch*)dgo, go );
            break;
        case GO_BSPLINE_MESH:
            changed = GeomObjectBSplineMeshProcessDep ( (GO_BSplineMesh*)dgo, go );
            break;
        case GO_BSPLINE_HOLE:
            changed = GeomObjectBSplineHoleProcessDep ( (GO_BSplineHole*)dgo, go );
            break;
        default:
            changed = false;
            break;
          }
          dgo->tag = changed ? TAG_BLACK : TAG_WHITE;
          result = true;
        }
      }
  }
  return result;
} /*GeomObjectProcessDependencies*/

void GeomObjectProcessDeletedDep ( geom_object *dgo )
{
  geom_object *go;
  int         i;

  for ( go = first_go;  go;  go = go->next )
    for ( i = 0; i < go->maxdn; i++ )
      if ( go->dependencies[i] == dgo ) {
        switch ( go->obj_type ) {
      case GO_BEZIER_CURVE:
          GeomObjectBezierCurveProcessDeletedDep ( (GO_BezierCurve*)go, dgo );
          break;
      case GO_BEZIER_PATCH:
          GeomObjectBezierPatchProcessDeletedDep ( (GO_BezierPatch*)go, dgo );
          break;
      case GO_BSPLINE_CURVE:
          GeomObjectBSplineCurveProcessDeletedDep ( (GO_BSplineCurve*)go, dgo );
          break;
      case GO_BSPLINE_PATCH:
          GeomObjectBSplinePatchProcessDeletedDep ( (GO_BSplinePatch*)go, dgo );
          break;
      case GO_BSPLINE_MESH:
          GeomObjectBSplineMeshProcessDeletedDep ( (GO_BSplineMesh*)go, dgo );
          break;
      case GO_BSPLINE_HOLE:
          GeomObjectBSplineHoleProcessDeletedDep ( (GO_BSplineHole*)go, dgo );
          break;
      default:
          break;
        }
      }
} /*GeomObjectProcessDeletedDep*/

