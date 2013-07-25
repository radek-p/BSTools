
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


boolean GeomObjectBSplineMeshMarkBetweenVertices ( GO_BSplineMesh *obj,
                                                   boolean mark )
{
  void        *sp;
  int         nv, nhe, nfac, *mvhei, *mfhei;
  BSMvertex   *mv;
  BSMhalfedge *mhe;
  BSMfacet    *mfac;
  byte        *mkcp;
  int         v0, v1, v2, deg, fhe, d, i, *dist;

  sp = pkv_GetScratchMemTop ();
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    goto failure;
  nv = obj->nv;
  nhe = obj->nhe;
  nfac = obj->nfac;
  v0 = obj->current_vertex[0];
  v1 = obj->current_vertex[1];
  if ( v0 < 0 || v0 >= nv || v1 < 0 || v1 >= nv )
    goto failure;
  dist = pkv_GetScratchMemi ( nv );
  if ( !dist )
    goto failure;
  mv = obj->meshv;
  mhe = obj->meshhe;
  mfac = obj->meshfac;
  mvhei = obj->meshvhei;
  mfhei = obj->meshfhei;
  mkcp = obj->mkcp;
  if ( !bsm_FindVertexDistances1 ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                   v0, dist ) )
    goto failure;
  if ( dist[v1] >= nv )
    goto failure;
  for (;;) {
    if ( mark )
      mkcp[v1] |= marking_mask;
    else
      mkcp[v1] &= ~marking_mask;
    if ( v1 == v0 )
      break;
    deg = mv[v1].degree;
    fhe = mv[v1].firsthalfedge;
    d = dist[v1]-1;
    for ( i = 0; i < deg; i++ ) {
      v2 = mhe[mvhei[fhe+i]].v1;
      if ( dist[v2] == d ) {
        v1 = v2;
        break;
      }
    }
    if ( i >= deg ) {
      i = mhe[mvhei[fhe]].facetnum;
      deg = mfac[i].degree;
      fhe = mfac[i].firsthalfedge;
      for ( i = 0; i < deg; i++ ) {
        v2 = mhe[mfhei[fhe+i]].v0;
        if ( dist[v2] == d ) {
          v1 = v2;
          break;
        }
      }
    }
    if ( dist[v2] != d )  /* just in case, should never be necessary */
      break;
  }
  obj->me.dlistmask &= ~BSM_DLM_CNET;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplineMeshMarkBetweenVertices*/

boolean GeomObjectBSplineMeshFilterPolyline ( GO_BSplineMesh *obj )
{
  void *sp;
  int         nv, nhe, nfac, *mvhei, *mfhei;
  BSMvertex   *mv;
  BSMhalfedge *mhe;
  BSMfacet    *mfac;
  int         v0, v1, v2, deg, fhe, d, i, j, *dist, npv, *ind, dim;
  double      *mvcp, *mvcps;

  sp = pkv_GetScratchMemTop ();
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    goto failure;
  nv = obj->nv;
  nhe = obj->nhe;
  nfac = obj->nfac;
  v0 = obj->current_vertex[0];
  v1 = obj->current_vertex[1];
  if ( v0 < 0 || v0 >= nv || v1 < 0 || v1 >= nv )
    goto failure;
  dist = pkv_GetScratchMemi ( 2*nv );
  if ( !dist )
    goto failure;
  ind = &dist[nv];
  mv = obj->meshv;
  mhe = obj->meshhe;
  mfac = obj->meshfac;
  mvhei = obj->meshvhei;
  mfhei = obj->meshfhei;
  if ( !bsm_FindVertexDistances1 ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                   v0, dist ) )
    goto failure;
  if ( dist[v1] >= nv )
    goto failure;
        /* find the vertex indices */
  for ( npv = 0; ;) {
    ind[npv++] = v1;
    if ( v1 == v0 )
      break;
    deg = mv[v1].degree;
    fhe = mv[v1].firsthalfedge;
    d = dist[v1]-1;
    for ( i = 0; i < deg; i++ ) {
      v2 = mhe[mvhei[fhe+i]].v1;
      if ( dist[v2] == d ) {
        v1 = v2;
        break;
      }
    }
    if ( i >= deg ) {
      i = mhe[mvhei[fhe]].facetnum;
      deg = mfac[i].degree;
      fhe = mfac[i].firsthalfedge;
      for ( i = 0; i < deg; i++ ) {
        v2 = mhe[mfhei[fhe+i]].v0;
        if ( dist[v2] == d ) {
          v1 = v2;
          break;
        }
      }
    }
    if ( i >= deg )
      goto failure;
  }
        /* now the actual filtering, simply midpoint */
  if ( npv >= 3 ) {
    if ( !GeomObjectBSplineMeshSaveCPoints ( obj ) )
      goto failure;
    mvcp = obj->meshvpc;
    mvcps = obj->savedcpoints;
    dim = obj->me.cpdimen;
    for ( i = 1; i < npv-1; i++ ) {
      v0 = ind[i-1];
      v1 = ind[i];
      v2 = ind[i+1];
      for ( j = 0; j < dim; j++ )
        mvcp[v1*dim+j] = 0.5*(mvcps[v0*dim+j]+mvcps[v2*dim+j]);
    }
  }
  obj->me.dlistmask = 0;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplineMeshFilterPolyline*/

