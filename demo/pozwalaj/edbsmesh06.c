
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


boolean GeomObjectBSplineMeshRefinement ( GO_BSplineMesh *obj )
{
  int            onv, onhe, onfac;
  BSMvertex      *omv;
  BSMhalfedge    *omhe;
  BSMfacet       *omfac;
  int            *omvhei, *omfhei;
  double         *omvpc;
  byte           *mkcp, *mkhe, *mkfac;

  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  if ( GeomObjectCopyCurrent () )
    current_go = (geom_object*)obj;
  mkcp = mkhe = mkfac = NULL;
  if ( bsm_RefineBSMeshd ( obj->me.cpdimen, obj->degree,
                           obj->nv, obj->meshv, obj->meshvhei,
                           obj->meshvpc, obj->nhe, obj->meshhe,
                           obj->nfac, obj->meshfac, obj->meshfhei,
                           &onv, &omv, &omvhei, &omvpc, &onhe, &omhe,
                           &onfac, &omfac, &omfhei ) ) {
    mkcp  = malloc ( onv );
    mkhe  = malloc ( onhe );
    mkfac = malloc ( onfac );
    if ( !mkcp || !mkhe || !mkfac )
      goto failure;
    {
      void *sp;
      char *s;

      sp = pkv_GetScratchMemTop ();
      s = pkv_GetScratchMem ( 160 );
      if ( s ) {
        sprintf ( s, "Refining: %d vertices, %d halfedges, %d facets",
                  onv, onhe, onfac );
        SetStatusText ( s, true );
      }
      pkv_SetScratchMemTop ( sp );
    }
    if ( !bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei ) )
      goto failure;
    memset ( mkcp, MASK_CP_MOVEABLE, onv );
    memset ( mkhe, 0, onhe );
    memset ( mkfac, 0, onfac );
    GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                  onv, omv, omvhei, omvpc, onhe, omhe,
                                  onfac, omfac, omfhei, mkcp, mkhe, mkfac );
    obj->me.dlistmask = 0;
    obj->spvlist_ok = false;
    obj->special_patches_ok = false;
    return true;
  }

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshRefinement*/

boolean GeomObjectBSplineMeshDoubling ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;

  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  bsm_DoublingNum ( obj->nv, obj->meshv, obj->meshvhei, obj->nhe, obj->meshhe,
                    obj->nfac, obj->meshfac, obj->meshfhei,
                    &onv, &onhe, &onfac );

  {
    void *sp;
    char *s;

    sp = pkv_GetScratchMemTop ();
    s = pkv_GetScratchMem ( 160 );
    if ( s ) {
      sprintf ( s, "Doubling: %d vertices, %d halfedges, %d facets",
                onv, onhe, onfac );
      SetStatusText ( s, true );
    }
    pkv_SetScratchMemTop ( sp );
  }
  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_Doublingd ( obj->me.cpdimen,
                        obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                        obj->nhe, obj->meshhe,
                        obj->nfac, obj->meshfac, obj->meshfhei,
                        &onv, omv, omvhei, omvpc, &onhe, omhe,
                        &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memset ( mkcp, MASK_CP_MOVEABLE, onv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshDoubling*/

boolean GeomObjectBSplineMeshAveraging ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;

  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  bsm_AveragingNum ( obj->nv, obj->meshv, obj->meshvhei, obj->nhe, obj->meshhe,
                     obj->nfac, obj->meshfac, obj->meshfhei,
                     &onv, &onhe, &onfac );

  {
    void *sp;
    char *s;

    sp = pkv_GetScratchMemTop ();
    s = pkv_GetScratchMem ( 160 );
    if ( s ) {
      sprintf ( s, "Averaging: %d vertices, %d halfedges, %d facets",
                onv, onhe, onfac );
      SetStatusText ( s, true );
    }
    pkv_SetScratchMemTop ( sp );
  }

  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_Averagingd ( obj->me.cpdimen,
                         obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                         obj->nhe, obj->meshhe,
                         obj->nfac, obj->meshfac, obj->meshfhei,
                         &onv, omv, omvhei, omvpc, &onhe, omhe,
                         &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memset ( mkcp, MASK_CP_MOVEABLE, onv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshAveraging*/

boolean GeomObjectBSplineMeshExtractSubmesh ( GO_BSplineMesh *obj )
{
  void        *sp;
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;
  boolean     *vtag;
  int         i;

  sp = pkv_GetScratchMemTop ();
  vtag = pkv_GetScratchMem ( obj->nv*sizeof(boolean) );
  if ( !vtag )
    goto failure2;
  for ( i = 0; i < obj->nv; i++ )
    vtag[i] = (obj->mkcp[i] & marking_mask) != 0;

  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    goto failure2;
  if ( !bsm_ExtractSubmeshVNum ( obj->nv, obj->meshv,
                     obj->meshvhei, obj->nhe, obj->meshhe,
                     obj->nfac, obj->meshfac, obj->meshfhei, vtag,
                     &onv, &onhe, &onfac ) )
    goto failure2;

  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure1;
  if ( !bsm_ExtractSubmeshVd ( obj->me.cpdimen,
                         obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                         obj->nhe, obj->meshhe,
                         obj->nfac, obj->meshfac, obj->meshfhei, vtag,
                         &onv, omv, omvhei, omvpc, &onhe, omhe,
                         &onfac, omfac, omfhei ) )
    goto failure1;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure1;
  memset ( mkcp, MASK_CP_MOVEABLE, onv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  pkv_SetScratchMemTop ( sp );
  return true;

failure1:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
failure2:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplineMeshExtractSubmesh*/

boolean GeomObjectBSplineMeshRemoveCurrentVertex ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;
  int         cvn;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  cvn = obj->current_vertex[0];
  if ( obj->nv < 4 || cvn < 0 || cvn >= obj->nv )
    return false;

  bsm_RemoveVertexNum ( obj->nv, obj->meshv, obj->meshvhei,
                        obj->nhe, obj->meshhe,
                        obj->nfac, obj->meshfac, obj->meshfhei,
                        cvn, &onv, &onhe, &onfac );
  if ( onv < 0 || onhe < 0 || onfac < 0 )
    return false;
  {
    void *sp;
    char *s;

    sp = pkv_GetScratchMemTop ();
    s = pkv_GetScratchMem ( 160 );
    if ( s ) {
      sprintf ( s, "Removing vertex: %d vertices, %d halfedges, %d facets",
                onv, onhe, onfac );
      SetStatusText ( s, true );
    }
    pkv_SetScratchMemTop ( sp );
  }

  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_RemoveVertexd ( obj->me.cpdimen,
                           obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                           obj->nhe, obj->meshhe,
                           obj->nfac, obj->meshfac, obj->meshfhei,
                           cvn, &onv, omv, omvhei, omvpc, &onhe, omhe,
                           &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  if ( cvn > 0 )
    memcpy ( mkcp, obj->mkcp, cvn );
  if ( cvn < obj->nv-1 )
    memcpy ( &mkcp[cvn], &obj->mkcp[cvn+1], obj->nv-cvn-1 );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshRemoveCurrentVertex*/

boolean GeomObjectBSplineMeshContractCurrentEdge ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;
  int         che;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  che = obj->current_edge[0];
  if ( obj->nhe < 4 || che < 0 || che >= obj->nhe )
    return false;

  bsm_ContractEdgeNum ( obj->nv, obj->meshv, obj->meshvhei, obj->nhe, obj->meshhe,
                       obj->nfac, obj->meshfac, obj->meshfhei,
                       che, &onv, &onhe, &onfac );
  if ( onv < 0 || onhe < 0 || onfac < 0 )
    return false;

  {
    void *sp;
    char *s;

    sp = pkv_GetScratchMemTop ();
    s = pkv_GetScratchMem ( 160 );
    if ( s ) {
      sprintf ( s, "Contracting edge: %d vertices, %d halfedges, %d facets",
                onv, onhe, onfac );
      SetStatusText ( s, true );
    }
    pkv_SetScratchMemTop ( sp );
  }

  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( bsm_ContractEdged ( obj->me.cpdimen,
                           obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                           obj->nhe, obj->meshhe,
                           obj->nfac, obj->meshfac, obj->meshfhei,
                           che, &onv, omv, omvhei, omvpc, &onhe, omhe,
                           &onfac, omfac, omfhei ) < 0 )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memset ( mkcp, MASK_CP_MOVEABLE, onv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshContractCurrentEdge*/

boolean GeomObjectBSplineMeshShrinkCurrentEdge ( GO_BSplineMesh *obj )
{
  BSMhalfedge *mhe;
  int         che, v0, v1, cpdimen, i;
  double      *mvpc;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  che = obj->current_edge[0];
  if ( che < 0 || che >= obj->nhe )
    return false;
  cpdimen = obj->me.cpdimen;
  mhe = obj->meshhe;
  mvpc = obj->meshvpc;
  v0 = mhe[che].v0*cpdimen;
  v1 = mhe[che].v1*cpdimen;
  for ( i = 0; i < cpdimen; i++ )
    mvpc[v0+i] = mvpc[v1+i] = 0.5*(mvpc[v0+i]+mvpc[v1+i]);
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;
} /*GeomObjectBSplineMeshShrinkCurrentEdge*/

boolean GeomObjectBSplineMeshGlueEdges ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;
  
  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  mkcp = mkhe = mkfac = NULL;
  if ( !bsm_GlueTwoHalfedgesNum ( obj->nv, obj->meshv, obj->meshvhei,
                        obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
                        obj->current_edge[0], obj->current_edge[1],
                        &onv, &onhe, &onfac ) )
    return false;
  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_GlueTwoHalfedgesd ( obj->me.cpdimen,
              obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
              obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
              obj->current_edge[0], obj->current_edge[1],
              &onv, omv, omvhei, omvpc, &onhe, omhe, &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memcpy ( mkcp, obj->mkcp, onv );
  memcpy ( mkhe, obj->mkhe, onhe );
  memcpy ( mkfac, obj->mkfac, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->current_edge[0] = obj->current_edge[1] = -1;
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshGlueEdges*/

boolean GeomObjectBSplineMeshGlueEdgeLoops ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac, llgt;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;
  
  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  mkcp = mkhe = mkfac = NULL;
  llgt = bsm_HalfedgeLoopLength ( obj->nv, obj->meshv, obj->meshvhei,
                                  obj->nhe, obj->meshhe, obj->current_edge[0] );
  if ( llgt < 0 )
    return false;
  onv = obj->nv-llgt;
  onhe = obj->nhe;
  onfac = obj->nfac;
  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_GlueHalfedgeLoopsd ( obj->me.cpdimen,
              obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
              obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
              obj->current_edge[0], obj->current_edge[1],
              &onv, omv, omvhei, omvpc, &onhe, omhe, &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memcpy ( mkcp, obj->mkcp, onv );
  memcpy ( mkhe, obj->mkhe, onhe );
  memcpy ( mkfac, obj->mkfac, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->current_edge[0] = obj->current_edge[1] = -1;
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe  );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshGlueEdgeLoops*/

boolean GeomObjectBSplineMeshSealHole ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;

  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  mkcp = mkhe = mkfac = NULL;
  if ( !bsm_SealMeshHoleNum ( obj->nv, obj->meshv, obj->meshvhei,
                              obj->nhe, obj->meshhe,
                              obj->nfac, obj->meshfac, obj->meshfhei,
                              obj->current_edge[0],
                              &onv, &onhe, &onfac ) )
    return false;
  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_SealMeshHoled ( obj->me.cpdimen,
              obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
              obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
              obj->current_edge[0],
              &onv, omv, omvhei, omvpc, &onhe, omhe, &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memcpy ( mkcp, obj->mkcp, onv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->current_edge[0] = obj->current_edge[1] = -1;
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshSealHole*/

boolean GeomObjectBSplineMeshSplitBoundaryEdge ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;

  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  onv = obj->nv+1;
  onhe = obj->nhe+1;
  onfac = obj->nfac;

  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_SplitBoundaryEdged ( obj->me.cpdimen,
              obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
              obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
              obj->current_edge[0],
              &onv, omv, omvhei, omvpc, &onhe, omhe, &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memcpy ( mkcp, obj->mkcp, onv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  mkcp[onv-1] = MASK_CP_MOVEABLE;
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->current_edge[0] = obj->current_edge[1] = -1;
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshSplitBoundaryEdge*/

boolean GeomObjectBSplineMeshRemoveCurrentFacet ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;
  int         cfn;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  cfn = obj->current_facet[0];
  if ( obj->nfac < 2 || cfn < 0 || cfn >= obj->nfac )
    return false;

  bsm_RemoveFacetNum ( obj->nv, obj->meshv, obj->meshvhei, obj->nhe, obj->meshhe,
                       obj->nfac, obj->meshfac, obj->meshfhei,
                       cfn, &onv, &onhe, &onfac );
  if ( onv < 0 || onhe < 0 || onfac < 0 )
    return false;

  {
    void *sp;
    char *s;

    sp = pkv_GetScratchMemTop ();
    s = pkv_GetScratchMem ( 160 );
    if ( s ) {
      sprintf ( s, "Removing facet: %d vertices, %d halfedges, %d facets",
                onv, onhe, onfac );
      SetStatusText ( s, true );
    }
    pkv_SetScratchMemTop ( sp );
  }

  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_RemoveFacetd ( obj->me.cpdimen,
                           obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                           obj->nhe, obj->meshhe,
                           obj->nfac, obj->meshfac, obj->meshfhei,
                           cfn, &onv, omv, omvhei, omvpc, &onhe, omhe,
                           &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memset ( mkcp, MASK_CP_MOVEABLE, onv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->spvlist_ok = false;
  obj->me.dlistmask = 0;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshRemoveCurrentFacet*/

boolean GeomObjectBSplineMeshDoubleCurrentFacEdges ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;
  int         cfn;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  cfn = obj->current_facet[0];
  if ( obj->nfac < 2 || cfn < 0 || cfn >= obj->nfac )
    return false;

  bsm_FacetEdgeDoublingNum ( obj->nv, obj->meshv, obj->meshvhei,
                             obj->nhe, obj->meshhe,
                             obj->nfac, obj->meshfac, obj->meshfhei,
                             cfn, &onv, &onhe, &onfac );
  if ( onv < 0 || onhe < 0 || onfac < 0 )
    return false;
  {
    void *sp;
    char *s;

    sp = pkv_GetScratchMemTop ();
    s = pkv_GetScratchMem ( 160 );
    if ( s ) {
      sprintf ( s, "Doubling facet edges: %d vertices, %d halfedges, %d facets",
                   onv, onhe, onfac );
      SetStatusText ( s, true );
    }
    pkv_SetScratchMemTop ( sp );
  }

  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  if ( !bsm_FacetEdgeDoublingd ( obj->me.cpdimen,
                           obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                           obj->nhe, obj->meshhe,
                           obj->nfac, obj->meshfac, obj->meshfhei,
                           cfn, &onv, omv, omvhei, omvpc, &onhe, omhe,
                           &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memset ( mkcp, MASK_CP_MOVEABLE, obj->nv );
  memset ( &mkcp[obj->nv], MASK_CP_MOVEABLE | marking_mask, onv-obj->nv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshDoubleCurrentFacEdges*/

boolean GeomObjectBSplineMeshDivideFacet ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp, *mkhe, *mkfac;
  int         v0, v1;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  v0 = obj->current_vertex[0];
  v1 = obj->current_vertex[1];
  if ( v0 < 0 || v0 >= obj->nv || v1 < 0 || v1 >= obj->nv || v0 == v1 )
    return false;

  onv = obj->nv;
  onhe = obj->nhe+2;
  onfac = obj->nfac+1;
  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  mkcp = malloc ( onv );
  mkhe = malloc ( onhe );
  mkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !mkcp || !mkhe || !mkfac )
    goto failure;
  {
    void *sp;
    char *s;

    sp = pkv_GetScratchMemTop ();
    s = pkv_GetScratchMem ( 160 );
    if ( s ) {
      sprintf ( s, "Dividing a facet: %d vertices, %d halfedges, %d facets",
                   onv, onhe, onfac );
      SetStatusText ( s, true );
    }
    pkv_SetScratchMemTop ( sp );
  }

  if ( !bsm_DivideFacetd ( obj->me.cpdimen, obj->nv, obj->meshv,
                           obj->meshvhei, obj->meshvpc,
                           obj->nhe, obj->meshhe, obj->nfac, obj->meshfac,
                           obj->meshfhei, v0, v1,
                           &onv, omv, omvhei, omvpc, &onhe, omhe,
                           &onfac, omfac, omfhei ) )
    goto failure;
  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;
  memcpy ( mkcp, obj->mkcp, onv );
  memset ( mkhe, 0, onhe );
  memset ( mkfac, 0, onfac );
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp, mkhe, mkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( mkcp ) free ( mkcp );
  if ( mkhe ) free ( mkhe );
  if ( mkfac ) free ( mkfac );
  return false;
} /*GeomObjectBSplineMeshDivideFacet*/

static void _MarkVerticesInsideLoop ( int nv, BSMvertex *mv, int *mvhei,
                                      int nhe, BSMhalfedge *mhe,
                                      int nfac, BSMfacet *mfac, int *mfhei,
                                      int inhe, int v0, byte *mkcp )
{
  pkv_queue *q;
  int       deg, fhe, i, v1, he;
  byte      mask;
  
  memset ( mkcp, MASK_CP_MOVEABLE, nv );
  mask = marking_mask & ~MASK_CP_MOVEABLE;
  if ( !mask )
    return;
  q = pkv_InitQueue ( nv, sizeof(int) );
  if ( q ) {
    mkcp[v0] |= mask;
    pkv_QueueInsert ( q, &v0 );
    do {
      pkv_QueueRemoveFirst ( q, &v0 );
      deg = mv[v0].degree;
      fhe = mv[v0].firsthalfedge;
      for ( i = 0; i < deg; i++ ) {
        he = mvhei[fhe+i];
        if ( he < inhe ) {
          v1 = mhe[he].v1;
          if ( !(mkcp[v1] & mask) ) {
            mkcp[v1] |= marking_mask;
            pkv_QueueInsert ( q, &v1 );
          }
        }
      }
    } while ( !pkv_QueueEmpty ( q ) );
    PKV_FREE ( q );
  }
} /*_MarkVerticesInsideLoop*/

boolean GeomObjectBSplineMeshDoubleEdgeLoop ( GO_BSplineMesh *obj )
{
  void        *sp;
  int         inv, inhe, infac, onv, onhe, onfac;
  int         *loop, loop_length;
  boolean     *mkv, bhe;
  byte        *imkhe, *omkcp, *omkhe, *omkfac;
  int         i, j, k, v0, v1, vd, vfhe, he, he1;
  BSMvertex   *imv, *omv;
  BSMhalfedge *imhe, *omhe;
  BSMfacet    *omfac;
  int         *imvhei, *omvhei, *omfhei;
  double      *omvpc;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;

  sp = pkv_GetScratchMemTop ();
  omv = NULL;  omhe = NULL;  omfac = NULL;  omvhei = omfhei = NULL;
  omvpc = NULL;  omkcp = NULL;  omkhe = NULL;  omkfac = NULL;

  inv = obj->nv;
  inhe = obj->nhe;
  infac = obj->nfac;
  loop = pkv_GetScratchMemi ( inhe );
  mkv = pkv_GetScratchMem ( inv*sizeof(boolean) );
  if ( !loop || !mkv )
    goto failure;

        /* find the halfedge loop, if there is one */
  imv = obj->meshv;
  imvhei = obj->meshvhei;
  imhe = obj->meshhe;
  imkhe = obj->mkhe;
  memset ( mkv, false, inv*sizeof(boolean) );
  loop_length = 0;
  for ( i = 0; i < inhe; i++ )
    if ( imkhe[i] & marking_mask ) {
      loop[0] = i;
      j = imhe[i].otherhalf;
      if ( j >= 0 )
        imkhe[j] &= ~marking_mask;
      loop_length = 1;
      break;
    }
  if ( loop_length == 0 )
    goto failure;
  v0 = imhe[loop[0]].v0;
  mkv[v0] = true;
  v1 = imhe[loop[0]].v1;
  do {
    mkv[v1] = true;
    vd = imv[v1].degree;
    vfhe = imv[v1].firsthalfedge;
    for ( i = 0; i < vd; i++ ) {
      he = imvhei[vfhe+i];
      if ( imkhe[he] & marking_mask ) {
        loop[loop_length++] = he;
        j = imhe[he].otherhalf;
        if ( j >= 0 )
          imkhe[j] &= ~marking_mask;
        v1 = imhe[he].v1;
        if ( mkv[v1] && v1 != v0 )
          goto failure;
        break;
      }
    }
    if ( i >= vd )
      goto failure;
  } while ( v1 != v0 );
  if ( loop_length < 3 )
    goto failure;
        /* reorder the loop, if necessary */
        /* if there is a boundary edge in the loop */
        /* then its orientation must be preserved  */
  bhe = false;
  for ( i = 0; i < loop_length; i++ )
    if ( imhe[loop[i]].otherhalf < 0 ) {
      bhe = true;
      break;
    }
  if ( bhe ) {  /* there is a boundary edge */
    he = loop[i];
    he1 = loop[(i+1) % loop_length];
    if ( imhe[he].v1 != imhe[he1].v0 && imhe[he].v1 != imhe[he1].v1 ) {
            /* reverse the entire loop */
      for ( i = 0, j = loop_length-1;  i < j;  i++, j-- )
        { k = loop[i];  loop[i] = loop[j];  loop[j] = k; }
    }
  }
  for ( i = 0; i < loop_length-1; i++ ) {
    he = loop[i];
    he1 = loop[i+1];
    if ( imhe[he].v1 != imhe[he1].v0 && imhe[he].v0 != imhe[he1].v1 ) {
      if ( (loop[i] = imhe[he].otherhalf) < 0 )
        goto failure;
    }
  }
  he = loop[loop_length-1];
  if ( imhe[he].v1 != imhe[loop[0]].v0 ) {
    if ( (loop[loop_length-1] = imhe[he].otherhalf) < 0 )
      goto failure;
  }

printf ( "halfedge loop: " );
for ( i = 0; i < loop_length; i++ )
  printf ( "%d (%d,%d), ", loop[i], imhe[loop[i]].v0, imhe[loop[i]].v1 );
printf ( "\n" );

        /* allocate the memory and call the loop doubling procedure */
  bsm_EdgeLoopDoublingNum ( inv, inhe, infac, loop_length, &onv, &onhe, &onfac );
  omv = malloc ( onv*sizeof(BSMvertex) );
  omhe = malloc ( onhe*sizeof(BSMhalfedge) );
  omfac = malloc ( onfac*sizeof(BSMfacet) );
  omvhei = malloc ( onhe*sizeof(int) );
  omfhei = malloc ( onhe*sizeof(int) );
  omvpc = malloc ( obj->me.cpdimen*onv*sizeof(double) );
  omkcp = malloc ( onv );
  omkhe = malloc ( onhe );
  omkfac = malloc ( onfac );
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc ||
       !omkcp || !omkhe || !omkfac )
    goto failure;
  if ( !bsm_EdgeLoopDoublingd ( obj->me.cpdimen,
                                inv, imv, imvhei, obj->meshvpc,
                                inhe, imhe, infac, obj->meshfac, obj->meshfhei,
                                loop_length, loop,
                                &onv, omv, omvhei, omvpc,
                                &onhe, omhe, &onfac, omfac, omfhei ) )
    goto failure;

  obj->integrity_ok = bsm_CheckMeshIntegrity ( onv, omv, omvhei, onhe, omhe,
                                 onfac, omfac, omfhei );
  if ( !obj->integrity_ok )
    goto failure;

        /* mark vertices and loop halfedges to facilitate further editing */
  memset ( omkcp, MASK_CP_MOVEABLE, onv );
  _MarkVerticesInsideLoop ( onv, omv, omvhei, onhe, omhe, onfac, omfac, omfhei,
                            inhe, onv-1, omkcp );
  memset ( omkhe, 0, onhe );
  for ( i = 0; i < loop_length; i++ )
    omkhe[loop[i]] = marking_mask;
  memset ( omkfac, 0, onfac );

  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, omkcp, omkhe, omkfac );
  obj->me.dlistmask = 0;
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( omv ) free ( omv );
  if ( omhe ) free ( omhe );
  if ( omfac ) free ( omfac );
  if ( omvhei ) free ( omvhei );
  if ( omfhei ) free ( omfhei );
  if ( omvpc ) free ( omvpc );
  if ( omkcp ) free ( omkcp );
  if ( omkhe ) free ( omkhe );
  if ( omkfac ) free ( omkfac );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplineMeshDoubleEdgeLoop*/

