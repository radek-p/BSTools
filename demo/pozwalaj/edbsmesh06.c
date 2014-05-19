
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
  byte           *mkcp;

  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  if ( GeomObjectCopyCurrent () )
    current_go = (geom_object*)obj;
  mkcp = NULL;
  if ( bsm_RefineBSMeshd ( obj->me.cpdimen, obj->degree,
                           obj->nv, obj->meshv, obj->meshvhei,
                           obj->meshvpc, obj->nhe, obj->meshhe,
                           obj->nfac, obj->meshfac, obj->meshfhei,
                           &onv, &omv, &omvhei, &omvpc, &onhe, &omhe,
                           &onfac, &omfac, &omfhei ) ) {
    mkcp = malloc ( onv );
    if ( !mkcp )
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
    GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                  onv, omv, omvhei, omvpc, onhe, omhe,
                                  onfac, omfac, omfhei, mkcp );
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
  byte        *mkcp;

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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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
  byte        *mkcp;

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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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
  byte        *mkcp;
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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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
  byte        *mkcp;
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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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
  byte        *mkcp;
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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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

boolean GeomObjectBSplineMeshGlueEdgeLoops ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac, llgt;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp;
  
  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  mkcp = NULL;
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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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
  byte        *mkcp;


  if ( !bsm_CheckMeshIntegrity ( obj->nv, obj->meshv, obj->meshvhei,
                                 obj->nhe, obj->meshhe,
                                 obj->nfac, obj->meshfac, obj->meshfhei ) )
    return false;
  mkcp = NULL;
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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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
  return false;
} /*GeomObjectBSplineMeshSealHole*/

boolean GeomObjectBSplineMeshRemoveCurrentFacet ( GO_BSplineMesh *obj )
{
  int         onv, onhe, onfac;
  BSMvertex   *omv;
  BSMhalfedge *omhe;
  BSMfacet    *omfac;
  int         *omvhei, *omfhei;
  double      *omvpc;
  byte        *mkcp;
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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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
  byte        *mkcp;
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
  if ( !omv || !omhe || !omfac || !omvhei || !omfhei || !omvpc || !mkcp )
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
  GeomObjectAssignBSplineMesh ( obj, obj->me.spdimen, obj->rational,
                                onv, omv, omvhei, omvpc, onhe, omhe,
                                onfac, omfac, omfhei, mkcp );
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
  return false;
} /*GeomObjectBSplineMeshDoubleCurrentFacEdges*/

