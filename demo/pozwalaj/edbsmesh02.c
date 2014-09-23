
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
#include "edcolours.h"
#include "editor_bsm.h"

/*#define DUMPIT*/

#ifdef DUMPIT
static boolean dumpit;
#endif
/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectDrawBSplineSubdivMesh ( GO_BSplineMesh *obj )
{
  int         nv1, nhe1, nfac1, nv2, nhe2, nfac2, i;
  BSMvertex   *mv1, *mv2;
  BSMhalfedge *mhe1, *mhe2;
  BSMfacet    *mfac1, *mfac2;
  int         *mvhei1, *mvhei2, *mfhei1, *mfhei2;
  double      *mvpc1, *mvpc2;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  mv1 = mv2 = NULL;
  mhe1 = mhe2 = NULL;
  mfac1 = mfac2 = NULL;
  mvhei1 = mvhei2 = mfhei1 = mfhei2 = NULL;
  mvpc1 = mvpc2 = NULL;
  if  ( !bsm_RefineBSMeshd ( obj->me.cpdimen, obj->degree,
                             obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
                             obj->nhe, obj->meshhe,
                             obj->nfac, obj->meshfac, obj->meshfhei,
                             &nv1, &mv1, &mvhei1, &mvpc1,
                             &nhe1, &mhe1, &nfac1, &mfac1, &mfhei1 ) )
    goto failure;
  for ( i = 1; i < obj->subdivl; i++ ) {
    if ( !bsm_RefineBSMeshd ( obj->me.cpdimen, obj->degree,
             nv1, mv1, mvhei1, mvpc1, nhe1, mhe1, nfac1, mfac1, mfhei1,
             &nv2, &mv2, &mvhei2, &mvpc2, &nhe2, &mhe2, &nfac2, &mfac2, &mfhei2 ) )
      goto failure;
    free ( mv1 );     mv1 = mv2;
    free ( mvhei1 );  mvhei1 = mvhei2;
    free ( mvpc1 );   mvpc1 = mvpc2;
    free ( mhe1 );    mhe1 = mhe2;
    free ( mfac1 );   mfac1 = mfac2;
    free ( mfhei1 );  mfhei1 = mfhei2;
    nv1 = nv2;  nhe1 = nhe2;  nfac1 = nfac2;
  }
  glColor3fv ( OBJC_MESH_SURF_BOUNDARY );
  GeomObjectDrawMeshEdges ( obj->me.cpdimen, obj->me.spdimen,
                            nv1, mv1, mvpc1, mvhei1,
                            nhe1, mhe1, nfac1, mfac1, mfhei1,
                            MASK_HE_BOUNDARY );
  glColor3fv ( OBJC_MESH_SURF_INNER );
  GeomObjectDrawMeshEdges ( obj->me.cpdimen, obj->me.spdimen,
                            nv1, mv1, mvpc1, mvhei1,
                            nhe1, mhe1, nfac1, mfac1, mfhei1,
                            MASK_HE_INNER );
failure:
  if ( mv1 ) free ( mv1 );
  if ( mvhei1 ) free ( mvhei1 );
  if ( mvpc1 ) free ( mvpc1 );
  if ( mhe1 ) free ( mhe1 );
  if ( mfac1 ) free ( mfac1 );
  if ( mfhei1 ) free ( mfhei1 );
} /*GeomObjectDrawBSplineSubdivMesh*/

/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectDrawMeshPatch ( int d, int *vertnum, int *mtab, void *usrptr )
{
  void           *sp;
  GO_BSplineMesh *obj;
  double         *mcp, *cp, *bp, *knu, *knv;
  int            dim, degree, lkn, i, j, ku, kv, pitch;

  obj = usrptr;
  sp = pkv_GetScratchMemTop ();
  dim = obj->me.cpdimen;
  degree = obj->degree;
  if ( degree+1 != d )
    return;
  lkn = 2*degree+1;
  knu = pkv_GetScratchMemd ( lkn+1 );
  knv = pkv_GetScratchMemd ( lkn+1 );
  cp = pkv_GetScratchMemd ( d*d*dim );
  bp = pkv_GetScratchMemd ( d*d*dim );
  if ( knu && knv && cp && bp ) {
    for ( i = 0; i <= lkn; i++ )
      knu[i] = knv[i] = (double)i;
    mcp = obj->meshvpc;
    for ( i = j = 0;  i < d*d;  i++, j += dim )
      memcpy ( &cp[j], &mcp[dim*vertnum[i]], dim*sizeof(double) );
    pitch = d*dim;
    mbs_BSPatchToBezd ( dim, degree, lkn, knu, degree, lkn, knv,
                        pitch, cp, &ku, NULL, NULL, &kv, NULL, NULL, pitch, bp );
    glColor3fv ( OBJC_MESH_SURF_BOUNDARY );
    DrawBezierPatchWF ( dim, obj->me.spdimen, degree, degree, pitch, bp,
                        1, 1, 8*obj->density, 8*obj->density,
                        true, true, true, true );
    glColor3fv ( OBJC_MESH_SURF_INNER );
    DrawBezierPatchWF ( dim, obj->me.spdimen, degree, degree, pitch, bp,
                        obj->density, obj->density, 8*obj->density, 8*obj->density,
                        false, false, false, false );
#ifdef DUMPIT
    if ( dumpit )
      GeomObjectDumpBezierPatch ( obj->me.spdimen, obj->me.cpdimen, obj->rational,
                                  degree, degree, bp );
#endif
  }
  pkv_SetScratchMemTop ( sp );
} /*GeomObjectDrawMeshPatch*/

void GeomObjectDrawBSplineMeshPatches ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;

  bsm_FindRegularSubnets ( obj->nv, obj->meshv, obj->meshvhei,
                           obj->nhe, obj->meshhe,
                           obj->nfac, obj->meshfac, obj->meshfhei,
                           obj->degree+1, obj, GeomObjectDrawMeshPatch );
} /*GeomObjectDrawBSplineMeshPatches*/

/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectBSplineMeshOutputBezierPatch ( int n, int m, const double *cp,
                                              void *usrptr )
{
  GO_BSplineMesh *obj;

  obj = (GO_BSplineMesh*)usrptr;
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  glColor3fv ( OBJC_MESH_SURF_BOUNDARY );
  DrawBezierPatchWF ( obj->me.cpdimen, obj->me.spdimen, n, m,
                      (n+1)*obj->me.cpdimen, cp,
                      1, 1, 8*obj->density, 8*obj->density,
                      true, true, true, true );
  if ( obj->fill_G1 )
    glColor3fv ( OBJC_MESH_SURF_G1 );
  else if ( obj->fill_G2 )
    glColor3fv ( OBJC_MESH_SURF_G2 );
  else
    glColor3fv ( OBJC_MESH_SURF_G1Q2 );
  DrawBezierPatchWF ( obj->me.cpdimen, obj->me.spdimen, n, m,
                      (n+1)*obj->me.cpdimen, cp,
                      obj->density, obj->density, 8*obj->density, 8*obj->density,
                      false, false, false, false );
#ifdef DUMPIT
  if ( dumpit )
    GeomObjectDumpBezierPatch ( obj->me.spdimen, obj->me.cpdimen, obj->rational,
                                n, m, cp );
#endif
} /*GeomObjectBSplineMeshOutputBezierPatch*/

boolean GeomObjectBSplineMeshFindSpecialVertices ( GO_BSplineMesh *obj, int d )
{
  bsm_special_elem_list *spvlist;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  if ( !obj->spvlist_ok ) {
    spvlist = &obj->spvlist;
    if ( spvlist->spel )   { free ( spvlist->spel );  spvlist->spel = NULL; }
    if ( spvlist->spvert ) { free ( spvlist->spvert );  spvlist->spvert = NULL; }
    if ( !bsm_CountSpecialVSubnets ( obj->nv, obj->meshv, obj->meshvhei,
              obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
              d, &spvlist->nspecials, &spvlist->nspvert ) )
      return false;
    spvlist->spel = malloc ( spvlist->nspecials*sizeof(bsm_special_el) );
    spvlist->spvert = malloc ( spvlist->nspvert*sizeof(int) );
    if ( !spvlist->spel || !spvlist->spvert )
      goto failure;
    if ( !bsm_FindSpecialVSubnetList ( obj->nv, obj->meshv, obj->meshvhei,
              obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
              d, false, spvlist ) )
      goto failure;
    obj->spvlist_ok = true;
  }
  return true;

failure:
  if ( spvlist->spel )   { free ( spvlist->spel );  spvlist->spel = NULL; }
  if ( spvlist->spvert ) { free ( spvlist->spvert );  spvlist->spvert = NULL; }
  return false;
} /*GeomObjectBSplineMeshFindSpecialVertices*/

void GeomObjectDrawBSplineMeshCubicHoleFillingPatches ( GO_BSplineMesh *obj )
{
  bsm_special_elem_list *spvlist;
  bsm_special_el        *spel;
  int                   i, deg, d, nsp, *spvert;
  double                *spv;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  if ( obj->degree != 3 || obj->subdivision )
    return;
  if ( obj->me.dlistmask & BSM_DLM_HOLEFILL )
    glCallList ( obj->me.displaylist+BSM_DL_HOLEFILL );
  else if ( obj->special_patches_ok ) {
    glNewList ( obj->me.displaylist+BSM_DL_HOLEFILL, GL_COMPILE_AND_EXECUTE );
    deg = obj->special_deg;
    d = 3*(deg+1)*(deg+1);
    for ( i = 0, spv = obj->special_patches;
          i < obj->nspecial_patches;
          i++, spv = &spv[d] )
      GeomObjectBSplineMeshOutputBezierPatch ( deg, deg, spv, obj );
    glEndList ();
    obj->me.dlistmask |= BSM_DLM_HOLEFILL;
  }
  else if ( GeomObjectBSplineMeshFindSpecialVertices ( obj, 2 ) ) {
    glNewList ( obj->me.displaylist+BSM_DL_HOLEFILL, GL_COMPILE_AND_EXECUTE );
    spvlist = &obj->spvlist;
    spel    = spvlist->spel;
    spvert  = spvlist->spvert;
    nsp     = spvlist->nspecials;
    for ( i = 0; i < nsp; i++ )
      if ( spel[i].el_type == 0 && spel[i].snet_rad == 2 )
        GeomObjectBSplineMeshFillBicubicHole ( obj, spel[i].degree,
                    &spvert[spel[i].first_snet_vertex],
                    GeomObjectBSplineMeshOutputBezierPatch );
    glEndList ();
    obj->me.dlistmask |= BSM_DLM_HOLEFILL;
  }
} /*GeomObjectDrawBSplineMeshCubicHoleFillingPatches*/

/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectDrawBSplineMeshSurf ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  if ( obj->me.dlistmask & BSM_DLM_SURF )
    glCallList ( obj->me.displaylist+BSM_DL_SURF );
  else {
    glNewList ( obj->me.displaylist+BSM_DL_SURF, GL_COMPILE_AND_EXECUTE );
    if ( obj->subdivision )
      GeomObjectDrawBSplineSubdivMesh ( obj );
    else
      GeomObjectDrawBSplineMeshPatches ( obj );
    glEndList ();
    obj->me.dlistmask |= BSM_DLM_SURF;
  }
} /*GeomObjectDrawBSplineMeshSurf*/

/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectDrawBSplineMNet ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  if ( obj->me.dlistmask & BSM_DLM_CNET )
    glCallList ( obj->me.displaylist+BSM_DL_CNET );
  else {
    glNewList ( obj->me.displaylist+BSM_DL_CNET, GL_COMPILE_AND_EXECUTE );
    if ( !obj->mkhe || !marking_mask ) {
      glColor3fv ( OBJC_MESH_EDGE_BOUNDARY );
      GeomObjectDrawMeshEdges ( obj->me.cpdimen, obj->me.spdimen,
                                obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                MASK_HE_BOUNDARY );
      glColor3fv ( OBJC_MESH_EDGE_INNER );
      GeomObjectDrawMeshEdges ( obj->me.cpdimen, obj->me.spdimen,
                                obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                MASK_HE_INNER );
    }
    else {
      glColor3fv ( OBJC_MESH_EDGE_BOUNDARY );
      GeomObjectDrawMeshHalfedges ( obj->me.cpdimen, obj->me.spdimen,
                                obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                MASK_HE_BOUNDARY,
                                obj->mkhe, marking_mask, true, true, true );
      GeomObjectDrawMeshHalfedges ( obj->me.cpdimen, obj->me.spdimen,
                                obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                MASK_HE_BOUNDARY,
                                obj->mkhe, marking_mask, false, true, false );
      glLineWidth ( 2.0 );
      glColor3fv ( OBJC_MESH_EDGE_MKBOUNDARY );
      GeomObjectDrawMeshHalfedges ( obj->me.cpdimen, obj->me.spdimen,
                                obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                MASK_HE_BOUNDARY,
                                obj->mkhe, marking_mask, false, false, true );
      glColor3fv ( OBJC_MESH_EDGE_MKINNER );
      GeomObjectDrawMeshHalfedges ( obj->me.cpdimen, obj->me.spdimen,
                                obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                MASK_HE_INNER,
                                obj->mkhe, marking_mask, false, false, true );
      glLineWidth ( 1.0 );
      glColor3fv ( OBJC_MESH_EDGE_INNER );
      GeomObjectDrawMeshHalfedges ( obj->me.cpdimen, obj->me.spdimen,
                                obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                MASK_HE_INNER,
                                obj->mkhe, marking_mask, true, false, true );
    }
    if ( obj->blending ) {
      glColor3fv ( OBJC_MESH_VERTEX_MARKED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv, obj->meshvpc,
                    obj->mkcp, marking_mask,
                    MASK_MARKED | MASK_CP_BOUNDARY | MASK_CP_MOVEABLE );
      glColor3fv ( OBJC_MESH_VERTEX_BOUNDARY );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv, obj->meshvpc,
                    obj->mkcp, 0, MASK_CP_BOUNDARY | MASK_CP_MOVEABLE );
      glColor3fv ( OBJC_MESH_VERTEX_BMARKED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv, obj->meshvpc,
                    obj->mkcp, marking_mask, MASK_MARKED | MASK_CP_MOVEABLE );
      glColor3fv ( OBJC_MESH_VERTEX_BOUNDARY );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv, obj->meshvpc,
                    obj->mkcp, 0, MASK_CP_MOVEABLE );
      glColor3fv ( OBJC_MESH_VERTEX_CHANGED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv, obj->meshvpc,
                    obj->mkcp, 0, MASK_CP_SPECIAL );
    }
    else {
      glColor3fv ( OBJC_MESH_VERTEX_MARKED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv, obj->meshvpc,
                    obj->mkcp, marking_mask, MASK_MARKED | MASK_CP_MOVEABLE );
      glColor3fv ( OBJC_MESH_VERTEX_UNMARKED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, obj->nv, obj->meshvpc,
                    obj->mkcp, 0, MASK_CP_MOVEABLE );
    }
    glEndList ();
    obj->me.dlistmask |= BSM_DLM_CNET;
  }
} /*GeomObjectDrawBSplineMNet*/

/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectDrawBSplineMElements ( GO_BSplineMesh *obj )
{
  int     i;
  GLfloat *vcolour[BSM_NCURRENT_VERTICES] =
            {OBJC_MESH_VERTEX_SEL0,OBJC_MESH_VERTEX_SEL1};
  GLfloat *ecolour[BSM_NCURRENT_EDGES] =
            {OBJC_MESH_EDGE_SEL0,OBJC_MESH_EDGE_SEL1};
  GLfloat *fcolour[BSM_NCURRENT_FACETS] =
            {OBJC_MESH_FACET_SEL0,OBJC_MESH_FACET_SEL1};

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  if ( obj->me.dlistmask & BSM_DLM_ELEM )
    glCallList ( obj->me.displaylist+BSM_DL_ELEM );
  else {
    glNewList ( obj->me.displaylist+BSM_DL_ELEM, GL_COMPILE_AND_EXECUTE );
    for ( i = 0; i < BSM_NCURRENT_FACETS; i++ )
      if ( obj->current_facet[i] >= 0 && obj->current_facet[i] < obj->nfac ) {
        glColor3fv ( fcolour[i] );
        GeomObjectDrawMeshFacet ( obj->me.cpdimen, obj->me.spdimen,
                                  obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                  obj->nhe, obj->meshhe,
                                  obj->nfac, obj->meshfac, obj->meshfhei,
                                  obj->current_facet[i] );
      }
    for ( i = 0; i < BSM_NCURRENT_EDGES; i++ )
      if ( obj->current_edge[i] >= 0 && obj->current_edge[i] < obj->nhe ) {
        glColor3fv ( ecolour[i] );
        GeomObjectDrawMeshHalfedge ( obj->me.cpdimen, obj->me.spdimen,
                                     obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                     obj->nhe, obj->meshhe,
                                     obj->nfac, obj->meshfac, obj->meshfhei,
                                     obj->current_edge[i] );
      }
    for ( i = 0; i < BSM_NCURRENT_VERTICES; i++ )
      if ( obj->current_vertex[i] >= 0 && obj->current_vertex[i] < obj->nv ) {
        glColor3fv ( vcolour[i] );
        GeomObjectDrawMeshVertex ( obj->me.cpdimen, obj->me.spdimen,
                                   obj->nv, obj->meshv, obj->meshvpc, obj->meshvhei,
                                   obj->nhe, obj->meshhe,
                                   obj->nfac, obj->meshfac, obj->meshfhei,
                                   obj->current_vertex[i] );
      }
    glEndList ();
    obj->me.dlistmask |= BSM_DLM_ELEM;
  }
} /*GeomObjectDrawBSplineMElements*/

/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectBSMDrawVSpecial ( int d, int k, int *vertnum, int *mtab,
                                 void *usrptr )
{
  GO_BSplineMesh *obj;
  double      *mvpc;
  int         i, j, l, dd, dd2, cpdim, spdim;

  obj = (GO_BSplineMesh*)usrptr;
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  mvpc  = obj->meshvpc;
  cpdim = obj->me.cpdimen;
  spdim = obj->me.spdimen;
  dd = d+d+1;
  dd2 = dd*dd;
  for ( i = 0; i < k; i++ ) {
    for ( j = 0; j <= dd; j += 2 ) {
      glBegin ( GL_LINE_STRIP );
        for ( l = 0; l <= dd; l += 2 )
          DrawAVertex ( cpdim, spdim, &mvpc[cpdim*mtab[j*dd+l]] );
      glEnd ();
    }
    for ( l = 0; l <= dd; l += 2 ) {
      glBegin ( GL_LINE_STRIP );
        for ( j = 0; j <= dd; j += 2 )
          DrawAVertex ( cpdim, spdim, &mvpc[cpdim*mtab[j*dd+l]] );
      glEnd ();
    }
    mtab = &mtab[dd2];
  }
} /*GeomObjectBSMDrawVSpecial*/

void GeomObjectBSMDrawFSpecial ( int d, int k, int *vertnum, int *mtab,
                                 void *usrptr )
{
  GO_BSplineMesh *obj;
  double      *mvpc;
  int         i, j, l, dd, dd2, cpdim, spdim;

  obj = (GO_BSplineMesh*)usrptr;
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  mvpc  = obj->meshvpc;
  cpdim = obj->me.cpdimen;
  spdim = obj->me.spdimen;
  dd = d+d+1;
  dd2 = dd*(dd+2);
  for ( i = 0; i < k; i++ ) {
    for ( j = 0; j <= dd+2; j += 2 ) {
      glBegin ( GL_LINE_STRIP );
        for ( l = 0; l <= dd; l += 2 )
          DrawAVertex ( cpdim, spdim, &mvpc[cpdim*mtab[j*dd+l]] );
      glEnd ();
    }
    for ( l = 0; l <= dd; l += 2 ) {
      glBegin ( GL_LINE_STRIP );
        for ( j = 0; j <= dd+2; j += 2 )
          DrawAVertex ( cpdim, spdim, &mvpc[cpdim*mtab[j*dd+l]] );
      glEnd ();
    }
    mtab = &mtab[dd2];
  }
} /*GeomObjectBSMDrawFSpecial*/

void GeomObjectDrawBSplineMSpecialElements ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  if ( obj->me.dlistmask & BSM_DLM_SPECIAL )
    glCallList ( obj->me.displaylist+BSM_DL_SPECIAL );
  else {
    glNewList ( obj->me.displaylist+BSM_DL_SPECIAL, GL_COMPILE_AND_EXECUTE );
    if ( !obj->subdivision ) {
      glColor3fv ( OBJC_MESH_SPECIAL_VSUBNET );
      bsm_FindSpecialVSubnets ( obj->nv, obj->meshv, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                obj->degree-1, obj, GeomObjectBSMDrawVSpecial );
      glColor3fv ( OBJC_MESH_SPECIAL_FSUBNET );
      bsm_FindSpecialFSubnets ( obj->nv, obj->meshv, obj->meshvhei,
                                obj->nhe, obj->meshhe,
                                obj->nfac, obj->meshfac, obj->meshfhei,
                                obj->degree-1, obj, GeomObjectBSMDrawFSpecial );
    }
    glEndList ();
    obj->me.dlistmask |= BSM_DLM_SPECIAL;
  }
} /*GeomObjectDrawBSplineMSpecialElements*/

/* ///////////////////////////////////////////////////////////////////////// */
void GeomObjectDisplayBSplineMesh ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
#ifdef DUMPIT
  if ( !obj->me.dlistmask )
    dumpit = GeomObjectOpenDumpFile ();
  else
    dumpit = false;
#endif
  if ( obj->view_special )
    GeomObjectDrawBSplineMSpecialElements ( obj );
  GeomObjectDrawBSplineMElements ( obj );
  if ( !obj->subdivision && obj->view_holefill )
    GeomObjectDrawBSplineMeshCubicHoleFillingPatches ( obj );
  if ( obj->view_surf )
    GeomObjectDrawBSplineMeshSurf ( obj );
  if ( obj->view_cnet )
    GeomObjectDrawBSplineMNet ( obj );
#ifdef DUMPIT
  if ( dumpit )
    GeomObjectCloseDumpFile ();
#endif
} /*GeomObjectDisplayBSplineMesh*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean GeomObjectBSplineMeshSetDensLevel ( GO_BSplineMesh *obj, int dens )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  if ( obj->subdivision ) {
    if ( dens > 0 && dens <= MAX_SUBDIV_LEVEL ) {
      obj->subdivl = dens;
      obj->me.dlistmask &= ~BSM_DLM_SURF;
      return true;
    }
    else
      return false;
  }
  else {
    if ( dens > 0 && dens <= MAX_PNET_DENSITY ) {
      obj->density = dens;
      obj->me.dlistmask &= ~(BSM_DLM_SURF | BSM_DLM_HOLEFILL);
      return true;
    }
    else
      return false;
  }
} /*GeomObjectBSplineMeshSetDensLevel*/

boolean GeomObjectBSplineMeshSetSubdiv ( GO_BSplineMesh *obj, boolean on )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  obj->subdivision = on;
  if ( on )
    obj->blending = false;
  obj->me.dlistmask &= ~BSM_DLM_SURF;
  obj->spvlist_ok = false;
  return true;
} /*GeomObjectBSplineMeshSetSubdiv*/

static boolean _GeomObjectBSMMarkBZVertBlending ( GO_BSplineMesh *obj )
{
  void        *sp;
  int         i, d, nv, nhe;
  byte        *mkcp;
  BSMvertex   *mv;
  BSMhalfedge *mhe;
  int         *mvhei;
  boolean     result;
  char        *vtag;

  sp = pkv_GetScratchMemTop ();
  result = false;
  nv = obj->nv;
  vtag = pkv_GetScratchMem ( nv );
  if ( vtag ) {
    mv = obj->meshv;
    mvhei = obj->meshvhei;
    nhe = obj->nhe;
    mhe = obj->meshhe;
    mkcp = obj->mkcp;
    d = obj->degree;
    bsm_TagBoundaryZoneVertices ( nv, mv, mvhei, nhe, mhe, d, vtag );
    for ( i = 0; i < nv; i++ )
      if ( vtag[i] < d ) {
        result = true;
        mkcp[i] |= MASK_CP_BOUNDARY;
      }
      else
        mkcp[i] &= ~MASK_CP_BOUNDARY;
  }
  pkv_SetScratchMemTop ( sp );
  return result;
} /*_GeomObjectSBMMarkBZVertBlending*/

static void _GeomObjectBSMMarkBZVertOther ( GO_BSplineMesh *obj )
{
  int  i, nv;
  byte *mkcp;

  nv = obj->nv;
  mkcp = obj->mkcp;
  for ( i = 0; i < nv; i++ )
    mkcp[i] &= ~MASK_CP_BOUNDARY;
} /*_GeomObjectBSMMarkBZVertOther*/

boolean GeomObjectBSplineMeshSetBlending ( GO_BSplineMesh *obj, boolean bl )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  obj->me.dlistmask &= !BSM_DLM_CNET;
  obj->blending = bl;
  if ( bl ) {
    _GeomObjectBSMMarkBZVertBlending ( obj );
    obj->subdivision = false;
    obj->degree = 3;
    obj->fill_Coons = false;  obj->fill_Bezier = true;
    obj->fill_G1 = obj->fill_G1Q2 = false;  obj->fill_G2 = true;
  }
  else
    _GeomObjectBSMMarkBZVertOther ( obj );
  obj->me.dlistmask &= ~BSM_DLM_SURF;
  obj->spvlist_ok = false;
  return bl;
} /*GeomObjectBSplineMeshSetBlending*/

