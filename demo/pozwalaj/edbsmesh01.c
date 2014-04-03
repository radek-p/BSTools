
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

void GeomObjectAssignBSplineMesh ( GO_BSplineMesh *obj,
                      int spdimen, boolean rational,
                      int nv, BSMvertex *mv, int *mvhei, double *mvpc,
                      int nhe, BSMhalfedge *mhe,
                      int nfac, BSMfacet *mfac, int *mfhei, byte *mkcp )
{
  int nnv;
  if ( obj->meshv ) { free ( obj->meshv );  nnv = obj->nv; }
    else nnv = -1;
  if ( obj->meshvhei ) free ( obj->meshvhei );
  if ( obj->meshvpc ) free ( obj->meshvpc );
  if ( obj->meshhe ) free ( obj->meshhe );
  if ( obj->meshfac ) free ( obj->meshfac );
  if ( obj->meshfhei ) free ( obj->meshfhei );
  if ( obj->mkcp ) free ( obj->mkcp );
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  if ( obj->special_patches ) {
    free ( obj->special_patches );
    obj->special_patches = NULL;
  }
  obj->nspecial_patches = 0;
  obj->special_patches_ok = false;
  obj->maxnv = obj->nv = nv;
  obj->maxnhe = obj->nhe = nhe;
  obj->maxnfac = obj->nfac = nfac;
  obj->meshv = mv;
  obj->meshvhei = mvhei;
  obj->meshvpc = mvpc;
  obj->meshhe = mhe;
  obj->meshfac = mfac;
  obj->meshfhei = mfhei;
  obj->mkcp = mkcp;
  obj->rational = rational;
  obj->me.spdimen = spdimen;
  obj->me.dlistmask = 0;
  if ( rational )
    obj->me.cpdimen = spdimen+1;
  else
    obj->me.cpdimen = spdimen;
  if ( nv != nnv ) {
    if ( nv > 10000 ) {
      obj->density = 1;
      obj->subdivl = 1;
    }
    else if ( nv > 1000 ) {
      obj->density = 2;
      obj->subdivl = 2;
    }
  }
  obj->spvlist_ok = false;
  GeomObjectBSplineMeshSetBlending ( obj, obj->blending );
} /*GeomObjectAssignBSplineMesh*/

boolean GeomObjectExtendBSplineMesh ( GO_BSplineMesh *obj,
                      int spdimen, boolean rational,
                      int nv, BSMvertex *mv, int *mvhei, double *mvpc, 
                      int nhe, BSMhalfedge *mhe,
                      int nfac, BSMfacet *mfac, int *mfhei, byte *mkcp )
{
  int         nnv, nnhe, nnfac, *nmvhei, *nmfhei;
  BSMvertex   *nmv;
  BSMhalfedge *nmhe;
  BSMfacet    *nmfac;
  double      *nmvpc;
  byte        *nmkcp;

  if ( obj->me.spdimen != spdimen || obj->rational != rational )
    return false;
  obj->special_patches_ok = false;
  nnv = nv + obj->nv;
  nnhe = nhe + obj->nhe;
  nnfac = nfac + obj->nfac;
  nmvhei = malloc ( nnhe*sizeof(int) );
  nmfhei = malloc ( nnhe*sizeof(int) );
  nmv = malloc ( nnv*sizeof(BSMvertex) );
  nmhe = malloc ( nnhe*sizeof(BSMhalfedge) );
  nmfac = malloc ( nnfac*sizeof(BSMfacet) );
  nmvpc = malloc ( nnv*obj->me.cpdimen*sizeof(double) );
  nmkcp = malloc ( nnv );
  if ( !nmvhei || !nmfhei || !nmv || !nmhe || !nmfac || !nmvpc || !nmkcp )
    goto failure;

  bsm_MergeMeshesd ( spdimen, obj->nv, obj->meshv, obj->meshvhei, obj->meshvpc,
         obj->nhe, obj->meshhe, obj->nfac, obj->meshfac, obj->meshfhei,
         nv, mv, mvhei, mvpc, nhe, mhe, nfac, mfac, mfhei,
         &nnv, nmv, nmvhei, nmvpc, &nnhe, nmhe, &nnfac, nmfac, nmfhei );
  free ( mv );
  free ( mhe );
  free ( mfac );
  free ( mvhei );
  free ( mfhei );
  free ( mkcp );
  free ( mvpc );
  GeomObjectAssignBSplineMesh ( obj, spdimen, rational, nnv, nmv, nmvhei, nmvpc,
                                nnhe, nmhe, nnfac, nmfac, nmfhei, nmkcp );
  memset ( nmkcp, MASK_CP_MOVEABLE, nnv-nv );
  memset ( &nmkcp[nnv-nv], MASK_CP_MOVEABLE | marking_mask, nv );
  obj->spvlist_ok = false;
  return true;

failure:
  free ( mv );
  free ( mhe );
  free ( mfac );
  free ( mvhei );
  free ( mfhei );
  free ( mkcp );
  free ( mvpc );
  if ( nmvhei ) free ( nmvhei );
  if ( nmfhei ) free ( nmfhei );
  if ( nmv ) free ( nmv );
  if ( nmhe ) free ( nmhe );
  if ( nmfac ) free ( nmfac );
  if ( nmvpc ) free ( nmvpc );
  if ( nmkcp ) free ( nmkcp );
  obj->spvlist_ok = false;
  return false;
} /*GeomObjectExtendBSplineMesh*/
 
boolean GeomObjectInitBSplineMesh ( GO_BSplineMesh *obj,
                                    char spdimen, boolean rational )
{
  obj->me.obj_type = GO_BSPLINE_MESH;
  obj->me.ident = -1;
  obj->me.spdimen = spdimen;
  obj->rational = rational;
  if ( rational )
    obj->me.cpdimen = spdimen+1;
  else
    obj->me.cpdimen = spdimen;
  obj->subdivision = obj->blending = false;
  obj->me.active = false;
  obj->me.name[0] = 0;
  if ( !GeomObjectBSplineMeshInitKGon ( obj, 4, false ) )
    return false;
  obj->degree = 3;
  obj->density = 6;
  obj->subdivl = 4;
  obj->savedsize = 0;
  obj->view_surf = obj->view_cnet = obj->view_holefill = true;
  obj->current_vertex[0] = obj->current_vertex[1] =
    obj->current_edge[0] = obj->current_edge[1] =
    obj->current_facet[0] = obj->current_facet[1] = -1;
  obj->fill_Coons = obj->fill_G2 = true;
  obj->fill_Bezier = obj->fill_G1 = obj->fill_G1Q2 = false;
  obj->sl_g1q2param = 0.5;
  obj->g1q2param = 10.0;
  obj->me.displaylist = glGenLists ( BSM_NDL );
                      /* default parameter values for blending mesh optimization */
  obj->bl_constr = obj->bl_shape_only = false;
  obj->bsm_bl_C = 0.01;
  obj->nkn1 = 6;
  obj->nkn2 = 8;
  obj->nlevels = 0;
  obj->nblocks = 1;
  obj->maxit = 100;
  obj->bl_use_coarse = false;
  obj->coarse_name[0] = 0;
  obj->me.colour[0] = obj->me.colour[1] = 1.0;
  obj->me.colour[2] = 0.0;
  obj->me.display_pretrans = false;
  IdentTrans3d ( &obj->me.pretrans );
  obj->spvlist_ok = false;
  obj->special_patches_ok = false;
  return true;
} /*GeomObjectInitBSplineMesh*/

geom_object *GeomObjectAddBSplineMesh ( const char *name,
                                        char spdimen, boolean rational )
{
  GO_BSplineMesh *obj;

  obj = malloc ( sizeof(GO_BSplineMesh) );
  if ( obj ) {
    memset ( obj, 0, sizeof(GO_BSplineMesh) );
    if ( !GeomObjectInitBSplineMesh ( obj, spdimen, rational ) ) {
      free ( obj );
      return NULL;
    }
    strncpy ( obj->me.name, name, MAX_NAME_LENGTH+1 );
    if ( !first_go )
      first_go = last_go = &obj->me;
    else {
      obj->me.prev = last_go;
      last_go->next = &obj->me;
      last_go = &obj->me;
    }
    current_go = &obj->me;
    return &obj->me;
  }
  else
    return NULL;
} /*GeomObjectAddBSplineMesh*/

geom_object *GeomObjectCopyBSplineMesh ( GO_BSplineMesh *obj )
{
  GO_BSplineMesh *copy;
  BSMvertex      *mv;
  BSMhalfedge    *mhe;
  BSMfacet       *mfac;
  int            *mvhei, *mfhei;
  double         *mvpc;
  byte           *mkcp;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return NULL;
  copy = malloc ( sizeof(GO_BSplineMesh) );
  if ( copy ) {
    memset ( copy, 0, sizeof(GO_BSplineMesh) );
    copy->me.obj_type = GO_BSPLINE_MESH;
    copy->me.ident = -1;
    copy->me.spdimen = obj->me.spdimen;
    copy->me.cpdimen = obj->me.cpdimen;
    copy->me.active = false;
    strcpy ( copy->me.name, obj->me.name );
    copy->me.display_pretrans = false;
    copy->me.pretrans = obj->me.pretrans;
    copy->maxnv = obj->maxnv;
    copy->maxnhe = obj->maxnhe;
    copy->maxnfac = obj->maxnfac;
    mv    = malloc ( obj->maxnv*sizeof(BSMvertex) );
    mvhei = malloc ( obj->maxnhe*sizeof(int) );
    mvpc  = malloc ( obj->maxnv*obj->me.cpdimen*sizeof(double) );
    mhe   = malloc ( obj->maxnhe*sizeof(BSMhalfedge) );
    mfac  = malloc ( obj->maxnfac*sizeof(BSMfacet) );
    mfhei = malloc ( obj->maxnhe*sizeof(int) );
    mkcp     = malloc ( obj->maxnv );
    if ( !mv || !mvhei || !mvpc || !mhe || !mfac || !mfhei || !mkcp ) {
      if ( mv )    free ( mv );
      if ( mvhei ) free ( mvhei );
      if ( mvpc )  free ( mvpc );
      if ( mhe )   free ( mhe );
      if ( mfac )  free ( mfac );
      if ( mfhei ) free ( mfhei );
      if ( mkcp )  free ( mkcp );
      free ( copy );
      return NULL;
    }
    copy->degree = obj->degree;
    memcpy ( mv, obj->meshv, obj->nv*sizeof(BSMvertex) );
    memcpy ( mvhei, obj->meshvhei, obj->nhe*sizeof(int) );
    memcpy ( mvpc, obj->meshvpc, obj->nv*obj->me.cpdimen*sizeof(double) );
    memcpy ( mhe, obj->meshhe, obj->nhe*sizeof(BSMhalfedge) );
    memcpy ( mfac, obj->meshfac, obj->nfac*sizeof(BSMfacet) );
    memcpy ( mfhei, obj->meshfhei, obj->nhe*sizeof(int) );
    memcpy ( mkcp, obj->mkcp, obj->nv*sizeof(byte) );
    copy->blending = obj->blending;
    GeomObjectAssignBSplineMesh ( copy,
                      obj->me.spdimen, obj->rational,
                      obj->nv, mv, mvhei, mvpc, obj->nhe, mhe,
                      obj->nfac, mfac, mfhei, mkcp );
    copy->subdivision = obj->subdivision;
    copy->density = obj->density;
    copy->subdivl = obj->subdivl;
    copy->view_surf = copy->view_cnet = copy->view_holefill = true;
    copy->fill_Coons = obj->fill_Coons;
    copy->fill_Bezier = obj->fill_Bezier;
    copy->fill_G1 = obj->fill_G1;
    copy->fill_G2 = obj->fill_G2;
    copy->fill_G1Q2 = obj->fill_G1Q2;
    copy->sl_g1q2param = obj->sl_g1q2param;
    copy->g1q2param = obj->g1q2param;
    copy->current_vertex[0] = copy->current_vertex[1] =
      copy->current_edge[0] = copy->current_edge[1] =
      copy->current_facet[0] = copy->current_facet[1] = -1;
    copy->me.displaylist = glGenLists ( BSM_NDL );
    copy->me.dlistmask = 0;
    copy->bl_constr = obj->bl_constr;
    copy->bl_shape_only = obj->bl_shape_only;
    copy->bsm_bl_C = obj->bsm_bl_C;
    copy->nkn1 = obj->nkn1;
    copy->nkn2 = obj->nkn2;
    copy->nlevels = obj->nlevels;
    copy->nblocks = obj->nblocks;
    copy->maxit = obj->maxit;
    copy->spvlist_ok = false;
    return &copy->me;
  }
  else
    return NULL;
} /*GeomObjectCopyBSplineMesh*/

void GeomObjectDeleteBSplineMesh ( GO_BSplineMesh *obj )
{
  glDeleteLists ( obj->me.displaylist, BSM_NDL );
  if ( obj->meshv )           free ( obj->meshv );
  if ( obj->meshvhei )        free ( obj->meshvhei );
  if ( obj->meshvpc )         free ( obj->meshvpc );
  if ( obj->meshhe )          free ( obj->meshhe );
  if ( obj->meshfac )         free ( obj->meshfac );
  if ( obj->meshfhei )        free ( obj->meshfhei );
  if ( obj->mkcp )            free ( obj->mkcp );
  if ( obj->savedcpoints )    free ( obj->savedcpoints );
  if ( obj->spvlist.spel )    free ( obj->spvlist.spel );
  if ( obj->spvlist.spvert )  free ( obj->spvlist.spvert );
  if ( obj->special_patches ) free ( obj->special_patches );
  free ( obj );
} /*GeomObjectDeleteBSplineMesh*/

