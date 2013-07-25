
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
#include "g1blendingd.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "render.h"
#include "editor.h"
#include "editor_bsp.h"

/* ////////////////////////////////////////////////////////////////////////// */
void GeomObjectBSplinePatchDestroyBlMat ( GO_BSplinePatch *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  if ( obj->blp_amat ) { free ( obj->blp_amat );  obj->blp_amat = NULL; }
  if ( obj->blp_arow ) { free ( obj->blp_arow );  obj->blp_arow = NULL; }
  if ( obj->blp_prof ) { free ( obj->blp_prof );  obj->blp_prof = NULL; }
  obj->blp_lknu = obj->blp_lknv = obj->blp_n = 0;
} /*GeomObjectBSplinePatchDestroyBlMat*/

boolean GeomObjectBSplinePatchIsClamped ( GO_BSplinePatch *obj )
{
  int    deg, lkn;
  double *knots;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  deg = obj->degree_u;  /* ought to be 2 or 3 */
  lkn = obj->lastknot_u;
  knots = obj->knots_u;
  if ( knots[deg] > knots[1] || knots[lkn-1] > knots[lkn-deg] )
    return false;
  deg = obj->degree_v;  /* ought to be the same as obj->degree_u */
  lkn = obj->lastknot_v;
  knots = obj->knots_v;
  if ( knots[deg] > knots[1] || knots[lkn-1] > knots[lkn-deg] )
    return false;
  return true;
} /*GeomObjectBSplinePatchIsClamped*/

boolean GeomObjectBSplinePatchInitBlG1 ( GO_BSplinePatch *obj )
{
  int     lknu, lknv, ncp, i, j, k;
  point3d *cp;
  double  *knu, *knv, x, y;
  byte    *mkcp;

  lknu = max ( obj->lastknot_u, 7 );
  lknv = max ( obj->lastknot_v, 7 );
  ncp = (lknu-2)*(lknv-2);
  cp = malloc ( ncp*sizeof(point3d) );
  knu = malloc ( (lknu+1)*sizeof(double) );
  knv = malloc ( (lknv+1)*sizeof(double) );
  mkcp = malloc ( ncp );
  if ( !cp || !knu || !knv || !mkcp ) {
    if ( cp )   free ( cp );
    if ( knu )  free ( knu );
    if ( knv )  free ( knv );
    if ( mkcp ) free ( mkcp );
    return false;
  }
  free ( obj->cpoints );
  free ( obj->knots_u );
  free ( obj->knots_v );
  free ( obj->mkcp );
  obj->cpoints = (double*)cp;
  obj->knots_u = knu;
  obj->knots_v = knv;
  obj->mkcp = mkcp;
  obj->lastknot_u = lknu;
  obj->lastknot_v = lknv;
  obj->degree_u = obj->degree_v = 2;
  obj->rational = obj->closed_u = obj->closed_v =
  obj->clamped = obj->nharmonic = false;
  obj->me.cpdimen = obj->me.spdimen = 3;
  for ( i = k = 0;  i < lknu-2;  i++ ) {
    x = (double)i/(double)(lknu-3);
    for ( j = 0;  j < lknv-2;  j++, k++ ) {
      y = (double)j/(double)(lknv-3);
      SetPoint3d ( &cp[k], x, y, 0.0 );
    }
  }
  GeomObjectBSplinePatchSetUniformKnotsU ( obj );
  GeomObjectBSplinePatchSetUniformKnotsV ( obj );
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  return true;
} /*GeomObjectBSplinePatchInitBlG1*/

boolean GeomObjectBSplinePatchAdjustBlG1 ( GO_BSplinePatch *obj )
{
  int    lkn;
  double *knots;

  if ( obj->bsp_type != BSP_TYPE_BLENDING_G1 ) {
    obj->bsp_type = BSP_TYPE_GENERAL;  /* to allow degree change */
    GeomObjectBSplinePatchDestroyBlMat ( obj );
    GeomObjectBSplinePatchSetRational ( obj, false );
    GeomObjectDeleteDependencies ( &obj->me );
    if ( obj->degree_u != 2 )
      GeomObjectBSplinePatchSetDegreeU ( obj, 2 );
    if ( obj->degree_v != 2 )
      GeomObjectBSplinePatchSetDegreeV ( obj, 2 );
    if ( obj->degree_u != 2 || obj->degree_v != 2 ||
         obj->lastknot_u < 7 || obj->lastknot_v < 7 ) {
      if ( !GeomObjectBSplinePatchInitBlG1 ( obj ) ) {
        GeomObjectBSplinePatchAdjustGeneral ( obj );
        return false;
      }
    }
    else {
      obj->clamped = GeomObjectBSplinePatchIsClamped ( obj );
      GeomObjectBSplinePatchSetUniformKnotsU ( obj );
      GeomObjectBSplinePatchSetUniformKnotsV ( obj );
      if ( obj->clamped ) {
        lkn = obj->lastknot_u;
        knots = obj->knots_u;
        knots[0] = knots[1] = knots[2];
        knots[lkn] = knots[lkn-1] = knots[lkn-2];
        lkn = obj->lastknot_v;
        knots = obj->knots_v;
        knots[0] = knots[1] = knots[2];
        knots[lkn] = knots[lkn-1] = knots[lkn-2];
      }
    }
    obj->blp_range[0] = obj->blp_range[2] = 2;
    obj->blp_range[1] = obj->lastknot_u-5;
    obj->blp_range[3] = obj->lastknot_v-5;
    GeomObjectBSplinePatchMarkBLRange ( obj );
    obj->closed_v = obj->nharmonic = false;
    obj->bsp_type = BSP_TYPE_BLENDING_G1;
  }
  return true;
} /*GeomObjectBSplinePatchAdjustBlG1*/

boolean GeomObjectBSplinePatchInitBlG2 ( GO_BSplinePatch *obj )
{
  int     lknu, lknv, ncp, i, j, k;
  point3d *cp;
  double  *knu, *knv, x, y;
  byte    *mkcp;

  lknu = max ( obj->lastknot_u, 10 );
  lknv = max ( obj->lastknot_v, 10 );
  ncp = (lknu-3)*(lknv-3);
  cp = malloc ( ncp*sizeof(point3d) );
  knu = malloc ( (lknu+1)*sizeof(double) );
  knv = malloc ( (lknv+1)*sizeof(double) );
  mkcp = malloc ( ncp );
  if ( !cp || !knu || !knv || !mkcp ) {
    if ( cp )   free ( cp );
    if ( knu )  free ( knu );
    if ( knv )  free ( knv );
    if ( mkcp ) free ( mkcp );
    return false;
  }
  free ( obj->cpoints );
  free ( obj->knots_u );
  free ( obj->knots_v );
  free ( obj->mkcp );
  obj->cpoints = (double*)cp;
  obj->knots_u = knu;
  obj->knots_v = knv;
  obj->mkcp = mkcp;
  obj->lastknot_u = lknu;
  obj->lastknot_v = lknv;
  obj->degree_u = obj->degree_v = 3;
  obj->rational = obj->closed_u = obj->closed_v =
  obj->clamped = obj->nharmonic = false;
  obj->me.cpdimen = obj->me.spdimen = 3;
  for ( i = k = 0;  i < lknu-3;  i++ ) {
    x = (double)i/(double)(lknu-4);
    for ( j = 0;  j < lknv-3;  j++, k++ ) {
      y = (double)j/(double)(lknv-4);
      SetPoint3d ( &cp[k], x, y, 0.0 );
    }
  }
  GeomObjectBSplinePatchSetUniformKnotsU ( obj );
  GeomObjectBSplinePatchSetUniformKnotsV ( obj );
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  return true;
} /*GeomObjectBSplinePatchInitBlG2*/

boolean GeomObjectBSplinePatchAdjustBlG2 ( GO_BSplinePatch *obj )
{
  int    lkn;
  double *knots;

  if ( obj->bsp_type != BSP_TYPE_BLENDING_G2 ) {
    obj->bsp_type = BSP_TYPE_GENERAL;  /* to allow degree change */
    GeomObjectBSplinePatchDestroyBlMat ( obj );
    GeomObjectBSplinePatchSetRational ( obj, false );
    GeomObjectDeleteDependencies ( &obj->me );
    if ( obj->degree_u != 3 )
      GeomObjectBSplinePatchSetDegreeU ( obj, 3 );
    if ( obj->degree_v != 3 )
      GeomObjectBSplinePatchSetDegreeV ( obj, 3 );
    if ( obj->degree_u != 3 || obj->degree_v != 3 ||
         obj->lastknot_u < 10 || obj->lastknot_v < 10 ) {
      if ( !GeomObjectBSplinePatchInitBlG2 ( obj ) ) {
        GeomObjectBSplinePatchAdjustGeneral ( obj );
        return false;
      }
    }
    else {
      obj->clamped = GeomObjectBSplinePatchIsClamped ( obj );
      GeomObjectBSplinePatchSetUniformKnotsU ( obj );
      GeomObjectBSplinePatchSetUniformKnotsV ( obj );
      if ( obj->clamped ) {
        lkn = obj->lastknot_u;
        knots = obj->knots_u;
        knots[0] = knots[1] = knots[2] = knots[3];
        knots[lkn] = knots[lkn-1] = knots[lkn-2] = knots[lkn-3];
        lkn = obj->lastknot_v;
        knots = obj->knots_v;
        knots[0] = knots[1] = knots[2] = knots[3];
        knots[lkn] = knots[lkn-1] = knots[lkn-2] = knots[lkn-3];
      }
    }
    obj->blp_range[0] = obj->blp_range[2] = 3;
    obj->blp_range[1] = obj->lastknot_u-7;
    obj->blp_range[3] = obj->lastknot_v-7;
    GeomObjectBSplinePatchMarkBLRange ( obj );
    obj->closed_v = obj->nharmonic = false;
    obj->bsp_type = BSP_TYPE_BLENDING_G2;
  }
  return true;
} /*GeomObjectBSplinePatchAdjustBlG2*/

void GeomObjectBSplinePatchSetEntireBlendingRange ( GO_BSplinePatch *obj )
{
  obj->blp_range[0] = obj->blp_range[2] = obj->degree_u;
  obj->blp_range[1] = obj->lastknot_u-2*obj->degree_u-1;
  obj->blp_range[3] = obj->lastknot_v-2*obj->degree_v-1;
  GeomObjectBSplinePatchMarkBLRange ( obj );
} /*GeomObjectBSplinePatchSetEntireBlendingRange*/

void GeomObjectBSplinePatchFreeToClamped ( GO_BSplinePatch *obj )
{
  void   *sp;
  double newkn[4];
  int    i, dim, degu, degv, lknu, lknv, pitch;

  sp = pkv_GetScratchMemTop ();
  degu = obj->degree_u;
  degv = obj->degree_v;
  if ( degu < 2 || degu > 3 || degu != degv )
    goto way_out;
  lknu = obj->lastknot_u;
  lknv = obj->lastknot_v;
  dim = obj->me.cpdimen;
  pitch = dim*(lknv-degv);
  newkn[0] = 0.0;
  for ( i = 1; i <= degv; i++ )
    newkn[i] = (double)degv;
  mbs_multiBSChangeLeftKnotsd ( lknu-degu, dim, degv, obj->knots_v,
                                pitch, obj->cpoints, newkn );
  for ( i = 0; i < degv; i++ )
    newkn[i] = (double)(lknv-degv);
  newkn[degv] = (double)lknv;
  mbs_multiBSChangeRightKnotsd ( lknu-degu, dim, degv, lknv, obj->knots_v,
                                 pitch, obj->cpoints, newkn );
  if ( !obj->closed_u ) {
    newkn[0] = 0.0;
    for ( i = 1; i <= degu; i++ )
      newkn[i] = (double)degu;
    mbs_multiBSChangeLeftKnotsd ( 1, pitch, degu, obj->knots_u,
                                  0, obj->cpoints, newkn );
    for ( i = 0; i < degu; i++ )
      newkn[i] = (double)(lknu-degu);
    newkn[degu] = (double)lknu;
    mbs_multiBSChangeRightKnotsd ( 1, pitch, degu, lknu, obj->knots_u,
                                   0, obj->cpoints, newkn );
  }
  obj->clamped = true;
  obj->me.dlistmask &= !BSP_DLM_CNET;

way_out:
  pkv_SetScratchMemTop ( sp );
} /*GeomObjectBSplinePatchFreeToClamped*/

void GeomObjectBSplinePatchClampedToFree ( GO_BSplinePatch *obj )
{
  void   *sp;
  double newkn[4];
  int    i, dim, degu, degv, lknu, lknv, pitch;

  sp = pkv_GetScratchMemTop ();
  degu = obj->degree_u;
  degv = obj->degree_v;
  if ( degu < 2 || degu > 3 || degu != degv )
    goto way_out;
  lknu = obj->lastknot_u;
  lknv = obj->lastknot_v;
  dim = obj->me.cpdimen;
  pitch = dim*(lknv-degv);
  if ( !obj->closed_u ) {
    for ( i = 0; i <= degv; i++ )
      newkn[i] = (double)i;
    mbs_multiBSChangeLeftKnotsd ( lknu-degu, dim, degv, obj->knots_v,
                                  pitch, obj->cpoints, newkn );
    for ( i = 0; i <= degv; i++ )
      newkn[i] = (double)(lknv-degv+i);
    mbs_multiBSChangeRightKnotsd ( lknu-degu, dim, degv, lknv, obj->knots_v,
                                   pitch, obj->cpoints, newkn );
  }
  for ( i = 0; i <= degu; i++ )
    newkn[i] = (double)i;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch, degu, obj->knots_u,
                                0, obj->cpoints, newkn );
  for ( i = 0; i <= degu; i++ )
    newkn[i] = (double)(lknu-degu+i);
  mbs_multiBSChangeRightKnotsd ( 1, pitch, degu, lknu, obj->knots_u,
                                 0, obj->cpoints, newkn );
  obj->clamped = false;
  obj->me.dlistmask &= !BSP_DLM_CNET;

way_out:
  pkv_SetScratchMemTop ( sp );
} /*GeomObjectBSplinePatchClampedToFree*/

void GeomObjectBSPatchMarkCPGeneral ( GO_BSplinePatch *obj )
{
  int  k, l, ncp;
  byte *mkcp;

  mkcp = obj->mkcp;
  k = obj->lastknot_u-obj->degree_u;
  l = obj->lastknot_v-obj->degree_v;
  ncp = k*l;
  for ( k = 0; k < ncp; k++ )
    mkcp[k] |= MASK_CP_BOUNDARY;
} /*GeomObjectBSPatchMarkCPGeneral*/

void GeomObjectBSplinePatchMarkBLRange ( GO_BSplinePatch *obj )
{
  int  i, j, k, l, ncp, *range;
  byte *mkcp;

  mkcp = obj->mkcp;
  k = obj->lastknot_u-obj->degree_u;
  l = obj->lastknot_v-obj->degree_v;
  ncp = k*l;
  for ( i = 0; i < ncp; i++ )
    mkcp[i] |= MASK_CP_BOUNDARY;
  range = obj->blp_range;
  for ( i = range[0]; i <= range[1]; i++ )
    for ( j = range[2]; j <= range[3]; j++ ) {
      k = i*l+j;
      if ( k >= 0 && k < ncp )
        mkcp[k] &= ~MASK_CP_BOUNDARY;
    }
  obj->me.dlistmask &= ~BSP_DLM_CNET;
} /*GeomObjectBSplinePatchMarkBLRange*/

boolean GeomObjectBSplinePatchRefine ( GO_BSplinePatch *obj )
{
  void    *sp;
  int     degu, degv, lknu, lknv, nlknu, nlknv, ncp, cpdimen, pitch1, pitch2;
  double  *newknotsu, *newknotsv, *newcpoints, *auxcp;
  byte    *mkcp;
  boolean clamped;
  int     i;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  sp = pkv_GetScratchMemTop ();
  degu = obj->degree_u;
  degv = obj->degree_v;
  lknu = obj->lastknot_u;
  lknv = obj->lastknot_v;
  cpdimen = obj->me.cpdimen;
  nlknu = 2*(lknu-degu);
  nlknv = 2*(lknv-degv);
  ncp = (nlknu-degu)*(nlknv-degv);
  pitch1 = (lknv-degv)*cpdimen;
  pitch2 = (nlknv-degv)*cpdimen;
  newknotsu = newknotsv = newcpoints = NULL;
  mkcp = NULL;
  if ( nlknu >= MAX_BSP_KNOTS || nlknv >= MAX_BSP_KNOTS )
    goto failure;
  newknotsu = malloc ( (nlknu+1)*sizeof(double) );
  newknotsv = malloc ( (nlknv+1)*sizeof(double) );
  newcpoints = malloc ( ncp*cpdimen*sizeof(double) );
  mkcp = malloc ( ncp );
  auxcp = pkv_GetScratchMemd ( (nlknu-degu)*pitch1 );
  if ( !newknotsu || !newknotsv || !newcpoints || !mkcp || !auxcp )
    goto failure;

  clamped = obj->clamped;
  if ( clamped )
    GeomObjectBSplinePatchClampedToFree ( obj );
  if ( !mbs_multiLaneRiesenfeldd ( pitch1, 1, degu, lknu,
                                   0, obj->cpoints, &nlknu, 0, auxcp ) )
    goto failure;
  if ( !mbs_multiLaneRiesenfeldd ( cpdimen, nlknu-degu, degv, lknv,
                                   pitch1, auxcp, &nlknv, pitch2, newcpoints ) )
    goto failure;
  for ( i = 0; i <= nlknu; i++ )
    newknotsu[i] = (double)i;
  for ( i = 0; i <= nlknv; i++ )
    newknotsv[i] = (double)i;
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  GeomObjectAssignBSPatch ( obj, obj->me.spdimen, obj->rational,
                            degu, nlknu, newknotsu, degv, nlknv, newknotsv,
                            newcpoints, mkcp, obj->closed_u, obj->closed_v );
  obj->blp_range[0] = 2*obj->blp_range[0]-degu;
  obj->blp_range[1] = 2*obj->blp_range[1]+1;
  obj->blp_range[2] = 2*obj->blp_range[2]-degu;
  obj->blp_range[3] = 2*obj->blp_range[3]+1;
  GeomObjectBSplinePatchMarkBLRange ( obj );
  if ( clamped )
    GeomObjectBSplinePatchFreeToClamped ( obj );
  obj->nharmonic = false;
  pkv_SetScratchMemTop ( sp );
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  return true;

failure:
  if ( newknotsu ) free ( newknotsu );
  if ( newknotsv ) free ( newknotsv );
  if ( newcpoints ) free ( newcpoints );
  if ( mkcp ) free ( mkcp );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplinePatchRefine*/

boolean GeomObjectBSplinePatchSetupNHBlMat ( GO_BSplinePatch *obj )
{
  void  *sp;
  int    bl_lknu, bl_lknv;

  sp = pkv_GetScratchMemTop ();
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    goto failure;

  GeomObjectBSplinePatchDestroyBlMat ( obj );
  switch ( obj->bsp_type ) {
case BSP_TYPE_BLENDING_G1:
    obj->blp_deg = 2;
    if ( obj->closed_u ) {
/* for now only the entire range of "u" parameter */
      bl_lknu = obj->lastknot_u;
      bl_lknv = obj->blp_range[3]-obj->blp_range[2]+7;
      if ( bl_lknu < 7 || bl_lknv < 7 )
        goto failure;
      if ( g1bl_SetupClosedBiharmAMatrixd ( bl_lknu, bl_lknv,
                 &obj->blp_n, &obj->blp_prof, &obj->blp_amat, &obj->blp_arow ) ) {
        if ( !pkn_NRBSymCholeskyDecompd ( obj->blp_n, obj->blp_prof,
                                          obj->blp_amat, obj->blp_arow, NULL ) )
          goto failure;
      }
      else
        goto failure;
    }
    else {
      bl_lknu = obj->blp_range[1]-obj->blp_range[0]+7;
      bl_lknv = obj->blp_range[3]-obj->blp_range[2]+7;
      if ( bl_lknu < 7 || bl_lknv < 7 )
        goto failure;
      if ( g1bl_SetupBiharmAMatrixd ( bl_lknu, bl_lknv,
                 &obj->blp_n, &obj->blp_prof, &obj->blp_amat, &obj->blp_arow ) ) {
        if ( !pkn_NRBSymCholeskyDecompd ( obj->blp_n, obj->blp_prof,
                                          obj->blp_amat, obj->blp_arow, NULL ) )
          goto failure;
      }
      else
        goto failure;
    }
    break;
case BSP_TYPE_BLENDING_G2:
    obj->blp_deg = 3;
    if ( obj->closed_u ) {
/* for now only the entire range of "u" parameter */
      bl_lknu = obj->lastknot_u;
      bl_lknv = obj->blp_range[3]-obj->blp_range[2]+10;
      if ( bl_lknu < 10 || bl_lknv < 10 )
        goto failure;
      if ( g2bl_SetupClosedTriharmAMatrixd ( bl_lknu, bl_lknv,
                 &obj->blp_n, &obj->blp_prof, &obj->blp_amat, &obj->blp_arow ) ) {
        if ( !pkn_NRBSymCholeskyDecompd ( obj->blp_n, obj->blp_prof,
                                          obj->blp_amat, obj->blp_arow, NULL ) )
          goto failure;
      }
      else
        goto failure;
    }
    else {
      bl_lknu = obj->blp_range[1]-obj->blp_range[0]+10;
      bl_lknv = obj->blp_range[3]-obj->blp_range[2]+10;
      if ( bl_lknu < 10 || bl_lknv < 10 )
        goto failure;
      if ( g2bl_SetupTriharmAMatrixd ( bl_lknu, bl_lknv,
                 &obj->blp_n, &obj->blp_prof, &obj->blp_amat, &obj->blp_arow ) ) {
        if ( !pkn_NRBSymCholeskyDecompd ( obj->blp_n, obj->blp_prof,
                                          obj->blp_amat, obj->blp_arow, NULL ) )
          goto failure;
      }
      else
        goto failure;
    }
    break;
default:
    goto failure;
  }
  obj->blp_closed_u = obj->closed_u;
  obj->blp_lknu = bl_lknu;
  obj->blp_lknv = bl_lknv;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  GeomObjectBSplinePatchDestroyBlMat ( obj );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplinePatchSetupNHBlMat*/

boolean GeomObjectBSplinePatchUpdNHarmonic ( GO_BSplinePatch *obj )
{
  void    *sp;
  double  *rhs, *cp;
  int     blp_n, bl_lknu, bl_lknv, pitch, fcp;
  boolean clamped;

  sp = pkv_GetScratchMemTop ();
  clamped = false;
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    goto failure;
  if ( obj->clamped ) {
    GeomObjectBSplinePatchClampedToFree ( obj );
    clamped = true;
  }
  pitch = 3*(obj->lastknot_v-obj->degree_v);
  switch ( obj->bsp_type ) {
case BSP_TYPE_BLENDING_G1:
    if ( obj->closed_u ) {
      bl_lknu = obj->lastknot_u;
      bl_lknv = obj->blp_range[3]-obj->blp_range[2]+7;
      if ( !obj->blp_amat ||
           obj->blp_lknu != bl_lknu || obj->blp_lknv != bl_lknv ||
           obj->blp_deg != 2 || !obj->closed_u ) {
        if ( !GeomObjectBSplinePatchSetupNHBlMat ( obj ) )
          goto failure;
      }
      blp_n = obj->blp_n;
      rhs = pkv_GetScratchMemd ( 3*blp_n );
      if ( !rhs )
        goto failure;
      fcp = obj->blp_range[2]-2;
      cp = &obj->cpoints[3*fcp];
      if ( !g1bl_SetupClosedBiharmRHSd ( bl_lknu, bl_lknv, 3, pitch, cp, rhs ) )
        goto failure;
      if ( !pkn_NRBLowerTrSolved ( blp_n, obj->blp_prof, obj->blp_amat,
                                   obj->blp_arow, 3, 3, rhs, 3, rhs ) )
        goto failure;
      if ( !pkn_NRBUpperTrSolved ( blp_n, obj->blp_prof, obj->blp_amat,
                                   obj->blp_arow, 3, 3, rhs, 3, rhs ) )
        goto failure;
      pkv_Selectd ( obj->lastknot_u-4, 3*(bl_lknv-6), 3*(bl_lknv-6), pitch,
                    rhs, &cp[6] );
      pkv_Selectd ( 2, pitch, pitch, pitch, obj->cpoints,
                    &obj->cpoints[(obj->lastknot_u-4)*pitch] );
    }
    else {
      bl_lknu = obj->blp_range[1]-obj->blp_range[0]+7;
      bl_lknv = obj->blp_range[3]-obj->blp_range[2]+7;
      if ( !obj->blp_amat ||
           obj->blp_lknu != bl_lknu || obj->blp_lknv != bl_lknv ||
           obj->blp_deg != 2 || obj->closed_u ) {
        if ( !GeomObjectBSplinePatchSetupNHBlMat ( obj ) )
          goto failure;
      }
      blp_n = obj->blp_n;
      rhs = pkv_GetScratchMemd ( 3*blp_n );
      if ( !rhs )
        goto failure;
      fcp = (obj->blp_range[0]-2)*(obj->lastknot_v-2)+(obj->blp_range[2]-2);
      cp = &obj->cpoints[3*fcp];
      if ( !g1bl_SetupBiharmRHSd ( bl_lknu, bl_lknv, 3, pitch, cp, rhs ) )
        goto failure;
      if ( !pkn_NRBLowerTrSolved ( blp_n, obj->blp_prof, obj->blp_amat,
                                   obj->blp_arow, 3, 3, rhs, 3, rhs ) )
        goto failure;
      if ( !pkn_NRBUpperTrSolved ( blp_n, obj->blp_prof, obj->blp_amat,
                                   obj->blp_arow, 3, 3, rhs, 3, rhs ) )
        goto failure;
      pkv_Selectd ( bl_lknu-6, 3*(bl_lknv-6), 3*(bl_lknv-6), pitch,
                    rhs, &cp[2*pitch+6] );
    }
    break;
case BSP_TYPE_BLENDING_G2:
    if ( obj->closed_u ) {
      bl_lknu = obj->lastknot_u;
      bl_lknv = obj->blp_range[3]-obj->blp_range[2]+10;
      if ( !obj->blp_amat ||
           obj->blp_lknu != bl_lknu || obj->blp_lknv != bl_lknv ||
           obj->blp_deg != 2 || !obj->closed_u ) {
        if ( !GeomObjectBSplinePatchSetupNHBlMat ( obj ) )
          goto failure;
      }
      blp_n = obj->blp_n;
      rhs = pkv_GetScratchMemd ( 3*blp_n );
      if ( !rhs )
        goto failure;
      fcp = obj->blp_range[2]-3;
      cp = &obj->cpoints[3*fcp];
      if ( !g2bl_SetupClosedTriharmRHSd ( bl_lknu, bl_lknv, 3, pitch, cp, rhs ) )
        goto failure;
      if ( !pkn_NRBLowerTrSolved ( blp_n, obj->blp_prof, obj->blp_amat,
                                   obj->blp_arow, 3, 3, rhs, 3, rhs ) )
        goto failure;
      if ( !pkn_NRBUpperTrSolved ( blp_n, obj->blp_prof, obj->blp_amat,
                                   obj->blp_arow, 3, 3, rhs, 3, rhs ) )
        goto failure;
      pkv_Selectd ( obj->lastknot_u-6, 3*(bl_lknv-9), 3*(bl_lknv-9), pitch,
                    rhs, &cp[9] );
      pkv_Selectd ( 3, pitch, pitch, pitch, obj->cpoints,
                    &obj->cpoints[(obj->lastknot_u-6)*pitch] );
    }
    else {
      bl_lknu = obj->blp_range[1]-obj->blp_range[0]+10;
      bl_lknv = obj->blp_range[3]-obj->blp_range[2]+10;
      if ( !obj->blp_amat ||
           obj->blp_lknu != bl_lknu || obj->blp_lknv != bl_lknv ||
           obj->blp_deg != 2 || obj->closed_u ) {
        if ( !GeomObjectBSplinePatchSetupNHBlMat ( obj ) )
          goto failure;
      }
      blp_n = obj->blp_n;
      rhs = pkv_GetScratchMemd ( 3*blp_n );
      if ( !rhs )
        goto failure;
      fcp = (obj->blp_range[0]-3)*(obj->lastknot_v-3)+(obj->blp_range[2]-3);
      cp = &obj->cpoints[3*fcp];
      if ( !g2bl_SetupTriharmRHSd ( bl_lknu, bl_lknv, 3, pitch, cp, rhs ) )
        goto failure;
      if ( !pkn_NRBLowerTrSolved ( blp_n, obj->blp_prof, obj->blp_amat,
                                   obj->blp_arow, 3, 3, rhs, 3, rhs ) )
        goto failure;
      if ( !pkn_NRBUpperTrSolved ( blp_n, obj->blp_prof, obj->blp_amat,
                                   obj->blp_arow, 3, 3, rhs, 3, rhs ) )
        goto failure;
      pkv_Selectd ( bl_lknu-9, 3*(bl_lknv-9), 3*(bl_lknv-9), pitch,
                    rhs, &cp[3*pitch+9] );
    }
    break;
default:
    goto failure;
  }
  if ( clamped )  /* restore clamped boundary */
    GeomObjectBSplinePatchFreeToClamped ( obj );
  obj->me.dlistmask = 0;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( clamped )
    GeomObjectBSplinePatchFreeToClamped ( obj );
  obj->nharmonic = false;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplinePatchUpdNHarmonic*/

boolean GeomObjectBSplinePatchSetNHarmonic ( GO_BSplinePatch *obj,
                                             boolean nharmonic )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  switch ( obj->bsp_type ) {
case BSP_TYPE_BLENDING_G1:
case BSP_TYPE_BLENDING_G2:
    obj->nharmonic = nharmonic;
    if ( nharmonic )
      return GeomObjectBSplinePatchUpdNHarmonic ( obj );
    else
      return false;  /* no change */
default:
    return false;
  }
} /*GeomObjectBSplinePatchSetNHarmonic*/

/* ////////////////////////////////////////////////////////////////////////// */
void GeomObjectBSplinePatchDisplayInfoText ( GO_BSplinePatch *obj )
{
  char s[160];

  sprintf ( s, "deg = (%d,%d), lkn = (%d,%d)", obj->degree_u, obj->degree_v,
            obj->lastknot_u, obj->lastknot_v );
  SetStatusText ( s, true );
} /*GeomObjectBSplinePatchDisplayInfoText*/

