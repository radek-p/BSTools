
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


boolean GeomObjectBSplinePatchInsertKnotU ( GO_BSplinePatch *obj,
                                            double knot )
{
  int    dim, ncp, degu, degv, lknu, lknv;
  double *cp, *knu;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( obj->bsp_type == BSP_TYPE_BLENDING_G1 ||
       obj->bsp_type == BSP_TYPE_BLENDING_G2 )
    return false;

  cp = knu = NULL;
  mkcp = NULL;
  dim = obj->me.cpdimen;
  degu = obj->degree_u;
  degv = obj->degree_v;
  lknu = obj->lastknot_u;
  lknv = obj->lastknot_v;
  ncp = (lknu-degu)*(lknv-degv);
  if ( knot < obj->knots_u[degu] || knot > obj->knots_u[lknu-degu] )
    goto failure;
  knu = malloc ( (lknu+2)*sizeof(double) );
  cp = malloc ( (ncp+lknv-degv)*dim*sizeof(double) );
  mkcp = malloc ( ncp+lknv-degv );
  if ( !knu || !cp || !mkcp )
    goto failure;
  memcpy ( knu, obj->knots_u, (lknu+1)*sizeof(double) );
  memcpy ( cp, obj->cpoints, ncp*dim*sizeof(double) );
  if ( !obj->closed_u )
    mbs_multiKnotInsd ( degu, &lknu, knu, 1, (lknv-degv)*dim,
                        0, 0, cp, knot );
  else
    mbs_multiKnotInsClosedd ( degu, &lknu, knu, 1, (lknv-degv)*dim,
                              0, 0, cp, knot );
  ncp = (lknu-degu)*(lknv-degv);
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  free ( obj->knots_u );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->knots_u = knu;
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  obj->lastknot_u = lknu;
  obj->maxknots_u = lknu+1;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  obj->me.dlistmask = 0;
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  return true;

failure:
  if ( cp ) free ( cp );
  if ( knu ) free ( knu );
  if ( mkcp ) free ( mkcp );
  return false;
} /*GeomObjectBSplinePatchInsertKnotU*/

boolean GeomObjectBSplinePatchRemoveKnotU ( GO_BSplinePatch *obj,
                                            int kni )
{
  int    dim, ncp, degu, degv, lknu, lknv;
  double *cp, *knu;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( obj->bsp_type == BSP_TYPE_BLENDING_G1 ||
       obj->bsp_type == BSP_TYPE_BLENDING_G2 )
    return false;

  cp = knu = NULL;
  mkcp = NULL;
  dim = obj->me.cpdimen;
  degu = obj->degree_u;
  degv = obj->degree_v;
  lknu = obj->lastknot_u;
  lknv = obj->lastknot_v;
  ncp = (lknu-degu)*(lknv-degv);
  if ( kni <= degu || kni >= lknu-degu )
    goto failure;
  knu = malloc ( (lknu+1)*sizeof(double) );
  cp = malloc ( ncp*dim*sizeof(double) );
  mkcp = malloc ( ncp );
  if ( !knu || !cp || !mkcp )
    goto failure;
  memcpy ( knu, obj->knots_u, (lknu+1)*sizeof(double) );
  memcpy ( cp, obj->cpoints, ncp*dim*sizeof(double) );
  if ( !obj->closed_u )
    mbs_multiKnotRemoved ( degu, &lknu, knu, 1, (lknv-degv)*dim,
                           0, 0, cp, kni );
  else
    mbs_multiKnotRemoveClosedd ( degu, &lknu, knu, 1, (lknv-degv)*dim,
                                 0, 0, cp, kni );
  ncp = (lknu-degu)*(lknv-degv);
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  free ( obj->knots_u );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->knots_u = knu;
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  obj->lastknot_u = lknu;
  obj->maxknots_u = lknu+1;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  obj->me.dlistmask = 0;
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  return true;

failure:
  if ( cp ) free ( cp );
  if ( knu ) free ( knu );
  if ( mkcp ) free ( mkcp );
  return false;
} /*GeomObjectBSplinePatchRemoveKnotU*/

boolean GeomObjectBSplinePatchInsertKnotV ( GO_BSplinePatch *obj,
                                            double knot )
{
  int    dim, ncp, degu, degv, lknu, lknv;
  double *cp, *knv;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( obj->bsp_type == BSP_TYPE_BLENDING_G1 ||
       obj->bsp_type == BSP_TYPE_BLENDING_G2 )
    return false;

  cp = knv = NULL;
  mkcp = NULL;
  dim = obj->me.cpdimen;
  degu = obj->degree_u;
  degv = obj->degree_v;
  lknu = obj->lastknot_u;
  lknv = obj->lastknot_v;
  ncp = (lknu-degu)*(lknv-degv);
  if ( knot < obj->knots_v[degv] || knot > obj->knots_v[lknv-degv] )
    goto failure;
  knv = malloc ( (lknv+2)*sizeof(double) );
  cp = malloc ( (ncp+lknu-degu)*dim*sizeof(double) );
  mkcp = malloc ( ncp+lknu-degu );
  if ( !knv || !cp || !mkcp )
    goto failure;
  memcpy ( knv, obj->knots_v, (lknv+1)*sizeof(double) );
  memcpy ( cp, obj->cpoints, ncp*dim*sizeof(double) );
  if ( !obj->closed_v )
    mbs_multiKnotInsd ( degv, &lknv, knv, lknu-degu, dim,
                        (lknv-degv)*dim, (lknv-degv+1)*dim, cp, knot );
  else
    mbs_multiKnotInsClosedd ( degv, &lknv, knv, lknu-degu, dim,
                              (lknv-degv)*dim, (lknv-degv+1)*dim, cp, knot );
  ncp = (lknu-degu)*(lknv-degv);
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  free ( obj->knots_v );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->knots_v = knv;
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  obj->lastknot_v = lknv;
  obj->maxknots_v = lknv+1;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  obj->me.dlistmask = 0;
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  return true;

failure:
  if ( cp ) free ( cp );
  if ( knv ) free ( knv );
  if ( mkcp ) free ( mkcp );
  return false;
} /*GeomObjectBSplinePatchInsertKnotV*/

boolean GeomObjectBSplinePatchRemoveKnotV ( GO_BSplinePatch *obj,
                                            int kni )
{
  int    dim, ncp, degu, degv, lknu, lknv;
  double *cp, *knv;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( obj->bsp_type == BSP_TYPE_BLENDING_G1 ||
       obj->bsp_type == BSP_TYPE_BLENDING_G2 )
    return false;

  cp = knv = NULL;
  mkcp = NULL;
  dim = obj->me.cpdimen;
  degu = obj->degree_u;
  degv = obj->degree_v;
  lknu = obj->lastknot_u;
  lknv = obj->lastknot_v;
  ncp = (lknu-degu)*(lknv-degv);
  if ( kni <= degv || kni >= lknv-degv )
    goto failure;
  knv = malloc ( (lknv+1)*sizeof(double) );
  cp = malloc ( ncp*dim*sizeof(double) );
  mkcp = malloc ( ncp );
  if ( !knv || !cp || !mkcp )
    goto failure;
  memcpy ( knv, obj->knots_v, (lknv+1)*sizeof(double) );
  memcpy ( cp, obj->cpoints, ncp*dim*sizeof(double) );
  if ( !obj->closed_v )
    mbs_multiKnotRemoved ( degv, &lknv, knv, lknu-degu, dim,
                           (lknv-degv)*dim, (lknv-degv-1)*dim, cp, kni );
  else
    mbs_multiKnotRemoveClosedd ( degv, &lknv, knv, lknu-degu, dim,
                                 (lknv-degv)*dim, (lknv-degv-1)*dim, cp, kni );
  ncp = (lknu-degu)*(lknv-degv);
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  free ( obj->knots_v );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->knots_v = knv;
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  obj->lastknot_v = lknv;
  obj->maxknots_v = lknv+1;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  obj->me.dlistmask = 0;
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  return true;

failure:
  if ( cp ) free ( cp );
  if ( knv ) free ( knv );
  if ( mkcp ) free ( mkcp );
  return false;
} /*GeomObjectBSplinePatchRemoveKnotV*/

void GeomObjectBSplinePatchMoveKnotU ( GO_BSplinePatch *obj, int kni )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  obj->me.dlistmask &= ~BSP_DLM_PATCH;
} /*GeomObjectBSplinePatchMoveKnotU*/

void GeomObjectBSplinePatchMoveKnotV ( GO_BSplinePatch *obj, int kni )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  obj->me.dlistmask &= ~BSP_DLM_PATCH;
} /*GeomObjectBSplinePatchMoveKnotV*/

boolean GeomObjectBSplinePatchSetUniformKnotsU ( GO_BSplinePatch *obj )
{
  int    i, lkn;
  double *kn;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  lkn = obj->lastknot_u;
  kn = obj->knots_u;
  for ( i = 0; i <= lkn; i++ )
    kn[i] = (double)i;
  obj->me.dlistmask &= ~BSP_DLM_PATCH;
  return true;
} /*GeomObjectBSplinePatchSetUniformKnotsU*/

boolean GeomObjectBSplinePatchSetUniformKnotsV ( GO_BSplinePatch *obj )
{
  int    i, lkn;
  double *kn;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  lkn = obj->lastknot_v;
  kn = obj->knots_v;
  for ( i = 0; i <= lkn; i++ )
    kn[i] = (double)i;
  obj->me.dlistmask &= ~BSP_DLM_PATCH;
  return true;
} /*GeomObjectBSplinePatchSetUniformKnotsV*/

boolean GeomObjectBSplinePatchSetClosedU ( GO_BSplinePatch *obj, boolean closed )
{
  int    deg, lkn, clcK, cpdim, pitch;
  double *kn, *cp, clcT;
  int    i;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( !closed ) {
    obj->closed_u = false;
    return true;
  }
  deg = obj->degree_u;
  lkn = obj->lastknot_u;
  kn = obj->knots_u;
  if ( lkn <= 3*deg )
    return false;
  clcK = lkn-2*deg;
  clcT = kn[clcK+deg]-kn[deg];
  for ( i = 1; i < 2*deg; i++ )
    kn[i] = 0.5*(kn[i]+kn[i+clcK]-clcT);
  for ( i = 1; i < 2*deg; i++ )
    kn[i+clcK] = kn[i]+clcT;
  kn[0] = kn[1];
  kn[lkn] = kn[lkn-1];
  cpdim = obj->me.cpdimen;
  pitch = cpdim*(obj->lastknot_v-obj->degree_v);
  cp = obj->cpoints;
  pkn_AddMatrixd ( 1, pitch*deg, 0, cp, 0, &cp[pitch*clcK], 0, cp );
  pkn_MultMatrixNumd ( 1, pitch*deg, 0, cp, 0.5, 0, cp );
  memcpy ( &cp[pitch*clcK], cp, pitch*deg*sizeof(double) );
  obj->me.dlistmask = 0;
  obj->closed_u = true;
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  return true;
} /*GeomObjectBSplinePatchSetClosedU*/

boolean GeomObjectBSplinePatchSetClosedV ( GO_BSplinePatch *obj, boolean closed )
{
  int    deg, lkn, clcK, cpdim, pitch, ncol;
  double *kn, *cp, clcT;
  int    i;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( !closed ) {
    obj->closed_v = false;
    return true;
  }
  deg = obj->degree_v;
  lkn = obj->lastknot_v;
  kn = obj->knots_v;
  if ( lkn <= 3*deg )
    return false;
  clcK = lkn-2*deg;
  clcT = kn[clcK+deg]-kn[deg];
  for ( i = 1; i < 2*deg; i++ )
    kn[i] = 0.5*(kn[i]+kn[i+clcK]-clcT);
  for ( i = 1; i < 2*deg; i++ )
    kn[i+clcK] = kn[i]+clcT;
  kn[0] = kn[1];
  kn[lkn] = kn[lkn-1];
  cpdim = obj->me.cpdimen;
  pitch = cpdim*(lkn-deg);
  cp = obj->cpoints;
  ncol = obj->lastknot_u-obj->degree_u;
  pkn_AddMatrixd ( ncol, cpdim*deg, pitch, cp, pitch, &cp[cpdim*clcK], pitch, cp );
  pkn_MultMatrixNumd ( ncol, cpdim*deg, pitch, cp, 0.5, pitch, cp );
  pkv_Moved ( ncol, cpdim*deg, pitch, cpdim*clcK, cp );
  obj->me.dlistmask = 0;
  obj->closed_v = true;
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  return true;
} /*GeomObjectBSplinePatchSetClosedV*/

boolean GeomObjectBSplinePatchFlipUV ( GO_BSplinePatch *obj )
{
  void    *sp, *mcp;
  double  *kn;
  int     ncp, n1, n2, pitch1, pitch2, els, i;
  boolean cl;

  sp = pkv_GetScratchMemTop ();
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    goto failure;
  if ( obj->bsp_type == BSP_TYPE_SPHERICAL )
    goto failure;

  n1 = obj->lastknot_v-obj->degree_v;
  n2 = obj->lastknot_u-obj->degree_u;
  ncp = n1*n2;
  els = obj->me.cpdimen*sizeof(double);
  mcp = pkv_GetScratchMem ( ncp*els );
  if ( !mcp )
    goto failure;
  pkv_TransposeMatrixc ( n2, n1, 1, n1, (char*)obj->mkcp, n2, mcp );
  memcpy ( obj->mkcp, mcp, ncp );
  pitch1 = n1*els;
  pitch2 = n2*els;
  pkv_TransposeMatrixc ( n2, n1, els, pitch1, (char*)obj->cpoints, pitch2, mcp );
  memcpy ( obj->cpoints, mcp, ncp*els );

  kn = obj->knots_u;    obj->knots_u = obj->knots_v;        obj->knots_v = kn;
  i = obj->maxknots_u;  obj->maxknots_u = obj->maxknots_v;  obj->maxknots_v = i;
  i = obj->degree_u;    obj->degree_u = obj->degree_v;      obj->degree_v = i;
  i = obj->lastknot_u;  obj->lastknot_u = obj->lastknot_v;  obj->lastknot_v = i;
  cl = obj->closed_u;   obj->closed_u = obj->closed_v;      obj->closed_v = cl;
  i = obj->dens_u;      obj->dens_u = obj->dens_v;          obj->dens_v = i;

  if ( obj->trimmed )
    GeomObjectBSplinePatchFlipTrimmedDomain ( obj );

  pkv_SetScratchMemTop ( sp );
  GeomObjectBSplinePatchDisplayInfoText ( obj );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplinePatchFlipUV*/

