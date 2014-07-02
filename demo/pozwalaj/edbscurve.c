
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
#include "mengerc.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "render.h"
#include "editor.h"
#include "editor_bsc.h"

void GeomObjectBSplineCurveInitMC ( GO_BSplineCurve *obj )
{
  static const double defpparam[5] =
    {1.420666814e+06, 1.106542879e+07, 1.063609655e+04,
     1.219877706e+05, 5.771308492e+03};

  obj->mc_exponent = 4.0;
  memcpy ( obj->mc_pparam, defpparam, 5*sizeof(double) );
  obj->mc_nqkn = 3;
  obj->mc_ppopt = 3;
} /*GeomObjectBSplineCurveInitMC*/

void GeomObjectBSplineCurveCopyMC ( GO_BSplineCurve *obj, GO_BSplineCurve *copy )
{
  copy->mc_exponent = obj->mc_exponent;
  memcpy ( copy->mc_pparam, obj->mc_pparam, 5*sizeof(double) );
  copy->mc_nqkn = obj->mc_nqkn;
  copy->mc_ppopt = obj->mc_ppopt;
} /*GeomObjectBSplineCurveCopyMC*/

boolean GeomObjectInitBSplineCurve ( GO_BSplineCurve *obj,
                                     char spdimen, boolean rational )
{
  static const double inikn[4] = { 0.0,0.0,1.0,1.0 };
  static const double inicp[8] = {-1.0,0.0,0.0,1.0, 1.0,0.0,0.0,1.0};

  obj->me.obj_type = GO_BSPLINE_CURVE;
  obj->me.ident = -1;
  obj->me.spdimen = spdimen;
  obj->rational = rational;
  if ( rational ) {
    obj->me.cpdimen = spdimen+1;
    obj->weightpoints = malloc ( obj->me.cpdimen*sizeof(double) );
    if ( !obj->weightpoints )
      return false;
  }
  else {
    obj->me.cpdimen = spdimen;
    obj->weightpoints = NULL;
  }
  obj->mengerc = false;
  obj->me.active = false;
  obj->me.name[0] = 0;
  obj->maxknots = 4;
  obj->savedsize = 0;
  obj->knots = malloc ( 4*sizeof(double) );
  obj->cpoints = malloc ( 2*obj->me.cpdimen*sizeof(double) );
  obj->mkcp = malloc ( 2 );
  if ( !obj->knots || !obj->cpoints || !obj->mkcp ) {
    if ( obj->knots )        free ( obj->knots );
    if ( obj->cpoints )      free ( obj->cpoints );
    if ( obj->weightpoints ) free ( obj->weightpoints );
    if ( obj->mkcp )         free ( obj->mkcp );
    return false;
  }
  memcpy ( obj->knots, inikn, 4*sizeof(double) );
  GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                             2, inicp, obj->cpoints );
  if ( obj->rational )
    GeomObjectSetupWeightPoints ( obj->me.cpdimen, 2, obj->cpoints,
                                  obj->weightpoints );
  memset ( obj->mkcp, MASK_CP_MOVEABLE, 2 );
  obj->lastknot = 3;
  obj->degree = 1;
  obj->closed = false;
  obj->view_curve = obj->view_cpoly = true;
  obj->view_bpoly = false;
  obj->view_curvature = obj->view_torsion = false;
  obj->graph_dens = 10;
  obj->me.displaylist = glGenLists ( BSC_NDL );
  obj->me.dlistmask = 0;
  obj->me.colour[0] = obj->me.colour[1] = 1.0;
  obj->me.colour[2] = 0.0;
  obj->me.display_pretrans = false;
  IdentTrans3d ( &obj->me.pretrans );
  obj->pipe_diameter = 0.1;
  GeomObjectBSplineCurveInitMC ( obj );
  return true;
} /*GeomObjectInitBSplineCurve*/

geom_object *GeomObjectAddBSplineCurve ( const char *name,
                                         char spdimen, boolean rational )
{
  GO_BSplineCurve *obj;

  obj = malloc ( sizeof(GO_BSplineCurve) );
  if ( obj ) {
    memset ( obj, 0, sizeof(GO_BSplineCurve) );
    if ( !GeomObjectInitBSplineCurve ( obj, spdimen, rational ) ) {
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
} /*GeomObjectAddBSplineCurve*/

geom_object *GeomObjectCopyBSplineCurve ( GO_BSplineCurve *obj )
{
  GO_BSplineCurve *copy;
  int             ncp;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return NULL;
  copy = malloc ( sizeof(GO_BSplineCurve) );
  if ( copy ) {
    ncp = obj->lastknot - obj->degree;
    memset ( copy, 0, sizeof(GO_BSplineCurve) );
    copy->me.obj_type = GO_BSPLINE_CURVE;
    copy->me.ident = -1;
    copy->me.spdimen = obj->me.spdimen;
    copy->me.cpdimen = obj->me.cpdimen;
    copy->me.active = false;
    strcpy ( copy->me.name, obj->me.name );
    copy->me.display_pretrans = false;
    copy->me.pretrans = obj->me.pretrans;
    copy->maxknots = obj->maxknots;
    if ( obj->rational ) {
      copy->weightpoints = malloc ( (ncp-1)*obj->me.cpdimen*sizeof(double) );
      if ( !copy->weightpoints ) {
        free ( copy );
        return NULL;
      }
    }
    else
      copy->weightpoints = NULL;
    copy->knots = malloc ( obj->maxknots*sizeof(double) );
    copy->cpoints = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
    copy->mkcp = malloc ( ncp );
    if ( !copy->knots || !copy->cpoints || !copy->mkcp ) {
      if ( copy->knots )        free ( copy->knots );
      if ( copy->cpoints )      free ( copy->cpoints );
      if ( copy->weightpoints ) free ( copy->weightpoints );
      if ( copy->mkcp )         free ( copy->mkcp );
      free ( copy );
      return NULL;
    }
    copy->lastknot = obj->lastknot;
    memcpy ( copy->knots, obj->knots, (obj->lastknot+1)*sizeof(double) );
    copy->degree = obj->degree;
    memcpy ( copy->cpoints, obj->cpoints,
             ncp*obj->me.cpdimen*sizeof(double) );
    memset ( copy->mkcp, MASK_CP_MOVEABLE, ncp );
    copy->rational = obj->rational;
    if ( copy->rational )
      GeomObjectSetupWeightPoints ( copy->me.cpdimen, obj->lastknot-obj->degree,
                                    copy->cpoints, copy->weightpoints );
    copy->closed = obj->closed;
    copy->view_curve = copy->view_cpoly = true;
    copy->view_bpoly = false;
    copy->view_curvature = copy->view_torsion = false;
    copy->graph_dens = 10;
    copy->pipe_diameter = obj->pipe_diameter;
    copy->me.displaylist = glGenLists ( BSC_NDL );
    copy->me.dlistmask = 0;
    GeomObjectBSplineCurveCopyMC ( obj, copy );
    return &copy->me;
  }
  else
    return NULL;
} /*GeomObjectCopyBSplineCurve*/

void GeomObjectDeleteBSplineCurve ( GO_BSplineCurve *obj )
{
  glDeleteLists ( obj->me.displaylist, BSC_NDL );
  if ( obj->knots )        free ( obj->knots );
  if ( obj->cpoints )      free ( obj->cpoints );
  if ( obj->savedcpoints ) free ( obj->savedcpoints );
  if ( obj->weightpoints ) free ( obj->weightpoints );
  if ( obj->mkcp )         free ( obj->mkcp );
  free ( obj );
} /*GeomObjectDeleteBSplineCurve*/

void GeomObjectAssignBSCurve ( GO_BSplineCurve *obj, int spdimen, boolean rational,  
                               int degree, int lastknot, double *knots,  
                               double *cpoints, byte *mkcp, boolean closed )
{
  if ( obj->cpoints ) free ( obj->cpoints );
  if ( obj->knots )   free ( obj->knots );
  if ( obj->mkcp )    free ( obj->mkcp );
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  obj->savedsize = 0;
  obj->me.spdimen = obj->me.cpdimen = spdimen;
  obj->rational = rational;
  if ( rational ) obj->me.cpdimen ++;
  obj->cpoints = cpoints;
  obj->knots = knots;
  obj->mkcp = mkcp;
  obj->degree = degree;
  obj->lastknot = lastknot;
  obj->closed = closed;
  obj->me.dlistmask = 0;
} /*GeomObjectAssignBSCurve*/

boolean GeomObjectBSplineCurveSetRational ( GO_BSplineCurve *obj )
{
  double *cp, *wp;
  int    i, j, ncp, dim, cdim;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  if ( !obj->rational ) {
    dim = obj->me.cpdimen;
    cdim = dim+1;
    ncp = obj->lastknot - obj->degree;
    cp = malloc ( ncp*cdim*sizeof(double) );
    wp = malloc ( (ncp-1)*cdim*sizeof(double) );
    if ( !cp || !wp ) {
      if ( cp ) free ( cp );
      if ( wp ) free ( wp );
      return false;
    }
    pkv_Selectd ( ncp, dim, dim, cdim, obj->cpoints, cp );
    for ( i = 0, j = dim;  i < ncp;  i++, j += cdim )
      cp[j] = 1.0;
    GeomObjectSetupWeightPoints ( cdim, ncp, cp, wp );
    free ( obj->cpoints );
    if ( obj->savedcpoints )
      { free ( obj->savedcpoints );  obj->savedcpoints = NULL; }
    obj->savedsize = 0;
    obj->cpoints = cp;
    obj->weightpoints = wp;
    obj->me.cpdimen = cdim;
    obj->me.dlistmask = 0;
    obj->rational = true;
    GeomObjectProcessDependencies ( (geom_object*)obj );
  }
  return true;
} /*GeomObjectBSplineCurveSetRational*/

boolean GeomObjectBSplineCurveSetNonRational ( GO_BSplineCurve *obj )
{
  double *cp, *rcp, w;
  int    i, j, k, l, ncp, dim, cdim;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  if ( obj->rational ) {
    dim = obj->me.cpdimen;
    cdim = dim-1;
    ncp = obj->lastknot - obj->degree;
    cp = malloc ( ncp*cdim*sizeof(double) );
    if ( !cp )
      return false;
    rcp = obj->cpoints;
    for ( i = j = k = 0;  i < ncp;  i++, j += dim, k += cdim ) {
      w = rcp[j+cdim];
      for ( l = 0; l < cdim; l++ )
        cp[k+l] = rcp[j+l]/w;
    }
    free ( obj->cpoints );
    if ( obj->savedcpoints )
      { free ( obj->savedcpoints );  obj->savedcpoints = NULL; }
    obj->savedsize = 0;
    obj->cpoints = cp;
    obj->me.cpdimen = cdim;
    obj->me.dlistmask = 0;
    obj->rational = false;
    GeomObjectProcessDependencies ( (geom_object*)obj );
  }
  return true;
} /*GeomObjectBSplineCurveSetNonRational*/

void GeomObjectDrawBSplineCurve ( GO_BSplineCurve *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  if (obj->me.dlistmask & BSC_DLM_CURVE )
    glCallList ( obj->me.displaylist+1 );
  else {
    glNewList ( obj->me.displaylist+1, GL_COMPILE_AND_EXECUTE );
    glColor3fv ( xglec_White );
    if ( obj->degree == 1 )
      DrawAPolyline ( obj->me.cpdimen, obj->me.spdimen, obj->lastknot-1, obj->cpoints );
    else
      DrawBSplineCurve ( obj->me.cpdimen, obj->me.spdimen, obj->degree, obj->lastknot,
                         obj->knots, obj->cpoints, 8*obj->degree );
    glEndList ();
    obj->me.dlistmask |= BSC_DLM_CURVE;
  }
} /*GeomObjectDrawBSplineCurve*/

void GeomObjectDrawBSplineCPoly ( GO_BSplineCurve *obj )
{
  int ncp;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  if ( obj->me.dlistmask & BSC_DLM_CPOLY )
    glCallList ( obj->me.displaylist );
  else {
    ncp = obj->lastknot-obj->degree;
    if ( obj->closed )
      ncp -= obj->degree-1;
    glNewList ( obj->me.displaylist, GL_COMPILE_AND_EXECUTE );
    glColor3fv ( xglec_Green );
    if ( obj->rational && obj->weightpoints ) {
      DrawLineSegments ( obj->me.cpdimen, obj->me.spdimen, ncp-1,
                         obj->cpoints, obj->weightpoints );
      glColor3fv ( xglec_Green4 );
      DrawLineSegments ( obj->me.cpdimen, obj->me.spdimen, ncp-1,
                         obj->weightpoints, &obj->cpoints[(int)obj->me.cpdimen] );
    }
    else
      DrawAPolyline ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints );
    glColor3fv ( xglec_Yellow );
    DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                  obj->mkcp, 0, MASK_CP_MOVEABLE );
    glColor3fv ( xglec_OrangeRed );
    DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                  obj->mkcp, marking_mask, MASK_MARKED | MASK_CP_MOVEABLE );
    if ( obj->rational ) {
      glColor3fv ( xglec_DarkOliveGreen2 );
      DrawPoints ( obj->me.cpdimen, obj->me.spdimen, ncp-1, obj->weightpoints );
    }
    glEndList ();
    obj->me.dlistmask |= BSC_DLM_CPOLY;
  }
} /*GeomObjectDrawBSplineCPoly*/

void GeomObjectDrawBSplineBPoly ( GO_BSplineCurve *obj )
{
  void   *sp;
  int    dim, deg, lkn, kpcs, olkn, i, k;
  double *bpoly;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  if ( obj->me.dlistmask & BSC_DLM_BPOLY )
    glCallList ( obj->me.displaylist+2 );
  else {
    sp = pkv_GetScratchMemTop ();
    dim = obj->me.cpdimen;
    deg = obj->degree;
    lkn = obj->lastknot;
    kpcs = mbs_NumKnotIntervalsd ( deg, obj->lastknot, obj->knots );
    bpoly = pkv_GetScratchMemd ( kpcs*dim*(deg+1) );
    if ( bpoly ) {
      mbs_multiBSCurvesToBezd ( dim, 1, deg, lkn, obj->knots, 0, obj->cpoints,
                                &kpcs, &olkn, NULL, 0, bpoly );
      glNewList ( obj->me.displaylist+2, GL_COMPILE_AND_EXECUTE );
      glColor3fv ( xglec_Green3 );
      for ( i = k = 0;  i < kpcs;  i++, k += dim*(deg+1) )
        DrawAPolyline ( dim, obj->me.spdimen, deg+1, &bpoly[k] );
      glEndList ();
      obj->me.dlistmask |= BSC_DLM_BPOLY;
    }
    pkv_SetScratchMemTop ( sp );
  }
} /*GeomObjectDrawBSplineBPoly*/

void GeomObjectDrawBSplineCurvature ( GO_BSplineCurve *obj )
{
  void     *sp;
  int      dim, deg, lkn, kpcs, olkn, i, j, k, dens;
  double   *bpoly;
  double   t;
  point2d  p2, q2;
  vector2d f2[2];
  point3d  p3, q3;
  vector3d f3[3];
  double   curv[2], scf;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  if ( obj->me.dlistmask & BSC_DLM_CURVATURE )
    glCallList ( obj->me.displaylist+3 );
  else {
    sp = pkv_GetScratchMemTop ();
    dim = obj->me.cpdimen;
    deg = obj->degree;
    lkn = obj->lastknot;
    kpcs = mbs_NumKnotIntervalsd ( deg, obj->lastknot, obj->knots );
    bpoly = pkv_GetScratchMemd ( kpcs*dim*(deg+1) );
    if ( bpoly ) {
      mbs_multiBSCurvesToBezd ( dim, 1, deg, lkn, obj->knots, 0, obj->cpoints,
                                &kpcs, &olkn, NULL, 0, bpoly );
      dens = obj->graph_dens;
      scf = xge_LogSlidebarValued ( 0.01, 100.0, obj->curvature_scale );
      glNewList ( obj->me.displaylist+3, GL_COMPILE_AND_EXECUTE );
      glColor3fv ( xglec_CornflowerBlue );
      glBegin ( GL_LINES );
      switch ( obj->me.spdimen ) {
    case 2:
        if ( dim == 2 ) {      /* non-rational */
          for ( i = k = 0;  i < kpcs;  i++, k += dim*(deg+1) ) {
            for ( j = 0; j <= dens; j++ ) {
              t = (double)j/(double)dens;
              mbs_BCFrenetC2d ( deg, (point2d*)&bpoly[k], t, &p2, f2, curv );
              AddVector2Md ( &p2, &f2[1], curv[0]*scf, &q2 );
              glVertex2dv ( &p2.x );
              glVertex2dv ( &q2.x );
            }
          }
        }
        else if ( dim == 3 ) { /* rational */
          for ( i = k = 0;  i < kpcs;  i++, k += dim*(deg+1) ) {
            for ( j = 0; j <= dens; j++ ) {
              t = (double)j/(double)dens;
              mbs_BCFrenetC2Rd ( deg, (point3d*)&bpoly[k], t, &p2, f2, curv );
              AddVector2Md ( &p2, &f2[1], curv[0]*scf, &q2 );
              glVertex2dv ( &p2.x );
              glVertex2dv ( &q2.x );
            }
          }
        }
        break;
    case 3:
        if ( dim == 3 ) {      /* non-rational */
          for ( i = k = 0;  i < kpcs;  i++, k += dim*(deg+1) ) {
            for ( j = 0; j <= dens; j++ ) {
              t = (double)j/(double)dens;
              mbs_BCFrenetC3d ( deg, (point3d*)&bpoly[k], t, &p3, f3, curv );
              AddVector3Md ( &p3, &f3[1], curv[0]*scf, &q3 );
              glVertex3dv ( &p3.x );
              glVertex3dv ( &q3.x );
            }
          }
        }
        else if ( dim == 3 ) { /* rational */
          for ( i = k = 0;  i < kpcs;  i++, k += dim*(deg+1) ) {
            for ( j = 0; j <= dens; j++ ) {
              t = (double)j/(double)dens;
              mbs_BCFrenetC3Rd ( deg, (point4d*)&bpoly[k], t, &p3, f3, curv );
              AddVector3Md ( &p3, &f3[1], curv[0]*scf, &q3 );
              glVertex3dv ( &p3.x );
              glVertex3dv ( &q3.x );
            }
          }
        }
        break;
    default:
        break;
      }
      glEnd ();
      glEndList ();
      obj->me.dlistmask |= BSC_DLM_CURVATURE;
    }
    pkv_SetScratchMemTop ( sp );
  }
} /*GeomObjectDrawBSplineCurvature*/

void GeomObjectDrawBSplineTorsion ( GO_BSplineCurve *obj )
{
  void     *sp;
  int      dim, deg, lkn, kpcs, olkn, i, j, k, dens;
  double   *bpoly;
  double   t;
  point3d  p3, q3;
  vector3d f3[3];
  double   curv[2], scf;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  if ( obj->me.dlistmask & BSC_DLM_TORSION )
    glCallList ( obj->me.displaylist+4 );
  else {
    sp = pkv_GetScratchMemTop ();
    dim = obj->me.cpdimen;
    deg = obj->degree;
    lkn = obj->lastknot;
    kpcs = mbs_NumKnotIntervalsd ( deg, obj->lastknot, obj->knots );
    bpoly = pkv_GetScratchMemd ( kpcs*dim*(deg+1) );
    if ( bpoly ) {
      mbs_multiBSCurvesToBezd ( dim, 1, deg, lkn, obj->knots, 0, obj->cpoints,
                                &kpcs, &olkn, NULL, 0, bpoly );
      dens = obj->graph_dens;
      scf = xge_LogSlidebarValued ( 0.01, 100.0, obj->torsion_scale );
      glNewList ( obj->me.displaylist+4, GL_COMPILE_AND_EXECUTE );
      glColor3fv ( xglec_SpringGreen2 );
      glBegin ( GL_LINES );
      switch ( obj->me.spdimen ) {
    case 3:
        if ( dim == 3 ) {      /* non-rational */
          for ( i = k = 0;  i < kpcs;  i++, k += dim*(deg+1) ) {
            for ( j = 0; j <= dens; j++ ) {
              t = (double)j/(double)dens;
              mbs_BCFrenetC3d ( deg, (point3d*)&bpoly[k], t, &p3, f3, curv );
              AddVector3Md ( &p3, &f3[2], curv[1]*scf, &q3 );
              glVertex3dv ( &p3.x );
              glVertex3dv ( &q3.x );
            }
          }
        }
        else if ( dim == 3 ) { /* rational */
          for ( i = k = 0;  i < kpcs;  i++, k += dim*(deg+1) ) {
            for ( j = 0; j <= dens; j++ ) {
              t = (double)j/(double)dens;
              mbs_BCFrenetC3Rd ( deg, (point4d*)&bpoly[k], t, &p3, f3, curv );
              AddVector3Md ( &p3, &f3[2], curv[1]*scf, &q3 );
              glVertex3dv ( &p3.x );
              glVertex3dv ( &q3.x );
            }
          }
        }
        break;
    default:
        break;
      }
      glEnd ();
      glEndList ();
      obj->me.dlistmask |= BSC_DLM_TORSION;
    }
    pkv_SetScratchMemTop ( sp );
  }
} /*GeomObjectDrawBSplineTorsion*/

void GeomObjectDisplayBSplineCurve ( GO_BSplineCurve *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  if ( obj->view_curve )
    GeomObjectDrawBSplineCurve ( obj );
  if ( obj->view_cpoly )
    GeomObjectDrawBSplineCPoly ( obj );
  if ( obj->view_bpoly )
    GeomObjectDrawBSplineBPoly ( obj );
  if ( obj->view_curvature )
    GeomObjectDrawBSplineCurvature ( obj );
  if ( obj->view_torsion )
    GeomObjectDrawBSplineTorsion ( obj );
} /*GeomObjectDisplayBSplineCurve*/

boolean GeomObjectBSplineCurveSetDegree ( GO_BSplineCurve *obj, int deg )
{
  int    k, nkn, ncp, nnlkn, nncp;
  double *cp, *kn, *acp, *wcp, *akn;
  byte   *mkcp;

  kn = cp = acp = akn = wcp = NULL;
  mkcp = NULL;
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  if ( deg < 1 || deg > MAX_DEGREE )
    return false;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  k = mbs_NumKnotIntervalsd ( obj->degree, obj->lastknot, obj->knots );
  if ( deg > obj->degree ) {
        /* upper bounds for the numbers of new knots and control points */
    nkn = obj->lastknot + (k+1)*(deg-obj->degree);
    ncp = nkn - deg;
    kn = malloc ( (nkn+1)*sizeof(double) );
    cp = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
    if ( !kn || !cp )
      goto failure;
    if ( !obj->closed ) {
      if ( !mbs_multiBSDegElevd ( 1, obj->me.cpdimen, obj->degree, obj->lastknot,
                            obj->knots, 0, obj->cpoints, deg-obj->degree,
                            &deg, &nnlkn, kn, 0, cp, true ) )
        goto failure;
    }
    else {
      if ( !mbs_multiBSDegElevClosedd ( 1, obj->me.cpdimen, obj->degree, obj->lastknot,
                                  obj->knots, 0, obj->cpoints, deg-obj->degree,
                                  &deg, &nnlkn, kn, 0, cp ) )
        goto failure;
    }
  }
  else if ( deg < obj->degree ) {
        /* upper bounds for the numbers of new knots and control points */
    nkn = obj->lastknot;
    ncp = nkn - deg;
    kn = malloc ( (nkn+1)*sizeof(double) );
    cp = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
    if ( !kn || !cp )
      goto failure;
    if ( !obj->closed )
      mbs_multiBSDegRedd ( 1, obj->me.cpdimen, obj->degree, obj->lastknot,
                           obj->knots, 0, obj->cpoints, obj->degree-deg,
                           &deg, &nnlkn, kn, 0, cp );
    else
      mbs_multiBSDegRedClosedd ( 1, obj->me.cpdimen, obj->degree, obj->lastknot,
                                 obj->knots, 0, obj->cpoints, obj->degree-deg,
                                 &deg, &nnlkn, kn, 0, cp );
  }
        /* copy the new representation to the final arrays and */
        /* assign them to the objects */
  nncp = nnlkn-deg;
  if ( obj->rational ) {
    wcp = malloc ( (nncp-1)*obj->me.cpdimen*sizeof(double) );
    if ( !wcp )
      goto failure;
  }
  else
    wcp = NULL;
  akn = malloc ( (nnlkn+1)*sizeof(double) );
  acp = malloc ( nncp*obj->me.cpdimen*sizeof(double) );
  mkcp = malloc ( nncp );
  if ( !akn || !acp || !mkcp || (obj->rational && !wcp) )
    goto failure;
  memcpy ( akn, kn, (nnlkn+1)*sizeof(double) );
  memcpy ( acp, cp, nncp*obj->me.cpdimen*sizeof(double) );
  memset ( mkcp, MASK_CP_MOVEABLE, nncp );
  if ( obj->weightpoints ) {
    free ( obj->weightpoints );
    obj->weightpoints = wcp;
    GeomObjectSetupWeightPoints ( obj->me.cpdimen, nncp, acp, wcp );
  }
  free ( obj->knots );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->knots = akn;
  obj->cpoints = acp;
  obj->mkcp = mkcp;
  obj->maxknots = nnlkn+1;
  obj->degree = deg;
  obj->lastknot = nnlkn;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  obj->me.dlistmask = 0;
  GeomObjectProcessDependencies ( (geom_object*)obj );
  GeomObjectBSplineCurveDisplayInfoText ( obj );
  return true;

failure:
  if ( mkcp ) free ( mkcp );
  if ( kn )  free ( kn );
  if ( cp )  free ( cp );
  if ( akn ) free ( akn );
  if ( acp ) free ( acp );
  if ( wcp ) free ( wcp );
  return false;
} /*GeomObjectBSplineCurveSetDegree*/

void GeomObjectBSplineCurveFindBBox ( GO_BSplineCurve *obj, Box3d *box )
{
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  GeomObjectFindBBox ( obj->me.cpdimen, obj->rational,
                       obj->lastknot-obj->degree, obj->cpoints, box );
} /*GeomObjectBSplineCurveFindBBox*/

boolean GeomObjectBSplineCurveFindCPoint ( GO_BSplineCurve *obj,
                                           CameraRecd *CPos, short x, short y,
                                           int *dist )
{
  boolean ok, ok1;
  int     ncp;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  if ( obj->closed )
    ncp = obj->lastknot - 2*obj->degree;
  else
    ncp = obj->lastknot - obj->degree;
  ok = GeomObjectFindNearestPoint ( obj->me.cpdimen, obj->me.spdimen,
                                    ncp, obj->cpoints,
                                    obj->mkcp, MASK_CP_MOVEABLE, CPos, x, y, dist );
  if ( obj->rational ) {
    if ( obj->closed )
      ncp ++;
    ok1 = GeomObjectFindNearestPoint ( obj->me.cpdimen, obj->me.spdimen,
                                       ncp-1, obj->weightpoints,
                                       NULL, 0, CPos, x, y, dist );
    if ( ok1 )
      current_point_ind += obj->lastknot - obj->degree;
  }
  else
    ok1 = false;
  return ok || ok1;
} /*GeomObjectBSplineCurveFindCPoint*/

void GeomObjectBSplineCurveSetCPoint ( GO_BSplineCurve *obj,
                                       CameraRecd *CPos, short x, short y )
{
  int c, k, ncp, clcK, dim;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  dim = obj->me.cpdimen;
  ncp = obj->lastknot - obj->degree;
  if ( current_point_ind >= 0 && current_point_ind < ncp ) {
    k = dim*current_point_ind;
    GeomObjectSetPoint ( CPos, x, y, obj->me.cpdimen, obj->me.spdimen,
                         &obj->cpoints[k] );
    if ( obj->closed ) {
      clcK = obj->lastknot-2*obj->degree;
      if ( current_point_ind < obj->degree )
        memcpy ( &obj->cpoints[k+clcK*dim], &obj->cpoints[k], dim*sizeof(double) );
      else if ( current_point_ind >= clcK )
        memcpy ( &obj->cpoints[k-clcK*dim], &obj->cpoints[k], dim*sizeof(double));
    }
  }
  else if (obj->rational ) {  /* a Farin point has been moved */
    c = current_point_ind-ncp;
    k = dim*c;
    GeomObjectSetWeightPoint ( CPos, x, y, dim, &obj->cpoints[k],
                               &obj->weightpoints[k] );
    if ( obj->closed ) {
      clcK = obj->lastknot-2*obj->degree;
      k += dim;
      if ( c < obj->degree-1 )
        memcpy ( &obj->cpoints[k+clcK*dim], &obj->cpoints[k], dim*sizeof(double) );
      else if ( c >= clcK-1 )
        memcpy ( &obj->cpoints[k-clcK*dim], &obj->cpoints[k], dim*sizeof(double) );
    }
  }
  if ( obj->rational )
    GeomObjectSetupWeightPoints ( dim, ncp, obj->cpoints, obj->weightpoints );
  obj->me.dlistmask = 0;
  GeomObjectProcessDependencies ( (geom_object*)obj );
} /*GeomObjectBSplineCurveSetCPoint*/

void GeomObjectBSplineCurveMarkCPoints ( GO_BSplineCurve *obj,
                                         CameraRecd *CPos, Box2s *box,
                                         char mask, int action )
{
  int deg, ncp;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  deg = obj->degree;
  ncp = obj->lastknot - deg;
  GeomObjectMarkPoints ( obj->me.cpdimen, obj->me.spdimen, ncp,
                         obj->mkcp, obj->cpoints, CPos, box, mask, action );
  if ( obj->closed )
    memcpy ( &obj->mkcp[ncp-deg], obj->mkcp, deg*sizeof(byte) );
  obj->me.dlistmask &= ~BSC_DLM_CPOLY;
} /*GeomObjectBSplineCurveMarkCPoints*/

void GeomObjectBSplineCurveMarkCPoint ( GO_BSplineCurve *obj,
                                        char mask, int action )
{
  int deg, ncp;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  deg = obj->degree;
  ncp = obj->lastknot - deg;
  GeomObjectMarkPoint ( ncp, obj->mkcp, mask, action );
  if ( obj->closed )
    memcpy ( &obj->mkcp[ncp-deg], obj->mkcp, deg*sizeof(byte) );
  obj->me.dlistmask &= ~BSC_DLM_CPOLY;
} /*GeomObjectBezierCurveMarkCPoint*/

boolean GeomObjectBSplineCurveSaveCPoints ( GO_BSplineCurve *obj )
{
  int ncp, size;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  if ( obj->savedcpoints ) free ( obj->savedcpoints );
  ncp = obj->lastknot - obj->degree;
  size = ncp*(obj->me.cpdimen)*sizeof(double);
  obj->savedcpoints = malloc ( size );
  if ( !obj->savedcpoints )
    return false;
  memcpy ( obj->savedcpoints, obj->cpoints, size );
  obj->savedsize = size;
  return true;
} /*GeomObjectBSplineCurveSaveCPoints*/

void GeomObjectBSplineCurveUndoLastTransformation ( GO_BSplineCurve *obj )
{
  int ncp, size;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  ncp = obj->lastknot - obj->degree;
  size = ncp*(obj->me.cpdimen)*sizeof(double);
  if ( obj->savedcpoints && obj->savedsize == size ) {
    memcpy ( obj->cpoints, obj->savedcpoints, size );
    if ( obj->rational )
      GeomObjectSetupWeightPoints ( obj->me.cpdimen, ncp,
                                    obj->cpoints, obj->weightpoints );
    obj->me.dlistmask = 0;
    GeomObjectProcessDependencies ( (geom_object*)obj );
  }
} /*GeomObjectBSplineCurveUndoLastTransformation*/

void GeomObjectBSplineCurveTransformCPoints ( GO_BSplineCurve *obj,
                                              trans3d *tr, char mask )
{
  int ncp;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  ncp = obj->lastknot - obj->degree;
  if ( obj->savedcpoints ) {
    GeomObjectTransformPoints ( obj->me.cpdimen, obj->me.spdimen, ncp,
                                obj->mkcp, mask,
                                obj->savedcpoints, obj->cpoints, tr );
    if ( obj->rational )
      GeomObjectSetupWeightPoints ( obj->me.cpdimen, ncp,
                                    obj->cpoints, obj->weightpoints );
    obj->me.dlistmask = 0;
    GeomObjectProcessDependencies ( (geom_object*)obj );
  }
} /*GeomObjectBSplineCurveTransformCPoints*/

boolean GeomObjectBSplineCurveGetPointCoord ( GO_BSplineCurve *obj, int p,
                                  int *spdimen, int *cpdimen, double **pc )
{
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  if ( p < 0 || p >= obj->lastknot-obj->degree )
    return false;
  *spdimen = obj->me.spdimen;
  *cpdimen = obj->me.cpdimen;
  *pc = &obj->cpoints[obj->me.cpdimen*p];
  return true;
} /*GeomObjectBSplineCurveGetPointCoord*/

void GeomObjectBSplineCurveMoveKnot ( GO_BSplineCurve *obj, int kni )
{
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  obj->me.dlistmask &= ~BSC_DLM_CURVE;
  GeomObjectProcessDependencies ( (geom_object*)obj );
} /*GeomObjectBSplineCurveMoveKnot*/

boolean GeomObjectBSplineCurveInsertKnot ( GO_BSplineCurve *obj,
                                           double knot )
{
  int    dim, ncp, i, deg, lkn;
  double *cp, *kn, *wcp;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;

  cp = kn = wcp = NULL;
  mkcp = NULL;
  dim = obj->me.cpdimen;
  deg = obj->degree;
  lkn = obj->lastknot;
  ncp = lkn-deg;
  if ( knot < obj->knots[deg] || knot > obj->knots[lkn-deg] )
    goto failure;
  kn = malloc ( (lkn+2)*sizeof(double) );
  cp = malloc ( (ncp+1)*dim*sizeof(double) );
  if ( obj->rational )
    wcp = malloc ( ncp*dim*sizeof(double) );
  mkcp = malloc ( ncp+1 );
  if ( !kn || !cp || !mkcp || (obj->rational && !wcp) )
    goto failure;
  memcpy ( kn, obj->knots, (lkn+1)*sizeof(double) );
  memcpy ( cp, obj->cpoints, ncp*dim*sizeof(double) );
  if ( !obj->closed )
    i = mbs_multiKnotInsd ( deg, &lkn, kn, 1, dim, 0, 0, cp, knot );
  else
    i = mbs_multiKnotInsClosedd ( deg, &lkn, kn, 1, dim, 0, 0, cp, knot );
  memcpy ( mkcp, obj->mkcp, i );
  mkcp[i] = MASK_CP_MOVEABLE;
  memcpy ( &mkcp[i+1], &obj->mkcp[i], ncp-i );
  if ( obj->closed )
    memcpy ( &mkcp[lkn-2*deg], &mkcp, deg*sizeof(byte) );
  free ( obj->knots );
  free ( obj->cpoints );
  free ( obj->mkcp );
  if ( obj->rational ) {
    free ( obj->weightpoints );
    obj->weightpoints = wcp;
    GeomObjectSetupWeightPoints ( dim, lkn-deg, cp, wcp );
  }
  obj->knots = kn;
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  obj->lastknot = lkn;
  obj->maxknots = lkn+1;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  obj->me.dlistmask = 0;
  GeomObjectProcessDependencies ( (geom_object*)obj );
  GeomObjectBSplineCurveDisplayInfoText ( obj );
  return true;

failure:
  if ( wcp )  free ( wcp );
  if ( cp )   free ( cp );
  if ( kn )   free ( kn );
  if ( mkcp ) free ( mkcp );
  return false;
} /*GeomObjectBSplineCurveInsertKnot*/

boolean GeomObjectBSplineCurveRemoveKnot ( GO_BSplineCurve *obj,
                                           int kni )
{
  int    dim, ncp, i, deg, lkn;
  double *cp, *kn;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;

  cp = kn = NULL;
  mkcp = NULL;
  dim = obj->me.cpdimen;
  deg = obj->degree;
  lkn = obj->lastknot;
  ncp = lkn-deg;
  if ( kni <= deg || kni >= lkn-deg )
    goto failure;
  kn = malloc ( (lkn+1)*sizeof(double) );
  cp = malloc ( ncp*dim*sizeof(double) );
  mkcp = malloc ( ncp-1 );
  if ( !kn || !cp || !mkcp )
    goto failure;
  memcpy ( kn, obj->knots, (lkn+1)*sizeof(double) );
  memcpy ( cp, obj->cpoints, ncp*dim*sizeof(double) );
  if ( !obj->closed )
    i = mbs_multiKnotRemoved ( deg, &lkn, kn, 1, dim, 0, 0, cp, kni );
  else
    i = mbs_multiKnotRemoveClosedd ( deg, &lkn, kn, 1, dim, 0, 0, cp, kni );
  if ( i > 1 )
    memcpy ( mkcp, obj->mkcp, i-1 );
  memcpy ( &mkcp[i-1], &obj->mkcp[i], ncp-i );
  if ( obj->closed )
    memcpy ( &mkcp[lkn-2*deg], &mkcp, deg*sizeof(byte) );
  free ( obj->knots );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->knots = kn;
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  obj->lastknot = lkn;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  if ( obj->rational )
    GeomObjectSetupWeightPoints ( dim, lkn-deg,
                                  obj->cpoints, obj->weightpoints );
  obj->me.dlistmask = 0;
  GeomObjectProcessDependencies ( (geom_object*)obj );
  GeomObjectBSplineCurveDisplayInfoText ( obj );
  return true;

failure:
  if ( cp ) free ( cp );
  if ( kn ) free ( kn );
  if ( mkcp ) free ( mkcp );
  return false;
} /*GeomObjectBSplineCurveRemoveKnot*/

boolean GeomObjectBSplineCurveSetUniformKnots ( GO_BSplineCurve *obj, boolean uniform )
{
  int    i, lkn;
  double *kn;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  obj->uniform = uniform;
  if ( uniform ) {
    lkn = obj->lastknot;
    kn = obj->knots;
    for ( i = 0; i <= lkn; i++ )
      kn[i] = (double)i;
    obj->me.dlistmask &= ~BSC_DLM_CURVE;
    GeomObjectProcessDependencies ( (geom_object*)obj );
  }
  return true;
} /*GeomObjectBSplineCurveSetUniformKnots*/

boolean GeomObjectBSplineCurveRefine ( GO_BSplineCurve *obj )
{
  void   *sp;
  int    i, lkn, deg;
  double *kn, *acp, *wcp;
  byte   *mkcp;

  sp = pkv_GetScratchMemTop ();
  kn = acp = wcp = NULL;
  mkcp = NULL;
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    goto failure;

  deg = obj->degree;
  lkn = obj->lastknot;
  kn = malloc ( (2*(lkn-deg)+1)*sizeof(double) );
  acp = malloc ( (2*lkn-3*deg)*obj->me.cpdimen*sizeof(double) );
  mkcp = malloc ( 2*lkn-3*deg );
  if ( !kn || !acp || !mkcp )
    goto failure;
  if ( obj->rational ) {
    wcp = malloc (  (2*lkn-3*deg-1)*obj->me.cpdimen*sizeof(double) );
    if ( !wcp )
      goto failure;
  }
  else
    wcp = NULL;
  if ( !mbs_multiLaneRiesenfeldd ( obj->me.cpdimen, 1, deg, lkn,
                                   0, obj->cpoints, &lkn, 0, acp ) )
    goto failure;

  free ( obj->knots );
  free ( obj->cpoints );
  if ( obj->weightpoints ) free ( obj->weightpoints );
  free ( obj->mkcp );
  obj->lastknot = lkn;
  obj->knots = kn;
  obj->cpoints = acp;
  obj->weightpoints = wcp;
  obj->mkcp = mkcp;
  for ( i = 0; i <= lkn; i++ )
    kn[i] = (double)i;
  for ( i = 0; i < lkn-deg; i++ )
    mkcp[i] = MASK_CP_MOVEABLE;
  if ( obj->rational )
    GeomObjectSetupWeightPoints ( obj->me.cpdimen, lkn-deg,
                                  obj->cpoints, obj->weightpoints );
  obj->me.dlistmask = 0;
  pkv_SetScratchMemTop ( sp );
  GeomObjectProcessDependencies ( (geom_object*)obj );
  return true;

failure:
  if ( kn )   free ( kn );
  if ( acp )  free ( acp );
  if ( wcp )  free ( wcp );
  if ( mkcp ) free ( mkcp );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBSplineCurveRefine*/

boolean GeomObjectBSplineCurveSetClosed ( GO_BSplineCurve *obj, boolean closed )
{
  int    deg, lkn, clcK, cpdim;
  double *kn, *cp, clcT;
  int    i;

  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  if ( !closed ) {
    obj->closed = false;
    return true;
  }
  deg = obj->degree;
  lkn = obj->lastknot;
  kn = obj->knots;
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
  cp = obj->cpoints;
  pkn_AddMatrixd ( 1, cpdim*deg, 0, cp, 0, &cp[cpdim*clcK], 0, cp );
  pkn_MultMatrixNumd ( 1, cpdim*deg, 0, cp, 0.5, 0, cp );
  memcpy ( &cp[cpdim*clcK], cp, cpdim*deg*sizeof(double) );
  memcpy ( &obj->mkcp[clcK], obj->mkcp, deg*sizeof(byte) );
  if ( obj->rational )
    GeomObjectSetupWeightPoints ( cpdim, lkn-deg,
                                  obj->cpoints, obj->weightpoints );
  obj->me.dlistmask = 0;
  obj->closed = true;
  GeomObjectProcessDependencies ( (geom_object*)obj );
  GeomObjectBSplineCurveDisplayInfoText ( obj );
  return true;
} /*GeomObjectBSplineCurveSetClosed*/

boolean GeomObjectBSplineCurveSetMengerc ( GO_BSplineCurve *obj, boolean mengerc )
{
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  if ( mengerc &&
       (obj->degree < 3 || !obj->closed ||
        obj->rational || obj->me.spdimen != 3) ) {
    obj->mengerc = false;
    return false;
  }
  if ( mengerc )
    GeomObjectBSplineCurveSetUniformKnots ( obj, true );
  obj->mengerc = mengerc;
  return true;
} /*GeomObjectBSplineCurveSetMengerc*/

void GeomObjectBSplineCurveSetCurvatureGraph ( GO_BSplineCurve *obj,
                                boolean view_curvature, double curvature_scale,
                                boolean view_torsion, double torsion_scale,
                                int graph_dens )
{
  if ( (view_curvature && !obj->view_curvature) ||
       curvature_scale != obj->curvature_scale ||
       graph_dens != obj->graph_dens )
    obj->me.dlistmask &= ~BSC_DLM_CURVATURE;
  obj->view_curvature = view_curvature;
  if ( (view_torsion && !obj->view_torsion) ||
       torsion_scale != obj->torsion_scale ||
       graph_dens != obj->graph_dens )
    obj->me.dlistmask &= ~BSC_DLM_TORSION;
  obj->curvature_scale = curvature_scale;
  obj->view_torsion = view_torsion && (obj->me.spdimen == 3);
  obj->torsion_scale = torsion_scale;
  obj->graph_dens = graph_dens;
} /*GeomObjectBSplineCurveSetCurvatureGraph*/

boolean GeomObjectWriteBSCAttributes ( GO_BSplineCurve *obj )
{
  return true;
} /*GeomObjectWriteBSCAttributes*/

boolean GeomObjectWriteBSplineCurve ( GO_BSplineCurve *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return false;
  return bsf_WriteBSplineCurved ( obj->me.spdimen, obj->me.cpdimen,
                                  obj->rational, obj->degree,
                                  obj->lastknot, obj->knots, obj->closed,
                                  obj->cpoints, obj->me.name, obj->me.ident,
                                  GeomObjectWriteAttributes,
                                  (void*)&obj->me );
} /*GeomObjectWriteBSplineCurve*/

boolean GeomObjectBSCResolveDependencies ( GO_BSplineCurve *obj )
{
  return true;
} /*GeomObjectBSCResolveDependencies*/

void GeomObjectReadBSplineCurve ( void *usrdata, const char *name, int ident,
                                  int degree, int lastknot,
                                  const double *knots, boolean closed,
                                  const point4d *cpoints, int spdimen,
                                  boolean rational )
{
  GO_BSplineCurve *obj;
  double          *cp, *wcp, *kn;
  byte            *mkcp;
  int             ncp, cpdimen;

  obj = (GO_BSplineCurve*)GeomObjectAddBSplineCurve ( name, spdimen, rational );
  if ( obj ) {
    obj->me.ident = ident;
    ncp = lastknot - degree;
    cpdimen = obj->me.cpdimen;
    cp = malloc ( ncp*cpdimen*sizeof(double) );
    kn = malloc ( (lastknot+1)*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( rational )
      wcp = malloc ( (ncp-1)*cpdimen*sizeof(double) );
    else
      wcp = NULL;
    if ( !cp || !kn || !mkcp || (rational && !wcp) ) {
      if ( wcp )  free ( wcp );
      if ( mkcp ) free ( mkcp );
      if ( kn )   free ( kn );
      if ( cp )   free ( cp );
      GeomObjectDeleteBSplineCurve ( obj );
      return;
    }
    if ( obj->weightpoints )
      free ( obj->weightpoints );
    free ( obj->cpoints );
    free ( obj->knots );
    free ( obj->mkcp );
    if ( obj->savedcpoints ) {
      free ( obj->savedcpoints );
      obj->savedcpoints = NULL;
    }
    obj->cpoints = cp;
    obj->weightpoints = wcp;
    obj->knots = kn;
    obj->mkcp = mkcp;
    obj->degree = degree;
    obj->lastknot = lastknot;
    obj->closed = closed;
    obj->maxknots = lastknot+1;
    GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                               ncp, (double*)cpoints, cp );
    if ( closed )
      memcpy ( &mkcp[lastknot-2*degree], mkcp, degree*sizeof(byte) );
    if ( rational )
      GeomObjectSetupWeightPoints ( cpdimen, ncp, cp, wcp );
    memcpy ( kn, knots, (lastknot+1)*sizeof(double) );
  }
} /*GeomObjectReadBSplineCurve*/

void GeomObjectBSplineCurveOutputToRenderer ( GO_BSplineCurve *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_CURVE )
    return;
  if ( obj->me.spdimen != 3 )
    return;
  if ( obj->rational )
    RendEnterBSCurve3Rd ( obj->degree, obj->lastknot, obj->knots,
                          (point4d*)obj->cpoints, 0.5*obj->pipe_diameter,
                          obj->me.colour );
  else
    RendEnterBSCurve3d ( obj->degree, obj->lastknot, obj->knots,
                         (point3d*)obj->cpoints, 0.5*obj->pipe_diameter,
                         obj->me.colour );
} /*GeomObjectBSplineCurveOutputToRenderer*/

void GeomObjectBSplineCurveDisplayInfoText ( GO_BSplineCurve *obj )
{
  char s[160];

  sprintf ( s, " deg = %d, lkn = %d", obj->degree, obj->lastknot );
  SetStatusText ( s, true );
} /*GeomObjectBSplineCurveDisplayInfoText*/

boolean GeomObjectBSplineCurveProcessDep ( GO_BSplineCurve *obj, geom_object *go )
{
  return false;
} /*GeomObjectBSplineCurveProcessDep*/

void GeomObjectBSplineCurveProcessDeletedDep ( GO_BSplineCurve *obj, geom_object *go )
{
} /*GeomObjectBSplineCurveProcessDeletedDep*/

