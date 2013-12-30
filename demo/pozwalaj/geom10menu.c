
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
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

#include "widgets.h"
#include "editor.h"  
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"
#include "pozwalaj.h"
#include "render.h"


xge_widget *InitGeom10Menu ( xge_widget *prev )
{
  xge_widget *menu;

  geom10win2D = xge_New2Dwind ( win1, NULL, GEOMWIN1_2D,
                                xge_WIDTH-SIDEMENUWIDTH,
                                xge_HEIGHT-TOPMENUHEIGHT,
                                SIDEMENUWIDTH, TOPMENUHEIGHT,
                                &g10win2D, DrawG10win2D );
        /* knots for B-spline curves */
  geom10knotwin = xge_NewKnotWind ( win1, NULL, GEOMWIN1_KN,
                                    xge_WIDTH-SIDEMENUWIDTH, 5*MAX_DEGREE+6,
                                    SIDEMENUWIDTH, TOPMENUHEIGHT,
                                    &g10knotwin, 0, NULL );
        /* knots and domain for B-spline patches */
  geom10t2knotwin = xge_NewT2KnotWind ( win1, NULL, GEOMWIN1_T2KN,
                                        xge_WIDTH-SIDEMENUWIDTH,
                                        xge_HEIGHT-TOPMENUHEIGHT,
                                        SIDEMENUWIDTH, TOPMENUHEIGHT,
                                        5*MAX_DEGREE,
                                        &g10t2knotwin, DrawG10winT2KN,
                                        0, NULL, 0, NULL );
        /* knots and windows for B-spline curves, which define the */
        /* spherical product */
  geom10win2Deqmer = xge_New2Dwind ( win1, NULL, GEOMWIN1_2DEQMER,
                                     xge_WIDTH-SIDEMENUWIDTH,
                                     xge_HEIGHT-TOPMENUHEIGHT-5*MAX_DEGREE-10,
                                     SIDEMENUWIDTH, TOPMENUHEIGHT,
                                     &g10win2Deqmer, DrawG10win2Deqmer );
  geom10knotwineqmer = xge_NewKnotWind ( win1, geom10win2Deqmer, GEOMWIN1_KNEQMER,
                                     xge_WIDTH-SIDEMENUWIDTH, 5*MAX_DEGREE+6,
                                     SIDEMENUWIDTH, TOPMENUHEIGHT+geom10win2Deqmer->h+4,
                                     &g10knotwineqmer, 0, NULL );
  geom10win = NULL;
  menu = xge_NewFMenu ( win1, prev, GEOMMENU1, xge_WIDTH-SIDEMENUWIDTH,
                       xge_HEIGHT-TOPMENUHEIGHT, SIDEMENUWIDTH, TOPMENUHEIGHT,
                       geom10win );
  return menu;
} /*InitGeom10Menu*/

void Geom10MenuResize ( void )
{
  int w, h;

  w = geom10menu->w;
  h = geom10menu->h;
  if ( geom10win == geom10knotwineqmer ) {
    geom10win2Deqmer->msgproc ( geom10win2Deqmer, xgemsg_RESIZE, 0,
                                w, h-5*MAX_DEGREE-10 );
    geom10knotwineqmer->y = TOPMENUHEIGHT+h-5*MAX_DEGREE-6;
    geom10knotwineqmer->msgproc ( geom10knotwineqmer, xgemsg_RESIZE, 0,
                                  w, 5*MAX_DEGREE+6 );
  }
  else if ( geom10win == geom10knotwin )
    geom10win->msgproc ( geom10win, xgemsg_RESIZE, 0, w, 5*MAX_DEGREE+6 );
  else if ( geom10win )
    geom10win->msgproc ( geom10win, xgemsg_RESIZE, 0, w, h );
} /*Geom10MenuResize*/

int Geom10MenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( er->id ) {
case GEOMWIN1_2D:
    return Geom10win2DCallBack ( er, msg, key, x, y );
case GEOMWIN1_KN:
    return Geom10winKNCallBack ( er, msg, key, x, y );
case GEOMWIN1_T2KN:
    return Geom10winT2KNCallBack ( er, msg, key, x, y );
case GEOMWIN1_2DEQMER:
    return Geom10win2DeqmerCallBack ( er, msg, key, x, y );
case GEOMWIN1_KNEQMER:
    return Geom10winKNeqmerCallBack ( er, msg, key, x, y );
default:
    return 0;
  }
} /*Geom10MenuCallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawG10win2D ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *_2Dwin;

  _2Dwin = er->data0;
  glXWaitX ();
  xgle_DrawGeomWinBackground ( er, 0 );
  if ( _2Dwin->display_coord && _2Dwin->inside )
    xgle_2DwindDrawCursorPos ( _2Dwin, xge_xx, xge_yy );
  xgle_SetGLCamerad ( &_2Dwin->CPos );

  glXWaitGL ();
  xgleCopyGLRect ( er->w, er->h, er->x, er->y );
  xge_DrawGeomWinSelectionRect ( er, &_2Dwin->selection_rect );
  xge_2DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*DrawG10win2D*/

int Geom10win2DCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
default:
    return 0;
  }
} /*Geom10win2DCallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void Geom10winKNSetKnots ( xge_KnotWind *kwin, int degree, int lastknot,
                           double *knots, boolean closed )
{
  kwin->degree = degree;
  kwin->lastknot = lastknot;
  kwin->maxknots = lastknot+2;
  kwin->knots = knots;
  kwin->closed = closed;
  if ( closed ) {
    kwin->clcK = lastknot-2*degree;
    kwin->clcT = knots[lastknot-degree]-knots[degree];
  }
} /*Geom10winKNSetKnots*/

static int _Geom10winKNCallBack ( xge_widget *er, int msg, int key,
                                  short x, short y, GO_BSplineCurve *obj )
{
  xge_KnotWind *knw;

  if ( !obj )
    return 0;
  knw = er->data0;
  switch ( msg ) {
case xgemsg_KNOTWIN_CHANGE_KNOT:
    switch ( er->id ) {
  case GEOMWIN1_KN:
  case GEOMWIN1_KNEQMER:
      GeomObjectBSplineCurveMoveKnot ( obj, knw->current_knot );
      if ( RenderingIsOn )
        StopRendering ();
      rendered_picture = false;
      xge_RedrawAll ();
      return 1;
  default:
      return 0;
    }

case xgemsg_KNOTWIN_INSERT_KNOT:
    switch ( er->id ) {
  case GEOMWIN1_KN:
  case GEOMWIN1_KNEQMER:
      if ( GeomObjectBSplineCurveInsertKnot ( obj, knw->newknot ) ) {
        Geom10winKNSetKnots ( knw, obj->degree, obj->lastknot,
                              obj->knots, obj->closed );
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_KNOTWIN_REMOVE_KNOT:
    switch ( er->id ) {
  case GEOMWIN1_KN:
  case GEOMWIN1_KNEQMER:
      if ( GeomObjectBSplineCurveRemoveKnot ( (GO_BSplineCurve*)current_go,
                                              knw->current_knot ) ) {
        Geom10winKNSetKnots ( knw, obj->degree, obj->lastknot,
                              obj->knots, obj->closed );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_KNOTWIN_CHANGE_MAPPING:
    switch ( er->id ) {
  case GEOMWIN1_KN:
  case GEOMWIN1_KNEQMER:
      return 1;
  default:
      return 0;
    }

case xgemsg_KNOTWIN_ERROR:
printf ( "error\n" );
    return 0;

default:
    return 0;
  }
} /*_Geom10winKNCallBack*/

int Geom10winKNCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  return _Geom10winKNCallBack ( er, msg, key, x, y,
                                (GO_BSplineCurve*)current_go );
} /*Geom10winKNCallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawG10winT2KNKnots ( xge_T2KnotWind *kwin )
{
  xge_widget *er;
  int        i;
  int        x, y, uxa, uxb, vya, vyb;
  double     uu;
  int        degu, lknu, degv, lknv;
  double     *knu, *knv;

  degu = kwin->degree_u;
  lknu = kwin->lastknot_u;
  knu  = kwin->knots_u;
  degv = kwin->degree_v;
  lknv = kwin->lastknot_v;
  knv  = kwin->knots_v;

  er = kwin->er;
  xgeSetForeground ( xgec_Green4 );
  x = kwin->CPos.xmin;
  y = kwin->CPos.ymin + kwin->CPos.height;
  xgeDrawLine ( kwin->CPos.xmin-5*MAX_DEGREE+5, y,
                kwin->CPos.xmin+kwin->CPos.width-5, y );
  xgeDrawLine ( x, kwin->CPos.ymin+5,
                x, kwin->CPos.ymin+kwin->CPos.height+5*MAX_DEGREE-5 );
  xgeSetForeground ( xgec_White );

  uxa = xge_T2KnotWindMapKnotU ( kwin, knu[degu] );
  uxa = max ( er->x+1, uxa );
  uxb = xge_T2KnotWindMapKnotU ( kwin, knu[lknu-degu] );
  uxb = min ( er->x+er->w-2, uxb );
  xgeDrawLine ( uxa, y, uxb, y );
  for ( uu = (double)(knu[1]-1.0), uxb = 0, i = 1;  i < lknu;  i++ ) {
    if ( knu[i] > uu ) {
      uxb = y + 5;
      uxa = xge_T2KnotWindMapKnotU ( kwin, (uu = knu[i]) );
    }
    else
      uxb += 5;
    xgeFillRectangle ( 3, 3, uxa-1, uxb-1 );
  }
  
  xgeSetForeground ( xgec_White );
  vya = xge_T2KnotWindMapKnotV ( kwin, knv[degv] );
  vya = max ( er->y+1, vya );
  vyb = xge_T2KnotWindMapKnotV ( kwin, knv[lknv-degv] );
  vyb = min ( er->y+er->h-2, vyb );
  xgeDrawLine ( x, vya, x, vyb );
  for ( uu = (double)(knv[1]-1.0), vyb = 0, i = 1;  i < lknv;  i++ ) {
    if ( knv[i] > uu ) {
      vyb = x - 5;
      vya = xge_T2KnotWindMapKnotV ( kwin, (uu = knv[i]) );
    }
    else
      vyb -= 5;
    xgeFillRectangle ( 3, 3, vyb-1, vya-1 );
  }
} /*DrawG10winT2KNKnots*/

void DrawG10winT2KNKnotLines ( xge_T2KnotWind *kwin )
{
  int    i;
  int    uxa, uxb, vya, vyb;
  double uu;
  int    degu, lknu, degv, lknv;
  double *knu, *knv;

  degu = kwin->degree_u;
  lknu = kwin->lastknot_u;
  knu  = kwin->knots_u;
  degv = kwin->degree_v;
  lknv = kwin->lastknot_v;
  knv  = kwin->knots_v;

  xgeSetForeground ( xgec_Grey2 );

  vya = xge_T2KnotWindMapKnotV ( kwin, knv[degv] );
  vyb = xge_T2KnotWindMapKnotV ( kwin, knv[lknv-degv] );
  for ( uu = (double)(knu[1]-1.0), i = degu;  i <= lknu-degu;  i++ ) {
    if ( knu[i] > uu ) {
      uxa = xge_T2KnotWindMapKnotU ( kwin, (uu = knu[i]) );
      xgeDrawLine ( uxa, vya, uxa, vyb );
    }
  }
  
  uxa = xge_T2KnotWindMapKnotU ( kwin, knu[degu] );
  uxb = xge_T2KnotWindMapKnotU ( kwin, knu[lknu-degu] );
  for ( uu = (double)(knv[1]-1.0), i = degv;  i <= lknv-degv;  i++ ) {
    if ( knv[i] > uu ) {
      vya = xge_T2KnotWindMapKnotV ( kwin, (uu = knv[i]) );
      xgeDrawLine ( uxa, vya, uxb, vya );
    }
  }
} /*DrawG10winT2KNKnotLines*/

void DrawG10winT2KNDomain ( xge_T2KnotWind *kwin )
{
  DrawG10winT2KNKnots ( kwin );
  DrawG10winT2KNKnotLines ( kwin );
} /*DrawG10winT2KNDomain*/

void DrawG10winT2KNDomainNet ( xge_T2KnotWind *kwin )
{
  void   *sp;
  int    i, j/*, pitch*/;
  int    xy, x0, x1, y0, y1;
  int    *xg, *yg;
  int    degu, lknu, degv, lknv;
  double *knu, *knv;

  degu = kwin->degree_u;
  lknu = kwin->lastknot_u;
  knu  = kwin->knots_u;
  degv = kwin->degree_v;
  lknv = kwin->lastknot_v;
  knv  = kwin->knots_v;

  x0 = kwin->CPos.xmin-5*MAX_DEGREE;
  xy = xge_T2KnotWindMapKnotU ( kwin, mbs_GrevilleAbscissad ( degu, knu, 0 ) );
  x0 = max ( x0, xy );
  x1 = kwin->CPos.xmin+kwin->CPos.width;
  xy = xge_T2KnotWindMapKnotU ( kwin, mbs_GrevilleAbscissad ( degu, knu,
                                lknu-degu-1 ) );
  x1 = min ( x1, xy );
  if ( x1 < x0 )
    return;
  y0 = kwin->CPos.ymin;
  xy = xge_T2KnotWindMapKnotV ( kwin, mbs_GrevilleAbscissad ( degv, knv,
                                lknv-degv-1 ) );
  y0 = max ( y0, xy );
  y1 = kwin->CPos.ymin+kwin->CPos.height+5*MAX_DEGREE;
  xy = xge_T2KnotWindMapKnotV ( kwin, mbs_GrevilleAbscissad ( degv, knv, 0 ) );
  y1 = min ( y1, xy );
  if ( y1 < y0 )
    return;

  sp = pkv_GetScratchMemTop ();
  xg = pkv_GetScratchMem ( (lknu-degu)*sizeof(int) );
  yg = pkv_GetScratchMem ( (lknv-degv)*sizeof(int) );
  if ( !xg || !yg )
    goto leave_it;

  xgeSetForeground ( xgec_Green1 );
  for ( i = 0; i < lknu-degu; i++ ) {
    xg[i] = xy = xge_T2KnotWindMapKnotU ( kwin, mbs_GrevilleAbscissad ( degu, knu, i ) );
    if ( xy >= x0 && xy <= x1 )
      xgeDrawLine ( xy, y0, xy, y1 );
  }
  for ( i = 0; i < lknv-degv; i++ ) {
    yg[i] = xy = xge_T2KnotWindMapKnotV ( kwin, mbs_GrevilleAbscissad ( degv, knv, i ) );
    if ( xy >= y0 && xy <= y1 )
      xgeDrawLine ( x0, xy, x1, xy );
  }

  /*pitch = lknv-degv;*/
  xgeSetForeground ( xgec_Yellow );
  for ( i = 0; i < lknu-degu; i++ )
    if ( xg[i] >= x0 && xg[i] <= x1 )
      for ( j = 0;  j < lknv-degv;  j++ ) {
        if ( yg[j] >= y0 && yg[j] <= y1 /*&&
             (mkcp[i*pitch+j] == MASK_CP_MOVEABLE)*/ )
          xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
      }
/*
  xgeSetForeground ( xgec_OrangeRed );
  for ( i = 0; i < lknu-degu; i++ )
    if ( xg[i] >= x0 && xg[i] <= x1 )
      for ( j = 0; j < lknv-degv; j++ ) {
        k = i*pitch+j;
        if ( yg[j] >= y0 && yg[j] <= y1 &&
             (mkcp[k] == (MASK_CP_MOVEABLE | MASK_CP_MARKED)) )
          xgeFillRectangle ( 3, 3, xg[i]-1, yg[j]-1 );
      }
*/
leave_it:
  pkv_SetScratchMemTop ( sp );
} /*DrawG10winT2KNDomainNet*/

void DrawG10winT2KN ( xge_widget *er, boolean onscreen )
{
  xge_T2KnotWind *kwin;

  kwin = er->data0;
  xge_DrawGeomWinBackground ( er );
  if ( kwin->display_coord && kwin->inside )
    xge_T2KnotWindDrawCursorPos ( kwin, xge_xx, xge_yy );
  DrawG10winT2KNDomain ( kwin );
  /* if ( ?? ) */
    DrawG10winT2KNDomainNet ( kwin );
  xge_DrawGeomWinSelectionRect ( er, &kwin->selection_rect );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*DrawG10winT2KN*/

void Geom10winT2SetKnots ( xge_T2KnotWind *kwin,
                           int degreeu, int lastknotu, double *knotsu, boolean closedu,
                           int degreev, int lastknotv, double *knotsv, boolean closedv )
{
  kwin->degree_u = degreeu;
  kwin->lastknot_u = lastknotu;
  kwin->maxknots_u = lastknotu+2;
  kwin->knots_u = knotsu;
  kwin->closed_u = closedu;
  if ( closedu ) {
    kwin->clcKu = lastknotu-2*degreeu;
    kwin->clcTu = knotsu[lastknotu-degreeu]-knotsu[degreeu];
  }
  kwin->degree_v = degreev;
  kwin->lastknot_v = lastknotv;
  kwin->maxknots_v = lastknotv+2;
  kwin->knots_v = knotsv;
  kwin->closed_v = closedv;
  if ( closedv ) {
    kwin->clcKv = lastknotv-2*degreev;
    kwin->clcTv = knotsv[lastknotv-degreev]-knotsv[degreev];
  }
} /*Geom10winT2SetKnots*/

int Geom10winT2KNCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  GO_BSplinePatch *obj;

  obj = (GO_BSplinePatch*)current_go;
  switch ( msg ) {
case xgemsg_T2KNOTWIN_RESIZE:
    return 0;

case xgemsg_T2KNOTWIN_PROJCHANGE:
    return 0;

case xgemsg_T2KNOTWIN_CHANGE_KNOT_U:
    switch ( er->id ) {
  case GEOMWIN1_T2KN:
      GeomObjectBSplinePatchMoveKnotU ( obj, g10t2knotwin.current_item );
      if ( RenderingIsOn )
        StopRendering ();
      rendered_picture = false;
      xge_RedrawAll ();
      return 1;
  default:
      return 0;
    }

case xgemsg_T2KNOTWIN_INSERT_KNOT_U:
    switch ( er->id ) {
  case GEOMWIN1_T2KN:
      if ( GeomObjectBSplinePatchInsertKnotU ( obj, g10t2knotwin.newknot ) ) {
        Geom10winT2SetKnots ( &g10t2knotwin,
                              obj->degree_u, obj->lastknot_u, obj->knots_u,
                              obj->closed_u,
                              obj->degree_v, obj->lastknot_v, obj->knots_v,
                              obj->closed_v );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_T2KNOTWIN_REMOVE_KNOT_U:
    switch ( er->id ) {
  case GEOMWIN1_T2KN:
      if ( GeomObjectBSplinePatchRemoveKnotU ( obj, g10t2knotwin.current_item ) ) {
        Geom10winT2SetKnots ( &g10t2knotwin,
                              obj->degree_u, obj->lastknot_u, obj->knots_u,
                              obj->closed_u,
                              obj->degree_v, obj->lastknot_v, obj->knots_v,
                              obj->closed_v );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_T2KNOTWIN_CHANGE_KNOT_V:
    switch ( er->id ) {
  case GEOMWIN1_T2KN:
      GeomObjectBSplinePatchMoveKnotV ( obj, g10t2knotwin.current_item );
      if ( RenderingIsOn )
        StopRendering ();
      rendered_picture = false;
      xge_RedrawAll ();
      return 1;
  default:
      return 0;
    }

case xgemsg_T2KNOTWIN_INSERT_KNOT_V:
    switch ( er->id ) {
  case GEOMWIN1_T2KN:
      if ( GeomObjectBSplinePatchInsertKnotV ( obj, g10t2knotwin.newknot ) ) {
        Geom10winT2SetKnots ( &g10t2knotwin,
                              obj->degree_u, obj->lastknot_u, obj->knots_u,
                              obj->closed_u,
                              obj->degree_v, obj->lastknot_v, obj->knots_v,
                              obj->closed_v );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_T2KNOTWIN_REMOVE_KNOT_V:
    switch ( er->id ) {
  case GEOMWIN1_T2KN:
      if ( GeomObjectBSplinePatchRemoveKnotV ( obj, g10t2knotwin.current_item ) ) {
        Geom10winT2SetKnots ( &g10t2knotwin,
                              obj->degree_u, obj->lastknot_u, obj->knots_u,
                              obj->closed_u,
                              obj->degree_v, obj->lastknot_v, obj->knots_v,
                              obj->closed_v );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_T2KNOTWIN_SELECT_POINTS:
    return 0;

case xgemsg_T2KNOTWIN_UNSELECT_POINTS:
    return 0;

case xgemsg_T2KNOTWIN_CHANGE_MAPPING:
    return 0;

case xgemsg_T2KNOTWIN_ERROR:
    return 0;

default:
    return 0;
  }
} /*Geom10winT2KNCallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawG10win2Deqmer ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *_2Dwin;

  _2Dwin = er->data0;
  glXWaitX ();
  xgle_DrawGeomWinBackground ( er, 0 );
  if ( _2Dwin->display_coord && _2Dwin->inside )
    xgle_2DwindDrawCursorPos ( _2Dwin, xge_xx, xge_yy );
  xgle_2DwindDrawAxes ( _2Dwin );
  xgle_SetGLCamerad ( &_2Dwin->CPos );
  if ( bsp_sproduct_eqmer ) {
    glPushMatrix ();
    GeomObjectDisplayBSplineCurve ( bsp_sproduct_eqmer );
    glPopMatrix ();
  }
  glXWaitGL ();
  xgleCopyGLRect ( er->w, er->h, er->x, er->y );
  xge_DrawGeomWinSelectionRect ( er, &_2Dwin->selection_rect );
  xge_2DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*DrawG10win2Deqmer*/

int Geom10win2DeqmerCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_2Dwind *ww;
  CameraRecd *CPos;
  boolean    found;
  int        spdimen, cpdimen;
  double     *pc;
  Box3d      box;

  ww = er->data0;
  CPos = &ww->CPos;
  switch ( msg ) {
case xgemsg_RESIZE:
    return 0;

case xgemsg_2DWIN_RESIZE:
case xgemsg_2DWIN_PROJCHANGE:
    return 0;

case xgemsg_2DWIN_PICK_POINT:
    found = GeomObjectFindObjCPoint ( 2, CPos, x, y,
                                      (geom_object*)bsp_sproduct_eqmer );
    if ( found ) {
      if ( GeomObjectGetPointCoord ( currentp_go, current_point_ind,
                                     &spdimen, &cpdimen, &pc ) )
         BottomDisplayPoint ( win0, spdimen, cpdimen, current_point_ind, pc, true );
    }
    return found;

case xgemsg_2DWIN_MOVE_POINT:
    GeomObjectSetCPoint ( CPos, x, y );
    if ( GeomObjectGetPointCoord ( currentp_go, current_point_ind,
                                   &spdimen, &cpdimen, &pc ) )
       BottomDisplayPoint ( win0, spdimen, cpdimen, current_point_ind, pc, true );
    if ( RenderingIsOn )
      StopRendering ();
    rendered_picture = false;
    xge_RedrawAll ();
    return 1;

case xgemsg_2DWIN_SELECT_POINTS:
    GeomObjectMarkCPoints ( 2, CPos, &ww->selection_rect, false );
    xge_SetClipping ( er );
    er->redraw ( er, true );
    return 1;

case xgemsg_2DWIN_UNSELECT_POINTS:
    GeomObjectMarkCPoints ( 2, CPos, &ww->selection_rect, true );
    xge_SetClipping ( er );
    er->redraw ( er, true );
    return 1;

case xgemsg_2DWIN_CHANGE_TRANS:
    Geom00Win2DShowTransOrigin ( ww );
    return 1;

case xgemsg_2DWIN_SAVE_POINTS:
    GeomObjectSaveCPoints ( 2 );
    return 1;

case xgemsg_2DWIN_TRANSFORM_POINTS:
    Geom00Win2DShowTransformation ( ww );
    GeomObjectTransformCPoints2D ( &ww->gwtrans, marking_mask );
    return 1;

case xgemsg_2DWIN_FIND_REFBBOX:
    GeomObjectBSplineCurveFindBBox ( bsp_sproduct_eqmer, &box );
    ww->RefBBox.x0 = box.x0;  ww->RefBBox.x1 = box.x1;
    ww->RefBBox.y0 = box.y0;  ww->RefBBox.y1 = box.y1;
    return 1;

case xgemsg_2DWIN_UNDO:
    GeomObjectUndoLastTransformation ( 3 );
    if ( RenderingIsOn )
      StopRendering ();
    rendered_picture = false;
    xge_RedrawAll ();
    return 1;

case xgemsg_2DWIN_ERROR:
    return 0;

default:
    return 0;
  }
} /*Geom10win2DeqmerCallBack*/

int Geom10winKNeqmerCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  return _Geom10winKNCallBack ( er, msg, key, x, y, bsp_sproduct_eqmer );
} /*Geom10winKNeqmerCallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void SetGeomWin10Empty ( void )
{
  if ( geom10win ) {
    geom10win = NULL;
    xge_SetMenuWidgets ( geom10menu, NULL, false );
    geom10menu->redraw = xge_DrawFMenu;
    bsp_sproduct_eqmer = NULL;
  }
} /*SetGeomWin10Empty*/

void SetGeomWin10Knotw ( GO_BSplineCurve *go )
{
  if ( geom10win != geom10knotwin ) {
    geom10win = geom10knotwin;
    xge_SetMenuWidgets ( geom10menu, geom10knotwin, false );
    Geom10MenuResize ();
    geom10menu->redraw = xge_DrawMenu;
    bsp_sproduct_eqmer = NULL;
  }
  xge_KnotWindFindMapping ( &g10knotwin );
} /*SetGeomWin10Knotw*/

void SetGeomWin10T2Knotw ( GO_BSplinePatch *go )
{
  if ( geom10win != geom10t2knotwin ) {
    geom10win = geom10t2knotwin;
    xge_SetMenuWidgets ( geom10menu, geom10t2knotwin, false );
    Geom10MenuResize ();
    geom10menu->redraw = xge_DrawMenu;
    bsp_sproduct_eqmer = NULL;
  }
  xge_T2KnotWindFindMapping ( &g10t2knotwin );
} /*SetGeomWin10T2Knotw*/

void SetGeomWin10Win2D ( GO_BSplineHole *go )
{
  if ( geom10win != geom10win2D ) {
    geom10win = geom10win2D;
    xge_SetMenuWidgets ( geom10menu, geom10win2D, false );
    Geom10MenuResize ();
    geom10menu->redraw = xge_DrawMenu;
    bsp_sproduct_eqmer = NULL;
  }
} /*SetGeomWin10Win2D*/

void SetGeomWin10SPrEqMer ( GO_BSplineCurve *eqmer )
{
  bsp_sproduct_eqmer = eqmer;
  if ( geom10win != geom10knotwineqmer ) {
    geom10win = geom10knotwineqmer;
    xge_SetMenuWidgets ( geom10menu, geom10knotwineqmer, false );
    Geom10MenuResize ();
    geom10menu->redraw = xge_DrawMenu;
  }
} /*SetGeomWin10SPrEqMer*/

