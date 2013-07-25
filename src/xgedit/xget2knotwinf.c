
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkgeom.h"
#include "multibs.h"

#include "xgedit.h"
#include "xgeprivate.h"


boolean xge_T2KnotWinfMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_T2KnotWinf *T2win;
  char           c;

  T2win = er->data0;
  switch ( msg ) {
case xgemsg_NULL: 
    return true;  

case xgemsg_ENTERING:
    if ( T2win->panning )
      xge_SetCurrentWindowCursor ( xgeCURSOR_FLEUR );
    else if ( T2win->selecting_mode )
      xge_SetCurrentWindowCursor ( xgeCURSOR_ARROW );
    else
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    T2win->inside = true;
    xge_xx = x;
    xge_yy = y;
    if ( T2win->display_coord )
      goto redraw_with_popups;
    return true;     

case xgemsg_EXITING: 
    if ( T2win->panning || T2win->selecting_mode )
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    T2win->inside = false;
    if ( T2win->display_coord ) {
redraw_with_popups:
      xge_SetClipping ( er );
      er->redraw ( er, false );
      xge_RedrawPopups ();
      xge_SetClipping ( er );
      xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
    }  
    return true;

case xgemsg_MMOVE:
    xge_xx = x;
    xge_yy = y;
    break;

case xgemsg_RESIZE:
    goto resize;   

case xgemsg_SPECIAL_KEY:
    return false;

default:
    break;
  }

  switch ( er->state ) {
case xgestate_T2KNOTWIN_PANNING:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      break;

  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      if ( xge_T2KnotWinfPan ( T2win, x, y ) )
        goto redraw_it;
      else
        break;

  default:
      break;
    }
    break;

case xgestate_T2KNOTWIN_ZOOMING:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_mode;
      break;

  case xgemsg_MMOVE:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_mode;
      if ( y != T2win->yy ) {
        xge_T2KnotWinfZoom ( T2win, y );
        goto redraw_it;
      }
      else  
        break;

  default:
      break;
    }
    break;

case xgestate_T2KNOTWIN_MOVINGKNOT_U:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      xge_BoundPoint ( er, &x, &y );
      if ( xge_T2KnotWinfSetKnotU ( T2win, x ) )
        goto redraw_it;
      else
        goto exit_mode;
      break;

  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      break;

  case xgemsg_KEY:
      if ( key == 'r' || key == 'R' ) {
        if ( xge_T2KnotWinfRemoveKnotU ( T2win ) ) {
          er->state = xgestate_NOTHING;
          xge_ReleaseFocus ( er );
          goto redraw_it;
        }
      }
      break;

  default:
      break;
    }
    break;

case xgestate_T2KNOTWIN_MOVINGKNOT_V:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      xge_BoundPoint ( er, &x, &y );
      if ( xge_T2KnotWinfSetKnotV ( T2win, y ) )
        goto redraw_it;
      else
        goto exit_mode;
      break;

  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
exit_mode:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        goto redraw_it;
      }
      break;

  case xgemsg_KEY:
      if ( key == 'r' || key == 'R' ) {
        if ( xge_T2KnotWinfRemoveKnotV ( T2win ) ) {
          er->state = xgestate_NOTHING;
          xge_ReleaseFocus ( er );
          goto redraw_it;
        }
      }
      break;

  default:
      break;
    }
    break;

case xgestate_T2KNOTWIN_SELECTING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      T2win->selection_rect.x1 = x;
      T2win->selection_rect.y1 = y;
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        xge_OrderSelectionRect ( &T2win->selection_rect );
        xge_callback ( er, xgemsg_T2KNOTWIN_SELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      goto redraw_it;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        T2win->selection_rect.x1 = x;
        T2win->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &T2win->selection_rect );
        xge_callback ( er, xgemsg_T2KNOTWIN_SELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_T2KNOTWIN_UNSELECTING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      T2win->selection_rect.x1 = x;
      T2win->selection_rect.y1 = y;
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
        T2win->selection_rect.x1 = x;
        T2win->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &T2win->selection_rect );
        xge_callback ( er, xgemsg_T2KNOTWIN_UNSELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      goto redraw_it;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
        T2win->selection_rect.x1 = x;
        T2win->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &T2win->selection_rect );
        xge_callback ( er, xgemsg_T2KNOTWIN_UNSELECT_POINTS, 0, x, y );
exit_select_mode:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        goto redraw_it;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( T2win->display_coord )
        goto redraw_it;
      break;

  case xgemsg_MCLICK:
      if ( T2win->panning ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_T2KNOTWIN_PANNING;
          xge_GrabFocus ( er, true );
          T2win->xx = x;
          T2win->yy = y;
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
             (key & xgemouse_RBUTTON_CHANGE) ) {
          er->state = xgestate_T2KNOTWIN_ZOOMING;
          xge_GrabFocus ( er, true );
          T2win->xx = T2win->yy = y;
        }
      }
      else if ( T2win->selecting_mode ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_T2KNOTWIN_SELECTING;
          goto continue_selecting_mode;
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          er->state = xgestate_T2KNOTWIN_UNSELECTING;
continue_selecting_mode:
          xge_GrabFocus ( er, true );
          T2win->xx = T2win->selection_rect.x0 = T2win->selection_rect.x1 = x;
          T2win->yy = T2win->selection_rect.y0 = T2win->selection_rect.y1 = y;
          goto redraw_it;
        }
      }
      else {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          c = xge_T2KnotWinfFindNearestKnot ( T2win, x, y );
          switch ( c ) {
        case 1:
            er->state = xgestate_T2KNOTWIN_MOVINGKNOT_U;
            xge_GrabFocus ( er, true );
            goto redraw_it;
            break;
        case 2:
            er->state = xgestate_T2KNOTWIN_MOVINGKNOT_V;
            xge_GrabFocus ( er, true );
            goto redraw_it;
            break;
        default:
            break;
          }
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          switch ( xge_T2KnotWinfFindDomWinRegion ( T2win, x, y ) ) {
        case 1:
            if ( xge_T2KnotWinfInsertKnotU ( T2win, x ) )
              goto redraw_it;
            break;
        case 2:
            if ( xge_T2KnotWinfInsertKnotV ( T2win, y ) )
              goto redraw_it;
            break;
       default:
            break;
          }
        }
      }
      break;

  case xgemsg_KEY:
      switch ( key ) {
    case 'f':  case 'F':
        xge_T2KnotWinfFindMapping ( T2win );
        goto redraw_it;
    case 'r':  case 'R':
        xge_T2KnotWinfResetMapping ( T2win );
        goto redraw_it;
    default:
        break;
      }
      break;

  case xgemsg_RESIZE:
resize:
      er->w = x;
      er->h = y;
      xge_T2KnotWinfSetupMapping ( T2win );
redraw_it:
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;

  default:
      break;
    }

default:
    break;
  }
  return true;
} /*xge_T2KnotWinfMsg*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_T2KnotWinfSetupMapping ( xge_T2KnotWinf *T2win )
{
  xge_widget *er;
  float      rbw, rbh;

  er = T2win->er;
  CameraInitFramef ( &T2win->CPos, true, true,
                     (short)(er->w-T2win->knot_margin),
                     (short)(er->h-T2win->knot_margin),
                     (short)(er->x+T2win->knot_margin), er->y, xge_aspect, 6 );
  T2win->CPos.position = T2win->centre;
  T2win->CPos.psi = T2win->CPos.theta = T2win->CPos.phi = 0.0;
  rbw = T2win->RefBBox.x1-T2win->RefBBox.x0;
  rbh = T2win->RefBBox.y1-T2win->RefBBox.y0;
  if ( T2win->CPos.width*rbh > T2win->CPos.height*rbw ) {
    T2win->CPos.vd.para.hgh = rbh;
    T2win->CPos.vd.para.dim_case = 2;
  }
  else {
    T2win->CPos.vd.para.wdt = rbw;
    T2win->CPos.vd.para.dim_case = 1;
  }
  CameraSetMappingf ( &T2win->CPos );
} /*xge_T2KnotWinfSetupMapping*/

void xge_T2KnotWinfInitMapping ( xge_T2KnotWinf *T2win,
                                 float umin, float umax, float vmin, float vmax )
{
  T2win->RefBBox.x0 = umin;  T2win->RefBBox.x1 = umax;
  T2win->RefBBox.y0 = vmin;  T2win->RefBBox.y1 = vmax;
  T2win->DefBBox = T2win->RefBBox;
  SetPoint3f ( &T2win->centre,
               (float)(0.5*(umin+umax)), (float)(0.5*(vmin+vmax)), 0.0 );
  xge_T2KnotWinfSetupMapping ( T2win );
} /*xge_T2KnotWinfInitMapping*/

void xge_T2KnotWinfFindMapping ( xge_T2KnotWinf *T2win )
{
  float umin, umax, vmin, vmax;

  if ( T2win->altknots_u && T2win->altlastknot_u > 1 ) {
    umin = min ( T2win->knots_u[1], T2win->altknots_u[1] );
    umax = max ( T2win->knots_u[T2win->lastknot_u-1],
                 T2win->altknots_u[T2win->altlastknot_u-1]);
  }
  else {
    umin = T2win->knots_u[1];
    umax = T2win->knots_u[T2win->lastknot_u-1];
  }
  if ( T2win->altknots_v && T2win->altlastknot_v > 1 ) {
    vmin = min ( T2win->knots_v[1], T2win->altknots_v[1] );
    vmax = max ( T2win->knots_v[T2win->lastknot_v-1],
                 T2win->altknots_v[T2win->altlastknot_v-1]);
  }
  else {
    vmin = T2win->knots_v[1];
    vmax = T2win->knots_v[T2win->lastknot_v-1];
  }
      /* extend by 5% */
  T2win->RefBBox.x0 = (float)(1.05*umin-0.05*umax);
  T2win->RefBBox.x1 = (float)(-0.05*umin+1.05*umax);
  T2win->RefBBox.y0 = (float)(1.05*vmin-0.05*vmax);
  T2win->RefBBox.y1 = (float)(-0.05*vmin+1.05*vmax);
  SetPoint3f ( &T2win->centre,
               (float)(0.5*(umin+umax)), (float)(0.5*(vmin+vmax)), 0.0 );
  xge_T2KnotWinfSetupMapping ( T2win );
} /*xge_T2KnotWinfFindMapping*/

void xge_T2KnotWinfResetMapping ( xge_T2KnotWinf *T2win )
{
  T2win->RefBBox = T2win->DefBBox;
  SetPoint3f ( &T2win->centre,
               (float)(0.5*(T2win->RefBBox.x0+T2win->RefBBox.x1)),
               (float)(0.5*(T2win->RefBBox.y0+T2win->RefBBox.y1)), 0.0 );
  xge_T2KnotWinfSetupMapping ( T2win );
} /*xge_T2KnotWinfResetMapping*/

void xge_T2KnotWinfZoom ( xge_T2KnotWinf *T2win, short y )
{
  xge_widget *er;
  float      f, w, h;

  er = T2win->er;
  f = (float)(y-T2win->yy)/(float)er->h;
  f = exp ( f );
  f = (float)max ( f, xge_KNOTWIN_MIN_SCALE );
  f = (float)min ( f, xge_KNOTWIN_MAX_SCALE );
  SetPoint3f ( &T2win->centre,
               (float)(0.5*(T2win->RefBBox.x0+T2win->RefBBox.x1)),
               (float)(0.5*(T2win->RefBBox.y0+T2win->RefBBox.y1)), 0.0 );
  w = 0.5*f*(T2win->RefBBox.x1-T2win->RefBBox.x0);
  h = 0.5*f*(T2win->RefBBox.y1-T2win->RefBBox.y0);
  T2win->RefBBox.x0 = T2win->centre.x-w;
  T2win->RefBBox.x1 = T2win->centre.x+w;
  T2win->RefBBox.y0 = T2win->centre.y-h;
  T2win->RefBBox.y1 = T2win->centre.y+h;
  xge_T2KnotWinfSetupMapping ( T2win );
  xge_callback ( er, xgemsg_T2KNOTWIN_PROJCHANGE, 0, 0, y );
  T2win->yy = y;
} /*xge_T2KnotWinfZoom*/

boolean xge_T2KnotWinfPan ( xge_T2KnotWinf *T2win, short x, short y )
{
  xge_widget *er;
  point2f    a, b;
  vector2f   shift;

  er = T2win->er;
  if ( x != T2win->xx || y != T2win->yy ) {
    SetPoint2f ( &a, (float)T2win->xx, (float)T2win->yy );
    CameraUnProjectPoint2f ( &T2win->CPos, &a, &a ); 
    SetPoint2f ( &b, (float)x, (float)y );
    CameraUnProjectPoint2f ( &T2win->CPos, &b, &b );
    SubtractPoints2f ( &a, &b, &shift );
    AddVector2f ( (point2f*)(void*)&T2win->centre.x, &shift,
                  (point2f*)(void*)&T2win->centre.x );
/* bound the centre --- to be added */
    a.x = (float)(0.5*(T2win->RefBBox.x1-T2win->RefBBox.x0));
    a.y = (float)(0.5*(T2win->RefBBox.y1-T2win->RefBBox.y0));
    T2win->RefBBox.x0 = T2win->centre.x-a.x;
    T2win->RefBBox.x1 = T2win->centre.x+a.x;
    T2win->RefBBox.y0 = T2win->centre.y-a.y;
    T2win->RefBBox.y1 = T2win->centre.y+a.y;
    xge_T2KnotWinfSetupMapping ( T2win );
    xge_callback ( er, xgemsg_T2KNOTWIN_PROJCHANGE, 0, x, y );
    T2win->xx = x;
    T2win->yy = y;
    return true;
  }
  else
    return false;
} /*xge_T2KnotWinfPan*/

/* ///////////////////////////////////////////////////////////////////////// */
char xge_T2KnotWinfFindDomWinRegion ( xge_T2KnotWinf *T2win, int x, int y )
{
  int  eta;
  char regcode;

  eta = T2win->CPos.ymin+T2win->CPos.height;
  if ( x >= T2win->CPos.xmin ) regcode = 1;  else regcode = 0;
  if ( y < eta ) regcode += 2;
  return regcode;
} /*xge_T2KnotWinfFindDomWinRegion*/

char xge_T2KnotWinfFindNearestKnot ( xge_T2KnotWinf *T2win, int x, int y )
{
  /* return value 0 - knot not found near x, y */
  /* 1 - found knot in the "u" sequence */
  /* 2 - found knot in the "v" sequence */
  char  regcode;
  int   i, r, s, d, e;
  short xi, eta;

  regcode = xge_T2KnotWinfFindDomWinRegion ( T2win, x, y );
  switch ( regcode ) {
case 1:
    if ( T2win->locked_u )
      return 0;
    /* search the "u" knot sequence */
    e = xge_MINDIST;
    for ( i = 1; i < T2win->lastknot_u; i++ ) {
      xi = xge_T2KnotWinfMapKnotU ( T2win, T2win->knots_u[i] );
      if ( (d = abs ( xi - x )) < e ) {
        T2win->current_item = i;
        e = d;
      }
    }
    if ( e < xge_MINDIST ) {
      eta = (short)(T2win->CPos.ymin+T2win->CPos.height);
      r = mbs_KnotMultiplicityf ( T2win->lastknot_u-2, &T2win->knots_u[1],
                                  T2win->knots_u[T2win->current_item] );
      s = max ( 0, (y-eta-5)/5 ) + 1;
      while ( T2win->current_item < T2win->lastknot_u-1 &&
              T2win->knots_u[T2win->current_item] ==
              T2win->knots_u[T2win->current_item+1] )
        T2win->current_item ++;
      T2win->current_mult = (unsigned char)(min ( r, s ));
      return 1;
    }
    break;

case 2:
    if ( T2win->locked_v )
      return 0;
    /* search the "v" knot sequence */
    e = xge_MINDIST;
    for ( i = 1; i < T2win->lastknot_v; i++ ) {
      eta = xge_T2KnotWinfMapKnotV ( T2win, T2win->knots_v[i] );
      if ( (d = abs ( eta - y )) < e ) {
        T2win->current_item = i;
        e = d;
      }
    }
    if ( e < xge_MINDIST ) {
      r = mbs_KnotMultiplicityf ( T2win->lastknot_v-2, &T2win->knots_v[1],
                                  T2win->knots_v[T2win->current_item] );
      s = max ( 0, (T2win->CPos.xmin-5-x)/5 ) + 1;
      while ( T2win->current_item < T2win->lastknot_v-1 &&
              T2win->knots_v[T2win->current_item] ==
              T2win->knots_v[T2win->current_item+1] )
        T2win->current_item ++;
      T2win->current_mult = (unsigned char)(min ( r, s ));
      return 2;
    }
    break;

default:
    break;
  }
  return 0;
} /*xge_T2KnotWinfFindNearestKnot*/

short xge_T2KnotWinfMapKnotU ( xge_T2KnotWinf *T2win, float u )
{
  point2f p, q;

  SetPoint2f ( &p, u, 0.0 );
  CameraProjectPoint2f ( &T2win->CPos, &p, &q );
  return (short)(q.x+0.5);
} /*xge_T2KnotWinfMapKnotU*/

float xge_T2KnotWinfUnmapKnotU ( xge_T2KnotWinf *T2win, short xi )
{
  point2f p, q;

  SetPoint2f ( &q, (float)xi, 0.0 );
  CameraUnProjectPoint2f ( &T2win->CPos, &q, &p );
  return p.x;
} /*xge_T2KnotWinfUnmapKnotU*/

short xge_T2KnotWinfMapKnotV ( xge_T2KnotWinf *T2win, float v )
{
  point2f p, q;

  SetPoint2f ( &p, 0.0, v );
  CameraProjectPoint2f ( &T2win->CPos, &p, &q );
  return (short)(q.y+0.5);
} /*xge_T2KnotWinfMapKnotV*/

float xge_T2KnotWinfUnmapKnotV ( xge_T2KnotWinf *T2win, short eta )
{
  point2f p, q;

  SetPoint2f ( &q, 0.0, (float)eta );
  CameraUnProjectPoint2f ( &T2win->CPos, &q, &p );
  return p.y;
} /*xge_T2KnotWinfUnmapKnotV*/

static boolean _xge_T2KnotWinfSetKnot ( xge_T2KnotWinf *T2win, boolean uu,
                      short xy, int degree, int lastknot, float *knots,
                      boolean closed, float *clcT,
                      float (*unmapknot)(xge_T2KnotWinf *T2win, short xy) )
{
  void    *sp;
  float   t, dt, *newknots, _clcT;
  int     i, r, c, cc, clcK;

  sp = pkv_GetScratchMemTop ();
  if ( !(newknots = pkv_GetScratchMemf ( lastknot+1 )) )
    goto failure;
  memcpy ( newknots, knots, (lastknot+1)*sizeof(float) );

  t = unmapknot ( T2win, xy );
  r = mbs_FindKnotIntervalf ( 1, lastknot-2, newknots, t, NULL );
  if ( fabs(t-newknots[r]) < xge_KNOT_EPS )
    t = newknots[r];
  else if ( fabs(t-newknots[r+1]) < xge_KNOT_EPS )
    t = newknots[r+1];
  r = mbs_KnotMultiplicityf ( lastknot-2, &newknots[1], t );
  c = cc = T2win->current_item;
  if ( r+T2win->current_mult > degree )
    goto success;

  dt = t - knots[c];
  if ( closed ) {
    _clcT = *clcT;
    clcK = lastknot-2*degree;
    if ( T2win->moving_many ) {
                   /* still to be written */
      if ( c == 0 )
        goto move_them;
      else if ( c <= degree ) {
        for ( i = 0; i <= lastknot; i++ )
          newknots[i] += dt;
      }
      else if ( dt > 0.0 || t > newknots[c-1] ) {
        _clcT += dt;
        for ( i = c; i <= lastknot; i++ )
          if ( i-c < clcK )
            newknots[i] += dt;
          else
            newknots[i] = newknots[i-clcK]+_clcT;
        for ( i = c-1; i >= 0; i-- )  
          if ( i+clcK <= lastknot )
            newknots[i] = newknots[i+clcK]-_clcT;
      }
      else {
                   /* still to be written */
        goto success;
      }
    }
    else
      c = mbs_SetKnotClosedf ( degree, lastknot, newknots, _clcT,
                               c, T2win->current_mult, t );
    if ( !mbs_ClosedKnotsCorrectf ( degree, lastknot, newknots,
                                    _clcT, clcK, xge_KNOT_EPS ) )
      goto failure;
    else
      *clcT = _clcT;
  }
  else {
    if ( T2win->moving_many ) {
      if ( dt > 0.0 || c == 0 )
        goto move_them;
      else if ( t > newknots[c-1] ) {
move_them:
        for ( i = c; i <= lastknot; i++ )
          newknots[i] += dt;
      }
      else {
                   /* still to be written */
        goto success;
      }
    }
    else
      c = mbs_SetKnotf ( lastknot, newknots, c, T2win->current_mult, t );
  }
  if ( c >= 0 && c <= lastknot ) {
    memcpy ( knots, newknots, (lastknot+1)*sizeof(float) );
    T2win->current_item = c;
  }
  else
    goto failure;

success:
  pkv_SetScratchMemTop ( sp );
  if ( uu ) {
    if ( T2win->switchknu )
      xge_callback ( T2win->er, xgemsg_T2KNOTWIN_CHANGE_ALTKNOT_U, cc, 0, 0 );
    else
      xge_callback ( T2win->er, xgemsg_T2KNOTWIN_CHANGE_KNOT_U, cc, 0, 0 );
  }
  else {
    if ( T2win->switchknv )
      xge_callback ( T2win->er, xgemsg_T2KNOTWIN_CHANGE_ALTKNOT_V, cc, 0, 0 );
    else
      xge_callback ( T2win->er, xgemsg_T2KNOTWIN_CHANGE_KNOT_V, cc, 0, 0 );
  }
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_xge_T2KnotWinfSetKnot*/

boolean xge_T2KnotWinfSetKnotU ( xge_T2KnotWinf *T2win, short x )
{
  return _xge_T2KnotWinfSetKnot ( T2win, true, x,
                        T2win->degree_u, T2win->lastknot_u,
                        T2win->knots_u, T2win->closed_u, &T2win->clcTu,
                        xge_T2KnotWinfUnmapKnotU );
} /*xge_T2KnotWinfSetKnotU*/

boolean xge_T2KnotWinfSetKnotV ( xge_T2KnotWinf *T2win, short y )
{
  return _xge_T2KnotWinfSetKnot ( T2win, false, y,
                        T2win->degree_v, T2win->lastknot_v,
                        T2win->knots_v, T2win->closed_v, &T2win->clcTv,
                        xge_T2KnotWinfUnmapKnotV );
} /*xge_T2KnotWinfSetKnotV*/

boolean xge_T2KnotWinfInsertKnotU ( xge_T2KnotWinf *T2win, short x )
{
  int     r;
  float   u;
  boolean correct;

  if ( T2win->lastknot_u < T2win->maxknots_u-1 ) {
    u = xge_T2KnotWinfUnmapKnotU ( T2win, x );
    r = mbs_KnotMultiplicityf ( T2win->lastknot_u, T2win->knots_u, u );
    correct = u > T2win->knots_u[T2win->degree_u] &&
              u < T2win->knots_u[T2win->lastknot_u-T2win->degree_u] &&
              r < T2win->degree_u;
    T2win->newknot = u;
    if ( T2win->switchknu ) {
      if ( xge_callback ( T2win->er, xgemsg_T2KNOTWIN_INSERT_ALTKNOT_U,
                          (int)correct, 0, 0 ) )
        return true;
    }
    else {
      if ( xge_callback ( T2win->er, xgemsg_T2KNOTWIN_INSERT_KNOT_U,
                          (int)correct, 0, 0 ) )
        return true;
    }
  }
  else
    xge_callback ( T2win->er, xgemsg_KNOTWIN_ERROR, 1, 0, 0 );
  return false;
} /*xge_T2KnotWinfInsertKnotU*/

boolean xge_T2KnotWinfInsertKnotV ( xge_T2KnotWinf *T2win, short y )
{
  int     r;
  float   v;
  boolean correct;

  if ( T2win->lastknot_v < T2win->maxknots_v-1 ) {
    v = xge_T2KnotWinfUnmapKnotV ( T2win, y );
    r = mbs_KnotMultiplicityf ( T2win->lastknot_v, T2win->knots_v, v );
    correct = v > T2win->knots_v[T2win->degree_v] &&
              v < T2win->knots_v[T2win->lastknot_v-T2win->degree_v] &&
              r < T2win->degree_v;
    T2win->newknot = v;
    if ( T2win->switchknv ) {
      if ( xge_callback ( T2win->er, xgemsg_T2KNOTWIN_INSERT_ALTKNOT_V,
                          (int)correct, 0, 0 ) )
        return true;
    }
    else {
      if ( xge_callback ( T2win->er, xgemsg_T2KNOTWIN_INSERT_KNOT_V,
                          (int)correct, 0, 0 ) )
        return true;
    }
  }
  else
    xge_callback ( T2win->er, xgemsg_KNOTWIN_ERROR, 3, 0, 0 );
  return false;
} /*xge_T2KnotWinfInsertKnotV*/

boolean xge_T2KnotWinfRemoveKnotU ( xge_T2KnotWinf *T2win )
{
  boolean correct;

  correct = T2win->current_item > T2win->degree_u &&
            T2win->current_item < T2win->lastknot_u-T2win->degree_u &&
            ((!T2win->closed_u && T2win->lastknot_u >= 2*T2win->degree_u+2) ||
              T2win->lastknot_u >= 3*T2win->degree_u+2);
  if ( T2win->switchknu ) {
    return xge_callback ( T2win->er, xgemsg_T2KNOTWIN_REMOVE_ALTKNOT_U,
                          (int)correct, 0, 0 );
  }
  else if ( xge_callback ( T2win->er, xgemsg_T2KNOTWIN_REMOVE_KNOT_U,
                           (int)correct, 0, 0 ) ) {
    if ( T2win->closed_u )
      T2win->clcKu = T2win->lastknot_u - 2*T2win->degree_u;
    return true;
  }
  else
    return false;
} /*xge_T2KnotWinfRemoveKnotU*/

boolean xge_T2KnotWinfRemoveKnotV ( xge_T2KnotWinf *T2win )
{
  boolean correct;

  correct = T2win->current_item > T2win->degree_v &&
            T2win->current_item < T2win->lastknot_v-T2win->degree_v &&
            ((!T2win->closed_v && T2win->lastknot_v >= 2*T2win->degree_v+2) ||
              T2win->lastknot_v >= 3*T2win->degree_v+2);
  if ( T2win->switchknv ) {
    return xge_callback ( T2win->er, xgemsg_T2KNOTWIN_REMOVE_ALTKNOT_V,
                          (int)correct, 0, 0 );
  }
  else if ( xge_callback ( T2win->er, xgemsg_T2KNOTWIN_REMOVE_KNOT_V,
                           (int)correct, 0, 0 ) ) {
    if ( T2win->closed_v )
      T2win->clcKv = T2win->lastknot_v - 2*T2win->degree_v;
    return true;
  }
  else
    return false;
} /*xge_T2KnotWinfRemoveKnotV*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_T2KnotWinfDrawKnots ( xge_T2KnotWinf *T2win )
{
} /*xge_T2KnotWinfDrawKnots*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_T2KnotWinfSetAltKnots ( xge_T2KnotWinf *T2win,
               int altmaxknu, int lastaltknu, int altdegu, float *altknotsu,
               int altmaxknv, int lastaltknv, int altdegv, float *altknotsv )
{
  if ( altknotsu ) {
    T2win->altmaxknots_u = altmaxknu;
    T2win->altlastknot_u = lastaltknu;
    T2win->altknots_u = altknotsu;
    T2win->altdeg_u = altdegu;
    T2win->altknu = true;
  }
  else {
    T2win->altknots_u = NULL;
    T2win->altknu = T2win->switchknu = false;
  }
  if ( altknotsv ) {
    T2win->altmaxknots_v = altmaxknv;
    T2win->altlastknot_v = lastaltknv;
    T2win->altknots_v = altknotsv;
    T2win->altdeg_v = altdegv;
    T2win->altknv = true;
  }
  else {
    T2win->altknots_v = NULL;
    T2win->altknv = T2win->switchknu = false;
  }
} /*xge_T2KnotWinfSetAltKnots*/

void xge_T2KnotWinfSwitchAltKnots ( xge_T2KnotWinf *T2win,
               boolean altu, boolean altv )
{
#define SWAP(a,b,c) c = a,  a = b,  b = c;
  int   ti;
  float *tf;

  if ( altu != T2win->switchknu ) {
    SWAP ( T2win->maxknots_u, T2win->altmaxknots_u, ti )
    SWAP ( T2win->lastknot_u, T2win->altlastknot_u, ti )
    SWAP ( T2win->degree_u, T2win->altdeg_u, ti )
    SWAP ( T2win->knots_u, T2win->altknots_u, tf )
    T2win->switchknu = altu;
  }
  if ( altv != T2win->switchknv ) {
    SWAP ( T2win->maxknots_v, T2win->altmaxknots_v, ti )
    SWAP ( T2win->lastknot_v, T2win->altlastknot_v, ti )
    SWAP ( T2win->degree_v, T2win->altdeg_v, ti )
    SWAP ( T2win->knots_v, T2win->altknots_v, tf )
    T2win->switchknu = altv;
  }
#undef SWAP
} /*xge_T2KnotWinfSwitchAltKnots*/
 
/* ///////////////////////////////////////////////////////////////////////// */
void xge_T2KnotWinfDrawCursorPos ( xge_T2KnotWinf *T2win, short x, short y )
{
  xge_widget *er;
  char       s[10];
  point2f    p, q;
  int        xx;

  if ( T2win->display_coord ) {
    er = T2win->er;
    SetPoint2f ( &p, x, y );
    CameraUnProjectPoint2f ( &T2win->CPos, &p, &q );
    xgeSetForeground ( xgec_Green5 );
    xgeDrawLine ( x, er->y, x, er->y+er->h );
    xgeDrawLine ( er->x, y, er->x+er->w, y );
    xgeSetForeground ( xgec_Green );
    sprintf ( s, "%5.3f", q.x );
    xx = er->x+er->w-strlen(s)*6;
    x = (short)(min ( x+2, xx ));
    xgeDrawString ( s, x, er->y+er->h-4 );
    sprintf ( s, "%5.3f", q.y );
    y = (short)(max ( y, er->y+10 ));
    xgeDrawString ( s, er->x, y );
  }
} /*xge_T2KnotWinfDrawCursorPos*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *xge_NewT2KnotWinf ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                short knot_margin,
                                xge_T2KnotWinf *T2win,
                                void (*redraw)(xge_widget*, boolean),
                                int maxknots_u, float *knots_u,
                                int maxknots_v, float *knots_v )
{
  xge_widget *er;

  er = xge_NewWidget ( window_num, prev, id, w, h, x, y,
                       T2win, NULL, xge_T2KnotWinfMsg, redraw );
  if ( er ) {
    T2win->er = er;
    T2win->panning = T2win->selecting_mode = T2win->display_coord = false;
    T2win->closed_u = T2win->closed_v = false;
    T2win->locked_u = T2win->locked_v = false;
    T2win->knot_margin = knot_margin;
    T2win->maxknots_u = maxknots_u;
    T2win->lastknot_u = 1;
    T2win->degree_u   = 1;
    T2win->knots_u = knots_u;
    T2win->maxknots_v = maxknots_v;
    T2win->lastknot_v = 1;
    T2win->degree_v   = 1;
    T2win->knots_v = knots_v;
    T2win->altknu = T2win->switchknu = T2win->altknv = T2win->switchknv = false;
    T2win->altmaxknots_u = T2win->altlastknot_u = T2win->altdeg_u = 0;
    T2win->altknots_u = NULL;
    T2win->altmaxknots_v = T2win->altlastknot_v = T2win->altdeg_v = 0;
    T2win->altknots_v = NULL;
    xge_T2KnotWinfInitMapping ( T2win, -0.05, 1.05, -0.05, 1.05 );
  }
  return er;
} /*xge_NewT2KnotWinf*/

