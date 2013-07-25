
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "splhole.h"

/* ////////////////////////////////////////////////////////////////////////// */
void DrawKnotWind ( xge_widget *er, boolean onscreen )
{
  ghKnotWind *knw;
  char i, j;
  short x0, x1, y;

  knw = er->data0;
  xge_DrawGeomWinBackground ( er );

  x0 = MapKnot ( knw, knw->umin );
  x0 = max ( x0, er->x);
  x1 = MapKnot ( knw, knw->umax );
  x1 = min ( x1, er->x+er->w);
  if ( x0 < x1 ) {
    xgeSetForeground ( xgec_Green2 );
    for ( i = 0; i < knw->hole_k; i++ ) {
      y = (short)(er->y+i*7+5);
      xgeDrawLine ( x0, y, x1, y );
    }
    xgeSetForeground ( xgec_White );
    for ( i = 0; i < knw->hole_k; i++ ) {
      y = (short)(er->y+i*7+4);
      for ( j = 0; j < 11; j++ ) {
        x0 = MapKnot ( knw, knots[11*i+j] );
        xgeFillRectangle ( 3, 3, x0-1, y );
      }
    }
  }
  if ( knw->display_coord )
    DrawKnotWindCursorPos ( knw );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*DrawKnotWind*/

void DrawKnotWindCursorPos ( ghKnotWind *knw )
{
  xge_widget *er;
  char       s[30];
  double     u;
  short      x;

  er = knw->er;
  u = UnmapKnot ( knw, knw->xx );
  sprintf ( s, "%5.3f", u );
  xgeSetForeground ( xgec_Green5 );
  xgeDrawLine ( knw->xx, er->y, knw->xx, er->y+er->h );
  xgeSetForeground ( xgec_Green );
  x = (short)(er->x+er->w-strlen(s)*6);
  x = min ( x, knw->xx+2 );
  xgeDrawString ( s, x, er->y+er->h-4 );
} /*DrawKnotWindCursorPos*/

boolean FindNearestKnot ( ghKnotWind *knw, short x, short y )
{
#define TOL 7
  xge_widget *er;
  char       i, j, k;
  short      xu, yu, d, dmin;

  j = k = -1;
  er = knw->er;
  dmin = TOL+1;
  for ( i = 0; i < knw->hole_k; i++ ) {
    yu = (short)(er->y+5+7*i);
    d = abs ( y-yu );
    if ( d < dmin ) {
      j = i;
      dmin = d;
    }
  }
  if ( dmin > TOL )
    return false;
  knw->current_seq = j;
  dmin = TOL+1;
  for ( i = 2; i <= 8; i++ ) {
    xu = MapKnot ( knw, knw->knots[11*j+i] );
    d = abs ( x-xu );
    if ( d < dmin ) {
      k = i;
      dmin = d;
    }
  }
  if ( dmin > TOL )
    return false;
  knw->current_knot = k;
  return true;
#undef TOL
} /*FindNearestKnot*/

boolean SetKnot ( ghKnotWind *knw, short x )
{
  double *knots, uk, vk;
  char hole_k, i, cs, css, ck, ckk;

  uk = UnmapKnot ( knw, x );
  vk = 1.0-uk;
  knots = knw->knots;
  hole_k = knw->hole_k;
  cs = knw->current_seq;
  ck = knw->current_knot;
  ckk = 10-ck;
  switch ( ck ) {
case 2:
case 3:
    css = (cs+hole_k-2) % hole_k;
    goto proceed;

case 7:
case 8:
    css = (cs+2) % hole_k;
proceed:
    if ( uk > knots[11*cs+ck-1] && uk < knots[11*cs+ck+1] &&
         vk > knots[11*css+ckk-1] && vk < knots[11*css+ckk+1] ) {
      knots[11*cs+ck] = uk;
      knots[11*css+ckk] = vk;
      CallBack ( knw->er, xgemsg_KNOTWIN_CHANGE_KNOT, 0, 0, 0 );
      return true;
    }
    break;

case 5:
    if ( hole_k % 4 == 0 ) {
      if ( uk > knots[11*cs+ck-1] && uk < knots[11*cs+ck+1] ) {
        css = (cs+2) % hole_k;
        for ( i = 0; i < hole_k; i += 4 ) {
          knots[11*cs+ck] = uk;
          knots[11*css+ckk] = vk;
          cs = (cs+4) % hole_k;
          css = (css+4) % hole_k;
        }
        CallBack ( knw->er, xgemsg_KNOTWIN_CHANGE_KNOT, 0, 0, 0 );
        return true;
      }
    }
    break;

case 4:
case 6:
    if ( hole_k & 0x01 ) {  /* odd */
      if ( uk > knots[11*cs+ck-1] && uk < knots[11*cs+ck+1] ) {
        for ( cs = 0; cs < hole_k; cs++ ) {
          knots[11*cs+ck] = uk;
          knots[11*cs+ckk] = vk;
        }
        CallBack ( knw->er, xgemsg_KNOTWIN_CHANGE_KNOT, 0, 0, 0 );
        return true;
      }
    }
    else {                  /* even */
      if ( uk > knots[11*cs+ck-1] && uk < knots[11*cs+ck+1] ) {
        for ( i = 0; i < hole_k; i += 2 ) {
          knots[11*cs+ck] = uk;
          knots[11*cs+ckk] = vk;
          cs = (cs+2) % hole_k;
        }
        CallBack ( knw->er, xgemsg_KNOTWIN_CHANGE_KNOT, 0, 0, 0 );
        return true;
      }
    }
    break;

default:
    break;
  }
  return false;
} /*SetKnot*/

boolean KnotWindMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  ghKnotWind *knw;
  double     scf;

  knw = er->data0;
  switch ( msg ) {
case xgemsg_NULL:
    return true;
case xgemsg_ENTERING:
    if ( knw->panning )
      xge_SetCurrentWindowCursor ( xgecursor[3] );
    else
      xge_SetCurrentWindowCursor ( None );
    knw->xx = x;
    if ( knw->display_coord )
      goto redraw;
    return true;
case xgemsg_EXITING:
    if ( knw->panning )
      xge_SetCurrentWindowCursor ( None );
    knw->xx = -1;
    if ( knw->display_coord )
      goto redraw;
    return true;
case xgemsg_MMOVE:
    knw->xx = x;
    break;
case xgemsg_SPECIAL_KEY:
    return false;
case xgemsg_RESIZE:
    er->w = x;
    er->h = y;
    KnotWindInitMapping ( knw, knw->umin, knw->umax );
    if ( key ) {
      xge_SetClipping ( er );
      er->redraw ( er, true );
    }
    return true;
default:
    break;
  }

  switch ( er->state ) {
case xgestate_KNOTWIN_MOVINGKNOT:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_state1;
      if ( !SetKnot ( knw, x ) )
        goto exit_state1;
      goto redraw;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_state1;
      break;
  default:
      break;
    }
    break;

case xgestate_KNOTWIN_PANNING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_state1;
      else {
        KnotWindPan ( knw, x-xge_xx );
        goto change_knot_mapping;
      }
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
exit_state1:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        goto redraw;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_KNOTWIN_ZOOMING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_state1;
      else {
        scf = 1.0+(double)(x-xge_xx)/(double)er->w;
        KnotWindZoom ( knw, scf );
change_knot_mapping:
        xge_xx = x;
        xge_yy = y;
        goto redraw;
      }
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_state1;
      break;
  default:
      break;
    }
    break;

case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( knw->display_coord )
        goto redraw;
      break;
  case xgemsg_MCLICK:
      if ( (key & xgemouse_LBUTTON_DOWN) && (key & xgemouse_LBUTTON_CHANGE) ) {
        if ( knw->panning ) {
          xge_xx = x;
          xge_yy = y;
          er->state = xgestate_KNOTWIN_PANNING;
          xge_GrabFocus ( er, true );
          goto redraw;
        }
        else if ( FindNearestKnot ( knw, x, y ) ) {
          er->state = xgestate_KNOTWIN_MOVINGKNOT;
          xge_GrabFocus ( er, true );
          goto redraw;
        }
      }
      else if ( (key & xgemouse_RBUTTON_DOWN) && (key & xgemouse_RBUTTON_CHANGE) ) {
        if ( knw->panning ) {
          xge_xx = x;
          xge_yy = y;
          er->state = xgestate_KNOTWIN_ZOOMING;
          xge_GrabFocus ( er, true );
          goto redraw;
        }
      }
      break;
  case xgemsg_KEY:
      switch ( key ) {
    case 'r':  case 'R':
    case 'f':  case 'F':
        knw->knotscf = 1.0;
        knw->knotshf = 0.0;
        goto redraw;
    default:
        break;
      }
      break;
  default:
      break;
    }
    break;

default:
    break;
  }
  return false;

redraw:
  xge_SetClipping ( er );
  er->redraw ( er, true );
  return true;
} /*KnotWindMsg*/

void KnotWindInitMapping ( ghKnotWind *knw, double umin, double umax )
{
  xge_widget *er;
  short      ximin, ximax;

  er = knw->er;
  ximin = (short)(er->x+5);
  ximax = (short)(er->x+er->w-6);
  knw->umin = umin;
  knw->umax = umax;
  knw->akm = ((double)(ximax-ximin))/(umax-umin);
  knw->bkm = (double)ximin - knw->akm*umin;
} /*KnotWindInitMapping*/

short MapKnot ( ghKnotWind *knw, double x )
{
  return (short)((knw->knotscf*knw->akm)*(x+knw->knotshf)+knw->bkm+0.5);
} /*MapKnot*/

double UnmapKnot ( ghKnotWind *knw, short xi )
{
  return ((double)xi-knw->bkm)/(knw->knotscf*knw->akm)-knw->knotshf;
} /*UnmapKnot*/

void KnotWindPan ( ghKnotWind *knw, short dxi )
{
  knw->knotshf += dxi/(knw->knotscf*knw->akm);
} /*KnotWindPan*/

void KnotWindZoom ( ghKnotWind *knw, double scf )
{
  knw->knotscf *= scf;
  if ( knw->knotscf < xge_KNOTWIN_MIN_SCALE )
    knw->knotscf = xge_KNOTWIN_MIN_SCALE;
  else if ( knw->knotscf > xge_KNOTWIN_MAX_SCALE )
    knw->knotscf = xge_KNOTWIN_MAX_SCALE;
} /*KnotWindZoom*/

void KnotWindSetHoleK ( ghKnotWind *knw, char hole_k, double *knots )
{
  xge_widget *er;

  er = knw->er;
  knw->hole_k = hole_k;
  knw->knots = knots;
  er->h = 7*hole_k+4;
} /*KnotindSetHoleK*/

xge_widget *NewGHKnotWind ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            char hole_k, ghKnotWind *knw, double *knots )
{
  xge_widget *er;

  er = xge_NewWidget ( window_num, prev, id, w, h, x, y,
                       knw, knots, KnotWindMsg, DrawKnotWind );
  if ( er ) {
    knw->er = er;
    knw->hole_k = hole_k;
    knw->knots = knots;
    knw->panning = knw->display_coord = 0;
    knw->current_seq = knw->current_knot = -1;
    knw->knotscf = 1.0;
    knw->knotshf = 0.0;
    KnotWindInitMapping ( knw, 0.0, 1.0 );
  }
  return er;
} /*NewGHKnotWind*/

