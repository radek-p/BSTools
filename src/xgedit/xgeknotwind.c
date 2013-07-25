
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2011                            */
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


void xge_DrawKnotWind ( xge_widget *er, boolean onscreen )
{
  xge_KnotWind *knw;

  knw = er->data0;
  xge_DrawGeomWinBackground ( er );
  if ( knw->display_coord && knw->xx >= 0 )
    xge_KnotWindDrawCursorPos ( knw );
  xge_KnotWindDrawAxis ( knw );
  xge_KnotWindDrawKnots ( knw );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*xge_DrawKnotWind*/

boolean xge_KnotWindMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_KnotWind *knw;
  double       scf;

  knw = er->data0;
  switch ( msg ) {
case xgemsg_NULL:
    return true;
case xgemsg_ENTERING:
    if ( knw->panning )
      xge_SetCurrentWindowCursor ( xgeCURSOR_FLEUR );
    else
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    knw->xx = x;
    if ( knw->display_coord )
      goto redraw;
    return true;
case xgemsg_EXITING:
    if ( knw->panning )
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    knw->xx = -1;
    if ( knw->display_coord )
      goto redraw;
    return true;
case xgemsg_MMOVE:
    knw->xx = x;
    break;
case xgemsg_SPECIAL_KEY:
    return false;
default:
    break;
  }

  switch ( er->state ) {
case xgestate_KNOTWIN_MOVINGKNOT:
    switch ( msg ) {
  case xgemsg_KEY:
      if ( key == 'r' || key == 'R' ) {
        if ( xge_KnotWindRemoveKnot ( knw ) ) {
          er->state = xgestate_NOTHING;
          xge_ReleaseFocus ( er );
          goto redraw;
        }
        break;
      }
      else
        return false;

  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_state;
      break;

  case xgemsg_MMOVE:
      if ( key & xgemouse_LBUTTON_DOWN ) {
        xge_BoundPoint ( er, &x, &y );
        if ( xge_KnotWindSetKnot ( knw, x ) )
          goto redraw;
        else
          goto exit_state;
      }
      else {
exit_state:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        goto redraw;
      }
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
        xge_KnotWindPan ( knw, x-xge_xx );
        goto change_knot_mapping;
      }
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
exit_state1:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
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
        xge_KnotWindZoom ( knw, scf );
change_knot_mapping:
        xge_xx = x;
        xge_yy = y;
        xge_callback ( er, xgemsg_KNOTWIN_CHANGE_MAPPING, 0, 0, 0 );
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
  case xgemsg_MCLICK:
      if ( (key & xgemouse_LBUTTON_DOWN) && (key & xgemouse_LBUTTON_CHANGE) ) {
        if ( knw->panning ) {
          xge_xx = x;
          xge_yy = y;
          er->state = xgestate_KNOTWIN_PANNING;
          xge_GrabFocus ( er, true );
        }
        else if ( knw->locked ) {
          knw->newknot = xge_KnotWindUnmapKnot ( knw, x );
          xge_callback ( knw->er, xgemsg_KNOTWIN_MCLICK, key, x, y );
        }
        else if ( xge_KnotWindFindNearestKnot ( knw, x, y ) ) {
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
        }
        else if ( knw->locked ) {
          knw->newknot = xge_KnotWindUnmapKnot ( knw, x );
          xge_callback ( knw->er, xgemsg_KNOTWIN_MCLICK, key, x, y );
        }
        else {
          if ( xge_KnotWindInsertKnot ( knw, x ) )
            goto redraw;
        }
      }
      break;

  case xgemsg_MMOVE:
      if ( knw->locked ) {
        knw->newknot = xge_KnotWindUnmapKnot ( knw, x );
        xge_callback ( knw->er, xgemsg_KNOTWIN_MMOVE, key, x, y );
      }
      if ( knw->display_coord )
        goto redraw;
      break;

  case xgemsg_RESIZE:
      er->w = x;
      er->h = y;
      xge_KnotWindInitMapping ( knw, knw->umin, knw->umax );
      if ( key )
        goto redraw;
      else
        break;

  case xgemsg_KEY:
      switch ( key ) {
    case 'f':  case 'F':
        xge_KnotWindFindMapping ( knw );
        goto redraw;
    case 'r':  case 'R':
        xge_KnotWindResetMapping ( knw );
        goto redraw;
    default:
        return false;
      }
      break;

  default:
      break;
    }
    break;

default:
    break;
  }
  return true;

redraw:
  xge_SetClipping ( er );
  er->redraw ( er, true );
  return true;
} /*xge_KnotWindMsg*/

void xge_KnotWindDrawCursorPos ( xge_KnotWind *knw )
{
  xge_widget *er;
  char       s[30];
  double     u;
  int        x;

  er = knw->er;
  u = xge_KnotWindUnmapKnot ( knw, knw->xx );
  sprintf ( s, "%5.3f", u );
  xgeSetForeground ( xgec_Green5 );
  xgeDrawLine ( knw->xx, er->y, knw->xx, er->y+er->h );
  xgeSetForeground ( xgec_Green );
  x = er->x+er->w-strlen(s)*6;
  x = min ( knw->xx+2, x );
  xgeDrawString ( s, x, er->y+er->h-4 );
} /*xge_KnotWindDrawCursorPos*/

void xge_KnotWindDrawAxis ( xge_KnotWind *knw )
{
  xge_widget *er;
  short      y;

  er = knw->er;
  y = (short)(er->y+4);
  xgeSetForeground ( xgec_Green );
  xgeDrawLine ( er->x+5, y, er->x+er->w-6, y );
} /*xge_KnotWindDrawAxis*/

void xge_KnotWindDrawKnots ( xge_KnotWind *knw )
{
  xge_widget *er;
  int        i;
  short      ux, uxa, y, uy;
  double     uu;

  if ( knw->lastknot < 3 )
    return;
  er = knw->er;
  y = (short)(er->y+4);
  xgeSetForeground ( xgec_White );
  ux = xge_KnotWindMapKnot ( knw, knw->knots[knw->degree] ); 
  uxa = xge_KnotWindMapKnot ( knw, knw->knots[knw->lastknot-knw->degree] );
  xgeDrawLine ( ux, y, uxa, y );
  uu = knw->knots[0]-1.0;
  for ( i = 1, uy = y; i < knw->lastknot; i++ ) {
    if ( knw->knots[i] > uu ) {
      uy = (short)(y+5);
      uu = knw->knots[i];
    }
    else
      uy += 5;
    ux = xge_KnotWindMapKnot ( knw, knw->knots[i] );
    xgeFillRectangle ( 3, 3, ux-1, uy-1 );
  }
} /*xge_KnotWindDrawKnots*/

void xge_KnotWindInitMapping ( xge_KnotWind *knw, double umin, double umax )
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
} /*xge_KnotWindInitMapping*/
 
void xge_KnotWindZoom ( xge_KnotWind *knw, double scf )
{
  knw->knotscf *= scf;
  if ( knw->knotscf < xge_KNOTWIN_MIN_SCALE )
    knw->knotscf = xge_KNOTWIN_MIN_SCALE;
  else if ( knw->knotscf > xge_KNOTWIN_MAX_SCALE )
    knw->knotscf = xge_KNOTWIN_MAX_SCALE;
} /*xge_KnotWindZoom*/

void xge_KnotWindPan ( xge_KnotWind *knw, int dxi )
{
  knw->knotshf += dxi/(knw->knotscf*knw->akm);
} /*xge_KnotWindPan*/

void xge_KnotWindFindMapping ( xge_KnotWind *knw )
{
  double a, b;

  knw->knotscf = 1.0;
  knw->knotshf = 0.0;
  a = knw->knots[1];
  b = knw->knots[knw->lastknot-1];
  if ( knw->altknots && knw->lastaltknot > 1 ) {
    a = min ( a, knw->altknots[1] );
    b = max ( b, knw->altknots[knw->lastaltknot-1] );
  }
  xge_KnotWindInitMapping ( knw, a, b );
} /*xge_KnotWindFindMapping*/

void xge_KnotWindResetMapping ( xge_KnotWind *knw )
{
  knw->knotscf = 1.0;
  knw->knotshf = 0.0;
  xge_KnotWindInitMapping ( knw, knw->umin, knw->umax );
} /*xge_KnotWindResetMapping*/

short xge_KnotWindMapKnot ( xge_KnotWind *knw, double u )
{
  return (short)((knw->knotscf*knw->akm)*(u+knw->knotshf)+knw->bkm+0.5);
} /*xge_KnotWindMapKnot*/

double xge_KnotWindUnmapKnot ( xge_KnotWind *knw, short xi )
{
  return ((double)xi-knw->bkm)/(knw->knotscf*knw->akm)-knw->knotshf;
} /*xge_KnotWindUnmapKnot*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean xge_KnotWindFindNearestKnot ( xge_KnotWind *knw, int x, int y )
{
  xge_widget *er;
  double     d, e;
  int        xi, i, r, s;

  er = knw->er;
  e = (double)(xge_MINDIST+1);
  i = 1;
  while ( i < knw->lastknot ) {
    xi = xge_KnotWindMapKnot ( knw, knw->knots[i] );
    d = (double)abs (xi-x);
    if ( d < e ) {
      knw->current_knot = i;
      e = d;
    }
    i++;
  }
  if ( e < xge_MINDIST ) {
    r = mbs_KnotMultiplicityd ( knw->lastknot-2, &knw->knots[1],
                                knw->knots[knw->current_knot] );
    s = max ( 1, (y-(er->y+4))/5 );
    while ( knw->current_knot < knw->lastknot-1 &&
            knw->knots[knw->current_knot] >= knw->knots[knw->current_knot+1] )
      knw->current_knot++;
    knw->current_mult = (unsigned char)(min ( r, s ));
    return true;
  }
  else
    return false;
} /*xge_KnotWindFindNearestKnot*/

boolean xge_KnotWindSetKnot ( xge_KnotWind *knw, short x )
{
  void    *sp;
  int     c, cc, r, i;
  double  u, du, _clcT, *newknots;

  sp = pkv_GetScratchMemTop ();
  if ( knw->locked )
    goto failure;
  newknots = pkv_GetScratchMem ( (knw->lastknot+1)*sizeof(double) );
  if ( !newknots )
    goto failure;
  memcpy ( newknots, knw->knots, (knw->lastknot+1)*sizeof(double) );

  u = xge_KnotWindUnmapKnot ( knw, x );
  r = mbs_FindKnotIntervald ( 1, knw->lastknot, newknots, u, NULL );
  if ( fabs(u-newknots[r]) < xge_KNOT_EPS )
    u = newknots[r];
  else if ( fabs(u-newknots[r+1]) < xge_KNOT_EPS )
    u = newknots[r+1];
  r = mbs_KnotMultiplicityd ( knw->lastknot-2, &newknots[1], u );
  c = cc = knw->current_knot;
  if ( r+knw->current_mult > knw->degree )
    goto success;
                             /* modify the knot sequence */
  du = u-knw->knots[knw->current_knot];
  if ( knw->closed ) {
    _clcT = knw->clcT;
    if ( knw->moving_many ) {
                   /* still to be written */
      if ( knw->current_knot == 0 )
        goto move_them;
      else if ( knw->current_knot <= knw->degree ) {
        for ( i = 0; i <= knw->lastknot; i++ )
          newknots[i] += du;
      }
      else if ( du > 0.0 || u > newknots[knw->current_knot-1] ) {
        _clcT += du;
        for ( i = knw->current_knot; i <= knw->lastknot; i++ )
          if ( i-knw->current_knot < knw->clcK )
            newknots[i] += du;
          else
            newknots[i] = newknots[i-knw->clcK]+_clcT;
        for ( i = knw->current_knot-1; i >= 0; i-- )
          if ( i+knw->clcK <= knw->lastknot )
            newknots[i] = newknots[i+knw->clcK]-_clcT;
      }
      else {
                   /* still to be written */
        goto success;
      }
    }
    else
      c = mbs_SetKnotClosedd ( knw->degree, knw->lastknot, newknots, knw->clcT,
                               knw->current_knot, knw->current_mult, u );
    if ( !mbs_ClosedKnotsCorrectd ( knw->degree, knw->lastknot, newknots,
                                    _clcT, knw->clcK, xge_KNOT_EPS ) )
      goto failure;
    else
      knw->clcT = _clcT;
  }
  else {
    if ( knw->moving_many ) {
      if ( du > 0.0 || knw->current_knot == 0 ) 
        goto move_them;
      else if ( u > newknots[knw->current_knot-1] ) {
move_them:
        for ( i = knw->current_knot; i <= knw->lastknot; i++ )
          newknots[i] += du;
      }
      else {
                   /* still to be written */
        goto success;
      }
    }
    else
      c = mbs_SetKnotd ( knw->lastknot, newknots,
                         knw->current_knot, knw->current_mult, u );
  }
  if ( c >= 0 && c <= knw->lastknot ) {
    memcpy ( knw->knots, newknots, (knw->lastknot+1)*sizeof(double) );
    knw->current_knot = c;
  }
  else
    goto failure;

success:
  pkv_SetScratchMemTop ( sp );
  if ( knw->switchkn )
    xge_callback ( knw->er, xgemsg_KNOTWIN_CHANGE_ALTKNOT, cc, 0, 0 );
  else
    xge_callback ( knw->er, xgemsg_KNOTWIN_CHANGE_KNOT, cc, 0, 0 );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*xge_KnotWindSetKnot*/

boolean xge_KnotWindInsertKnot ( xge_KnotWind *knw, short x )
{
  int    r;
  double u;

  if ( !knw->locked && knw->lastknot < knw->maxknots-1 ) {
    u = xge_KnotWindUnmapKnot ( knw, x );
    r = mbs_KnotMultiplicityd ( knw->lastknot, knw->knots, u );
    if ( u > knw->knots[knw->degree] &&
         u < knw->knots[knw->lastknot-knw->degree] &&
         r < knw->degree ) {
      knw->newknot = u;
      if ( knw->switchkn ) {
        if ( xge_callback ( knw->er, xgemsg_KNOTWIN_INSERT_ALTKNOT, 0, 0, 0 ) )
          return true;
      }
      else {
        if ( xge_callback ( knw->er, xgemsg_KNOTWIN_INSERT_KNOT, 0, 0, 0 ) )
          return true;
      }
    }
    else
      xge_callback ( knw->er, xgemsg_KNOTWIN_ERROR, 0, 0, 0 );
  }
  else
    xge_callback ( knw->er, xgemsg_KNOTWIN_ERROR, 1, 0, 0 );
  return false;
} /*xge_KnotWindInsertKnot*/

boolean xge_KnotWindRemoveKnot ( xge_KnotWind *knw )
{
  boolean correct;

  if ( !knw->locked ) {
    correct = knw->current_knot > knw->degree &&
         knw->current_knot < knw->lastknot-knw->degree &&
         ((!knw->closed && knw->lastknot >= 2*knw->degree+2) ||
           knw->lastknot >= 3*knw->degree+2);
    if ( knw->switchkn )
      return xge_callback ( knw->er, xgemsg_KNOTWIN_REMOVE_ALTKNOT,
                            (int)correct, 0, 0 );
    else {
      if ( xge_callback ( knw->er, xgemsg_KNOTWIN_REMOVE_KNOT,
                          (int)correct, 0, 0 ) ) {
        if ( knw->closed )
          knw->clcK = knw->lastknot - 2*knw->degree;
        return true;
      }
    }
  }
  else
    xge_callback ( knw->er, xgemsg_KNOTWIN_ERROR, 2, 0, 0 );
  return false;
} /*xge_KnotWindRemoveKnot*/

void xge_KnotWindSetAltKnots ( xge_KnotWind *knw,
               int altmaxkn, int lastaltkn, int altdeg, double *altknots ) 
{
  if ( altknots ) {
    knw->maxaltknots = altmaxkn;
    knw->lastaltknot = lastaltkn;
    knw->altdegree = altdeg;
    knw->altknots = altknots;
    knw->altkn = true;
  }
  else {
    knw->altknots = NULL;
    knw->altkn = knw->switchkn = false;
  }
} /*xge_KnotWindSetAltKnots*/

void xge_KnotWindSwitchAltKnots ( xge_KnotWind *knw )
{
#define SWAP(a,b,c) c = a,  a = b,  b = c;
  int    ti;
  double *tf;

  SWAP ( knw->maxknots, knw->maxaltknots, ti );
  SWAP ( knw->lastknot, knw->lastaltknot, ti );
  SWAP ( knw->degree, knw->altdegree, ti );
  SWAP ( knw->knots, knw->altknots, tf );
#undef SWAP
} /*xge_KnotWindSwitchAltKnots*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *xge_NewKnotWind ( char window_num, xge_widget *prev, int id,
                              short w, short h, short x, short y,
                              xge_KnotWind *knw, int maxknots, double *knots )
{
  xge_widget *er;

  er = xge_NewWidget ( window_num, prev, id, w, h, x, y,
                       knw, NULL, xge_KnotWindMsg, xge_DrawKnotWind );
  if ( er ) {
    knw->er = er;
    knw->panning = knw->display_coord = knw->moving_many = knw->locked =
    knw->closed = false;
    knw->knots = knots;
    knw->maxknots = maxknots;
    knw->lastknot = 0;
    knw->degree = 1;
    knw->knotscf = 1.0;
    knw->knotshf = 0.0;
    knw->altkn = knw->switchkn = false;
    knw->maxaltknots = knw->lastaltknot = knw->altdegree = 0;
    knw->altknots = NULL;
    xge_KnotWindInitMapping ( knw, 0.0, 1.0 );
  }
  return er;
} /*xge_NewKnotWind*/

