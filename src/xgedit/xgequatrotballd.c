
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
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

#include "xgedit.h"
#include "xgeprivate.h"


void xge_DrawQuatRotBalld ( xge_widget *er, boolean onscreen )
{
  xge_quatrotballd *qball;
  short            xc, yc, x, y, r1, r2;
  char             *title;
  trans3f          tr;

  qball = (xge_quatrotballd*)er->data0;
  xgeSetForeground ( xgec_MENU_BACKGROUND );
  xgeFillRectangle ( er->x+er->w-1, er->y+er->h-1, er->x, er->y );
  xc = qball->xc;
  yc = qball->yc;
  r1 = qball->r1;
  r2 = qball->r2;
  xgeSetForeground ( xgec_Black );
  xgeFillArc ( 2*r2+1, 2*r2+1, xc-r2, yc-r2, 0, 360*64 );
  if ( er->state == xgestate_NOTHING )
    xgeSetForeground ( xgec_Blue3 );
  else
    xgeSetForeground ( xgec_Blue6 );
  xgeFillArc ( 2*r1+1, 2*r1+1, xc-r1, yc-r1, 0, 360*64 );
  xgeSetForeground ( xgec_Grey3 );
  xgeDrawArc ( 2*r1+1, 2*r1+1, xc-r1, yc-r1, 0, 360*64 );
  xgeSetForeground ( xgec_White );
  Trans3dTof ( qball->tr, &tr );
  _xge_QuatRotBallDrawCircles ( xc, yc, r1, &tr );
  xgeDrawArc ( 2*r2+1, 2*r2+1, xc-r2, yc-r2, 0, 360*64 );
  if ( (title = er->data1) ) {
    if ( er->w > er->h )  /* title aside */
      { x = (short)(xc+r2+4);  y = (short)(yc+4); }
    else                  /* title below */
      { x = (short)(xc-6*strlen(title)/2+1 );  y = (short)(yc+r2+13); }
    xgeDrawString ( title, x, y );
  }

  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawQuatRotBalld*/

static boolean xge_QuatRotBalldRotateIt ( xge_widget *er, short x, short y )
{
  xge_quatrotballd *qball;
  short            dx, dy, xc, yc, R;
  vector3d         a, b, v;
  double           d, e, phi, s, c;
  quaterniond      q;

  dx = x - xge_xx;
  dy = y - xge_yy;
  if ( !dx && !dy )
    return false;
  qball = (xge_quatrotballd*)er->data0;
  R  = qball->R;
  xc = qball->xc;
  yc = qball->yc;
  switch ( er->state ) {
case xgestate_QUATROT_TURNING1:
    SetVector3d ( &a, (double)(x-xc)/(double)R,
                      (double)(yc-y)/(double)R, 0.0 );
    break;
case xgestate_QUATROT_TURNING2:
    SetVector3d ( &a, 0.0, 0.0, 0.0 );
    break;
case xgestate_QUATROT_TURNING3:
    if ( !qball->axis_ok ) {
      SetVector3d ( &qball->axis, (double)(x-xc)/(double)R,
                        (double)(yc-y)/(double)R, 0.0 );
      qball->axis_ok = true;
    }
    a = qball->axis;
    break;
default:
    return false;
  }
  SetVector3d ( &b, (double)(x-xge_xx)/(double)R,
                    (double)(xge_yy-y)/(double)R, 0.0 );
  OrtVector3d ( &b, &a, &a );
  d = a.x*a.x + a.y*a.y;
  e = b.x*b.x + b.y*b.y;
  if ( d >= 1.0 ) {
    phi = atan2 ( sqrt(e), sqrt(d) );
    NormalizeVector3d ( &a );
    if ( a.x*b.y - a.y*b.x >= 0.0 )
      SetVector3d ( &v, 0.0, 0.0, 1.0 );
    else
      SetVector3d ( &v, 0.0, 0.0, -1.0 );
  }
  else {
    a.z = sqrt ( 1.0 - d );
    CrossProduct3d ( &a, &b, &v );
    NormalizeVector3d ( &v );
    phi = sqrt ( e );
  }
  c = cos ( 0.5*phi );
  s = sin ( 0.5*phi );
  SetQuaterniond ( &q, c, s*v.x, s*v.y, s*v.z );
  QuaternionMultd ( &q, qball->q, qball->q );
  IdentTrans3d ( qball->tr );
  QuaternionRotTrans3d ( qball->tr, qball->q );
  xge_xx = x;
  xge_yy = y;
  return true;
} /*xge_QuatRotBalldRotateIt*/

static boolean xge_QuatRotBalldRotateWheel ( xge_widget *er, int key )
{
#define ANG (PI/180.0)
  xge_quatrotballd *qball;
  double           c, s;
  quaterniond      q;

  qball = (xge_quatrotballd*)er->data0;
  c = cos ( 0.5*ANG );
  s = sin ( 0.5*ANG );
  switch ( er->state ) {
case xgestate_NOTHING:
    if ( key & xgemouse_WHEELFW_CHANGE )
      SetQuaterniond ( &q, c, -s, 0.0, 0.0 );
    else
      SetQuaterniond ( &q, c, s, 0.0, 0.0 );
    break;
case xgestate_QUATROT_TURNING1:
    if ( key & xgemouse_WHEELFW_CHANGE )
      SetQuaterniond ( &q, c, 0.0, s, 0.0 );
    else
      SetQuaterniond ( &q, c, 0.0, -s, 0.0 );
    break;
case xgestate_QUATROT_TURNING2:
    if ( key & xgemouse_WHEELFW_CHANGE )
      SetQuaterniond ( &q, c, 0.0, 0.0, s );
    else
      SetQuaterniond ( &q, c, 0.0, 0.0, -s );
    break;
default:
    return false;
  }
  if ( qball->insert )
    QuaternionMultd ( qball->q, &q, qball->q );
  else
    QuaternionMultd ( &q, qball->q, qball->q );
  IdentTrans3d ( qball->tr );
  QuaternionRotTrans3d ( qball->tr, qball->q );
  return true;
#undef ANG
} /*xge_QuatRotBalldRotateWheel*/

static void _xge_QuatRotBalldDrawSpecial ( xge_widget *er, boolean onscreen )
{
  xge_quatrotballd *qball;
  short            xc, yc, R;
  trans3f          tr;

  qball = (xge_quatrotballd*)er->data0;
  xc = qball->xc;
  yc = qball->yc;
  R =  qball->R;
  xgeSetForeground ( xgec_Grey3 );
  xgeDrawArc ( 2*R+1, 2*R+1, xc-R, yc-R, 0, 360*64 );
  xgeSetForeground ( xgec_White );
  Trans3dTof ( qball->tr, &tr );
  _xge_QuatRotBallDrawCircles ( xc, yc, R, &tr );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*_xge_QuatRotBalldDrawSpecial*/

static void _xge_QuatRotBalldOpenSpecial ( xge_widget *er )
{
  xge_quatrotballd *qball;
  short            R;

  qball = (xge_quatrotballd*)er->data0;
  R = qball->R;
  _xge_special_widget->msgproc = xge_EmptyMsg;
  _xge_special_widget->redraw = _xge_QuatRotBalldDrawSpecial;
  _xge_special_widget->id = -1;
  _xge_special_widget->w = _xge_special_widget->h = 2*R+2;
  _xge_special_widget->x = qball->xc-R;
  _xge_special_widget->y = qball->yc-R;
  _xge_special_widget->data0 = qball;
  _xge_special_widget->data1 = er;
  _xge_special_widget->xofs = _xge_special_widget->yofs = 0;
  _xge_special_widget->rpos = 0;
  _xge_special_widget->window_num = er->window_num;
  _xge_special_widget->state = xgestate_NOTHING;
  _xge_special_widget->next = _xge_special_widget->prev =
  _xge_special_widget->up = NULL;
  xge_AddPopup ( _xge_special_widget );
} /*_xge_QuatRotBalldOpenSpecial*/

static void _xge_QuatRotBalldCloseSpecial ( void )
{
  xge_RemovePopup ( false );
  _xge_special_widget->msgproc = xge_EmptyMsg;
  _xge_special_widget->redraw = xge_DrawEmpty;
  _xge_special_widget->id = -1;
  _xge_special_widget->w = _xge_special_widget->h = 0;
  _xge_special_widget->x = _xge_special_widget->y = -1;
  _xge_special_widget->data0 = _xge_special_widget->data1 = NULL;
  _xge_special_widget->xofs = _xge_special_widget->yofs = 0;
  _xge_special_widget->rpos = 0;
  _xge_special_widget->window_num = -1;
  _xge_special_widget->state = xgestate_NOTHING;
  _xge_special_widget->next = _xge_special_widget->prev =
  _xge_special_widget->up = NULL;
} /*_xge_uatRotBalldCloseSpecial*/

boolean xge_QuatRotBalldMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_quatrotballd *qball;
  int              xc, yc, r, d, dx, dy;

  qball = (xge_quatrotballd*)er->data0;
  xc = qball->xc;
  yc = qball->yc;
  r = qball->r2+4;
  switch ( er->state ) {
case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_ENTERING:
      goto setup_cursor;
  case xgemsg_EXITING:
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
      return true;
  case xgemsg_MMOVE:
      d = (x-xc)*(x-xc) + (y-yc)*(y-yc);
      if ( d <= r*r )
        goto setup_cursor;
      else
        xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
      return true;
  case xgemsg_MCLICK:
      d = (x-xc)*(x-xc) + (y-yc)*(y-yc);
      if ( d <= r*r ) {
        if ( key & (xgemouse_WHEELFW_CHANGE | xgemouse_WHEELBK_CHANGE) ) {
          if ( xge_QuatRotBalldRotateWheel ( er, key ) ) {
            xge_callback ( er, xgemsg_QUATROTBALL_COMMAND, 0, x, y );
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
        else if ( (key & xgemouse_LBUTTON_DOWN) &&
                  (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_QUATROT_TURNING1;
          goto enter_turning_state;
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          er->state = xgestate_QUATROT_TURNING2;
enter_turning_state:
          qball->axis_ok = false;
          dx = xc + (int)(((double)(x-xc)*qball->R)/(double)qball->r1 + 0.5);
          dy = yc + (int)(((double)(y-yc)*qball->R)/(double)qball->r1 + 0.5);
          XWarpPointer ( xgedisplay, None, xgewindow, 0, 0, 0, 0, dx, dy );
          xge_xx = dx;
          xge_yy = dy;
          _xge_QuatRotBalldOpenSpecial ( er );
          xge_GrabFocus ( er, true );
          goto redraw_special;
        }
      }
      return true;
  case xgemsg_KEY:
      goto process_key;
  case xgemsg_SPECIAL_KEY:
      goto toggle_insert;
  default:
      return false;
    }

case xgestate_QUATROT_TURNING1:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_turning_state;
      if ( key & xgemouse_RBUTTON_DOWN ) {
        er->state = xgestate_QUATROT_TURNING3;
        qball->axis_ok = false;
      }
      goto rotate_it;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_turning_state;
      else if ( key & (xgemouse_WHEELFW_CHANGE | xgemouse_WHEELBK_CHANGE) ) {
        if ( xge_QuatRotBalldRotateWheel ( er, key ) )
          goto turn_and_redraw_special;
      }
      else if ( key & xgemouse_RBUTTON_DOWN ) {
        er->state = xgestate_QUATROT_TURNING3;
        qball->axis_ok = false;
      }
      return true;
  case xgemsg_KEY:
      goto process_key;
  case xgemsg_SPECIAL_KEY:
      goto toggle_insert;
  default:
      return false;
    }

case xgestate_QUATROT_TURNING3:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        er->state = xgestate_QUATROT_TURNING2;
      else if ( !(key & xgemouse_RBUTTON_DOWN) )
        er->state = xgestate_QUATROT_TURNING1;
      goto rotate_it;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        er->state = xgestate_QUATROT_TURNING2;
      else if ( !(key & xgemouse_RBUTTON_DOWN) )
        er->state = xgestate_QUATROT_TURNING1;
      return true;
  case xgemsg_KEY:
      goto process_key;
  case xgemsg_SPECIAL_KEY:
      goto toggle_insert;
  default:
      return false;
    }

case xgestate_QUATROT_TURNING2:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_turning_state;
      if ( key & xgemouse_LBUTTON_DOWN ) {
        er->state = xgestate_QUATROT_TURNING3;
        qball->axis_ok = false;
      }
rotate_it:
      if ( xge_QuatRotBalldRotateIt ( er, x, y ) )
        goto turn_and_redraw_special;
      return true;
  case xgemsg_MCLICK:
      if ( (key & xgemouse_LBUTTON_DOWN) ) {
        er->state = xgestate_QUATROT_TURNING3;
        qball->axis_ok = false;
      }
      else if ( !(key & xgemouse_RBUTTON_DOWN) ) {
exit_turning_state:
        _xge_QuatRotBalldCloseSpecial ();
        dx = xc + (int)(((double)(x-xc)*qball->r1)/(double)qball->R + 0.5);
        dy = yc + (int)(((double)(y-yc)*qball->r1)/(double)qball->R + 0.5);
        XWarpPointer ( xgedisplay, None, xgewindow, 0, 0, 0, 0, dx, dy );
        xge_xx = dx;
        xge_yy = dy;
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        goto redraw_special;
      }
      else if ( key & (xgemouse_WHEELFW_CHANGE | xgemouse_WHEELBK_CHANGE) ) {
        if ( xge_QuatRotBalldRotateWheel ( er, key ) )
          goto turn_and_redraw_special;
      }
      return true;
  case xgemsg_KEY:
process_key:
      switch ( key ) {
    case 'r':  case 'R':
        SetQuaterniond ( qball->q, 1.0, 0.0, 0.0, 0.0 );
        IdentTrans3d ( qball->tr );
turn_and_redraw_special:
        xge_callback ( er, xgemsg_QUATROTBALL_COMMAND, 0, x, y );
redraw_special:
        xge_Redraw ();
        return true;
    default:
        return false;
      }
  case xgemsg_SPECIAL_KEY:
toggle_insert:
      switch ( key ) {
    case 0xFF63:  /* insert - notebook Toshiba */
    case 0xFF9E:  /* insert - Vobis, Dell */
        qball->insert = !qball->insert;
setup_cursor:
        if ( qball->insert )
          xge_SetCurrentWindowCursor ( xgeCURSOR_ARROW );
        else
          xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
        return true;
    default:
        return false;
      }
  default:
      return false;
    }

default:
    return false;
  }
} /*xge_QuatRotBalldMsg*/

xge_widget *xge_NewQuatRotBalld ( char window_num, xge_widget *prev, int id,
                                  short w, short h, short x, short y, short R,
                                  xge_quatrotballd *qball, quaterniond *q, trans3d *tr,
                                  char *title )
{
  int d;

  qball->q = q;
  qball->tr = tr;
  d = min ( w, h );
  if ( !(d & 0x0001) )
    d --;
  qball->r2 = (short)(d/2);
  qball->r1 = (short)(qball->r2-5);
  qball->R  = R;
  qball->xc = (short)(x + qball->r2);
  qball->yc = (short)(y + qball->r2);
  SetQuaterniond ( qball->q, 1.0, 0.0, 0.0, 0.0 );
  IdentTrans3d ( qball->tr );
  qball->insert = qball->axis_ok = false;
  return (qball->er = xge_NewWidget ( window_num, prev, id, w, h, x, y,
                        qball, title, xge_QuatRotBalldMsg, xge_DrawQuatRotBalld ));
} /*xge_NewQuatRotBalld*/

