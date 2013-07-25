
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

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"

#include "xgedit.h"
#include "xgeprivate.h"


/* ///////////////////////////////////////////////////////////////////////// */
boolean xge_2DwinfMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_2Dwinf *_2Dwin;
  char       c;

  _2Dwin = er->data0;
  switch ( msg ) {
case xgemsg_NULL:
    return true;

case xgemsg_ENTERING:
    if ( _2Dwin->panning )
      xge_SetCurrentWindowCursor ( xgeCURSOR_FLEUR );
    else if ( _2Dwin->selecting_mode )
      xge_SetCurrentWindowCursor ( xgeCURSOR_ARROW );
    else
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
      _2Dwin->inside = true;
    xge_xx = x;
    xge_yy = y;
    if ( _2Dwin->display_coord ) {
      xge_SetClipping ( er );
      er->redraw ( er, true );
    }
    return true;

case xgemsg_EXITING:
    if ( _2Dwin->panning || _2Dwin->selecting_mode )
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    _2Dwin->inside = false;
    if ( _2Dwin->display_coord ) {
      xge_SetClipping ( er );
      er->redraw ( er, true );
    }
    return true;

case xgemsg_MMOVE:
    xge_xx = x;
    xge_yy = y;
    break;

case xgemsg_RESIZE:
    if ( x != er->w || y != er->h ) {
      er->w = x;
      er->h = y;
      xge_2DwinfSetupProjection ( _2Dwin );
      xge_callback ( er, xgemsg_2DWIN_RESIZE, 0, x, y );
      if ( key ) {
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
    }
    return true;

case xgemsg_SPECIAL_KEY:
    return false;

default:
    break;
  }

  switch ( er->state ) {
case xgestate_2DWIN_MOVINGPOINT:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      xge_callback ( er, xgemsg_2DWIN_MOVE_POINT, 0, x, y );
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      break;
  default:
      break;
    }
    break;

case xgestate_2DWIN_PANNING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      xge_2DwinfPan ( er, x, y );
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      break;
  default:
      break;
    }
    break;

case xgestate_2DWIN_ZOOMING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_mode;
      if ( y != _2Dwin->yy ) {
        xge_2DwinfZoom ( er, y );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
exit_mode:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;
  default:
      break;
    }
    break;

case xgestate_2DWIN_SELECTING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      _2Dwin->selection_rect.x1 = x;
      _2Dwin->selection_rect.y1 = y;
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        xge_OrderSelectionRect ( &_2Dwin->selection_rect );
        xge_callback ( er, xgemsg_2DWIN_SELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        _2Dwin->selection_rect.x1 = x;
        _2Dwin->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &_2Dwin->selection_rect );
        xge_callback ( er, xgemsg_2DWIN_SELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_2DWIN_UNSELECTING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      _2Dwin->selection_rect.x1 = x;
      _2Dwin->selection_rect.y1 = y;
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
        xge_OrderSelectionRect ( &_2Dwin->selection_rect );
        xge_callback ( er, xgemsg_2DWIN_UNSELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
        _2Dwin->selection_rect.x1 = x;
        _2Dwin->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &_2Dwin->selection_rect );
        xge_callback ( er, xgemsg_2DWIN_UNSELECT_POINTS, 0, x, y );
exit_select_mode:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;
  default:
      break;
    }
    break;

case xgestate_2DWIN_MOVING_GEOM_WIDGET:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & (xgemouse_LBUTTON_DOWN | xgemouse_RBUTTON_DOWN)) )
        goto exit_widget_mode;
      xge_2DwinfMoveGeomWidget ( _2Dwin, x, y );
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_widget_mode;
      break;
  case xgemsg_KEY:
      switch ( key ) {
    case 'r':  case 'R':
        xge_2DwinfResetGeomWidget ( _2Dwin );
        break;
    default:
        break;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_2DWIN_USING_GEOM_WIDGET:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_widget_mode;
      if ( xge_2DwinfApplyGeomWidget ( _2Dwin, x, y, false ) ) {
        xge_callback ( er, xgemsg_2DWIN_TRANSFORM_POINTS, 0, x, y );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_widget_mode;
      break;
  default:
      break;
    }
    break;

case xgestate_2DWIN_ALTUSING_GEOM_WIDGET:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_widget_mode;
      if ( xge_2DwinfApplyGeomWidget ( _2Dwin, x, y, true ) ) {
        xge_callback ( er, xgemsg_2DWIN_TRANSFORM_POINTS, 0, x, y );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
exit_widget_mode:
        xge_2DwinfExitWidgetMode ( _2Dwin );
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;
  default:
      break;
    }
    break;

case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( _2Dwin->panning ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_2DWIN_PANNING;
          xge_GrabFocus ( er, true );
          _2Dwin->xx = x;
          _2Dwin->yy = y;
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          er->state = xgestate_2DWIN_ZOOMING;
          xge_GrabFocus ( er, true );
          _2Dwin->xx = _2Dwin->yy = y;
          _2Dwin->zoom = (float)(0.5*(_2Dwin->RefBBox.x1-_2Dwin->RefBBox.x0));
        }
      }
      else if ( _2Dwin->selecting_mode ) {
        if (  (key & xgemouse_LBUTTON_DOWN) &&   
              (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_2DWIN_SELECTING;
          goto continue_selecting_mode;
        }
        else if  ( (key & xgemouse_RBUTTON_DOWN) &&   
                   (key & xgemouse_RBUTTON_CHANGE) ) {
          er->state = xgestate_2DWIN_UNSELECTING;
continue_selecting_mode:
          xge_GrabFocus ( er, true );
          _2Dwin->xx = _2Dwin->selection_rect.x0 = _2Dwin->selection_rect.x1 = x;
          _2Dwin->yy = _2Dwin->selection_rect.y0 = _2Dwin->selection_rect.y1 = y;
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
      }
      else if ( _2Dwin->current_tool != xge_2DWIN_NO_TOOL ) {
        c = xge_2DwinfIsItAGeomWidget ( _2Dwin, key, x, y );
        switch ( c ) {
      case 1:
          er->state = xgestate_2DWIN_USING_GEOM_WIDGET;
          goto enter_using;
      case 2:
          er->state = xgestate_2DWIN_ALTUSING_GEOM_WIDGET;
enter_using:
          xge_callback ( er, xgemsg_2DWIN_SAVE_POINTS, 0, x, y );
          _2Dwin->xx = x;
          _2Dwin->yy = y;
          xge_GrabFocus ( er, true );
          break;
      case 3:
          er->state = xgestate_2DWIN_MOVING_GEOM_WIDGET;
          xge_GrabFocus ( er, true );
          _2Dwin->xx = x;
          _2Dwin->yy = y;
          xge_SetClipping ( er );
          er->redraw ( er, true );
          break;
     default:
          break;
        }
      }
      else {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
            /* find out, whether the user has picked a control point */
          if ( xge_callback ( er, xgemsg_2DWIN_PICK_POINT, 0, x, y ) ) {
              /* yes */
            er->state = xgestate_2DWIN_MOVINGPOINT;
            xge_GrabFocus ( er, true );
            _2Dwin->xx = x;
            _2Dwin->yy = y;
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
      }
      break;

  case xgemsg_MMOVE:
      if ( _2Dwin->display_coord ) {
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;

  case xgemsg_KEY:
      switch ( key ) {
    case 'r':  case 'R':  /* reset projection */
        _2Dwin->RefBBox = _2Dwin->DefBBox;
        xge_2DwinfSetupProjection ( _2Dwin );
        xge_callback ( er, xgemsg_2DWIN_PROJCHANGE, key, x, y );
        xge_SetClipping ( er );
        er->redraw ( er, true );
        break;
    case 'f':  case 'F':  /* find bounding box */
        if ( xge_callback ( er, xgemsg_2DWIN_FIND_REFBBOX, key, x, y ) ) {
          xge_2DwinfSetupProjection ( _2Dwin );
          xge_callback ( er, xgemsg_2DWIN_PROJCHANGE, key, x, y );
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
        break;
    case 'u':  case 'U':  /* undo, if possible */
        if ( xge_callback ( er, xgemsg_2DWIN_UNDO, 0, x, y ) ) {
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
        break;
    default:  /* any other key to be consumed by the application */
        return xge_callback ( er, xgemsg_2DWIN_KEY, key, x, y );
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
} /*xge_2DwinfMsg*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_2DwinfSetDefBBox ( xge_2Dwinf *_2Dwin,
                            float x0, float x1, float y0, float y1 )
{
  _2Dwin->DefBBox.x0 = x0;  _2Dwin->DefBBox.x1 = x1;
  _2Dwin->DefBBox.y0 = y0;  _2Dwin->DefBBox.y1 = y1;
} /*xge_2DwinfSetDefBBox*/

void xge_2DwinfSetupProjection ( xge_2Dwinf *_2Dwin )
{
  xge_widget *er;
  CameraRecf *CPos;
  float      w, h, dx, dy;

  er = _2Dwin->er;
  CPos = &_2Dwin->CPos;
  CameraInitFramef ( CPos, true, true,
                     er->w, er->h, er->x, er->y, xge_aspect, 0 );
  dx = (float)er->w;  dy = (float)er->h/xge_aspect;
  w = _2Dwin->RefBBox.x1 - _2Dwin->RefBBox.x0;
  h = _2Dwin->RefBBox.y1 - _2Dwin->RefBBox.y0;
  SetPoint3f ( &CPos->position, (float)(_2Dwin->RefBBox.x0+0.5*w),
               (float)(_2Dwin->RefBBox.y0+0.5*h), 0.0 );
  CPos->psi = CPos->theta = CPos->phi = 0.0;
  if ( dx/w < dy/h ) {
    CPos->vd.para.wdt = w;
    CPos->vd.para.dim_case = 1;
  }
  else {
    CPos->vd.para.hgh = h;
    CPos->vd.para.dim_case = 2;
  }
  CameraSetMappingf ( CPos );
} /*xge_2DwinfSetupProjection*/

void xge_2DwinfPan ( xge_widget *er, short x, short y )
{
  xge_2Dwinf *_2Dwin;
  point2f     a, b, centre;
  vector2f    shift;

  _2Dwin = er->data0;
  if ( x != _2Dwin->xx || y != _2Dwin->yy ) {
    SetPoint2f ( &a, (float)_2Dwin->xx, (float)_2Dwin->yy );
    CameraUnProjectPoint2f ( &_2Dwin->CPos, &a, &a );
    SetPoint2f ( &b, (float)x, (float)y );
    CameraUnProjectPoint2f ( &_2Dwin->CPos, &b, &b );
    SubtractPoints2f ( &a, &b, &shift );
    SetPoint2f ( &centre,
                 (float)(0.5*(_2Dwin->RefBBox.x0+_2Dwin->RefBBox.x1)),
                 (float)(0.5*(_2Dwin->RefBBox.y0+_2Dwin->RefBBox.y1)) );
    AddVector2f ( &centre, &shift, &centre );
/* bound the centre --- to be added */
    SetPoint2f ( &a,
                 (float)(0.5*(_2Dwin->RefBBox.x1-_2Dwin->RefBBox.x0)),
                 (float)(0.5*(_2Dwin->RefBBox.y1-_2Dwin->RefBBox.y0)) );
    _2Dwin->RefBBox.x0 = centre.x-a.x;  _2Dwin->RefBBox.x1 = centre.x+a.x;
    _2Dwin->RefBBox.y0 = centre.y-a.y;  _2Dwin->RefBBox.y1 = centre.y+a.y;
    xge_2DwinfSetupProjection ( _2Dwin );
    xge_callback ( er, xgemsg_2DWIN_PROJCHANGE, 0, x, y );
    _2Dwin->xx = x;
    _2Dwin->yy = y;
  }
} /*xge_2DwinfPan*/

void xge_2DwinfZoom ( xge_widget *er, short y )
{
  xge_2Dwinf *_2Dwin;
  float      f;
  point2f    centre;

  _2Dwin = er->data0;
  f = (float)(y-_2Dwin->xx)/(float)er->h;
  f = (float)(_2Dwin->zoom*exp(f));
  f = (float)max ( f, xge_2DWIN_MIN_ZOOM );
  f = (float)min ( f, xge_2DWIN_MAX_ZOOM );
  SetPoint2f ( &centre,
               (float)(0.5*(_2Dwin->RefBBox.x0+_2Dwin->RefBBox.x1)),
               (float)(0.5*(_2Dwin->RefBBox.y0+_2Dwin->RefBBox.y1)) );
  _2Dwin->RefBBox.x0 = centre.x-f;  _2Dwin->RefBBox.x1 = centre.x+f;
  _2Dwin->RefBBox.y0 = centre.y-f;  _2Dwin->RefBBox.y1 = centre.y+f;
  xge_2DwinfSetupProjection ( _2Dwin );
  xge_callback ( er, xgemsg_2DWIN_PROJCHANGE, 0, 0, y );
  _2Dwin->yy = y;
} /*xge_2DwinfZoom*/

void xge_2DwinfInitProjection ( xge_2Dwinf *_2Dwin,
                                float x0, float x1, float y0, float y1 )
{
  _2Dwin->RefBBox.x0 = x0;  _2Dwin->RefBBox.x1 = x1;
  _2Dwin->RefBBox.y0 = y0;  _2Dwin->RefBBox.y1 = y1;
  CameraSetDepthRangef ( &_2Dwin->CPos, -1.0, 1.0 );
  xge_2DwinfSetupProjection ( _2Dwin );
} /*xge_2DwinfInitProjection*/

void xge_2DwinfResetGeomWidgets ( xge_2Dwinf *_2Dwin )
{
  SetVector2f ( &_2Dwin->scaling_factors, 1.0, 1.0 );
  _2Dwin->scaling_size = 1;
  _2Dwin->tool_mode = 0;
  _2Dwin->rotating_radius = 0;
  _2Dwin->shear_radius = 0.0;
} /*xge_2DwinfResetGeomWidgets*/

void xge_2DwinfResetGeomWidgetPos ( xge_2Dwinf *_2Dwin )
{
  SetPoint2f ( &_2Dwin->scaling_centre, 0.0, 0.0 );
  SetPoint2f ( &_2Dwin->rotating_centre, 0.0, 0.0 );
  SetPoint2f ( &_2Dwin->shear_centre, 0.0, 0.0 );
  SetVector2f ( &_2Dwin->shear_axis[0], 1.0, 0.0 );
  SetVector2f ( &_2Dwin->shear_axis[1], 0.0, 1.0 );
  xge_2DwinfResetGeomWidgets ( _2Dwin );
} /*xge_2DwinfResetGeomWidgetPos*/

void xge_2DwinfEnableGeomWidget ( xge_2Dwinf *_2Dwin, char tool )
{
  _2Dwin->selecting_mode = _2Dwin->panning =
  _2Dwin->moving_tool = _2Dwin->scaling_tool =
  _2Dwin->rotating_tool = _2Dwin->shear_tool = false;
  switch ( tool ) {
case xge_2DWIN_MOVING_TOOL:
    _2Dwin->moving_tool = true;
    goto go_on;
case xge_2DWIN_SCALING_TOOL:
    _2Dwin->scaling_tool = true;
    goto go_on;
case xge_2DWIN_ROTATING_TOOL:
    _2Dwin->rotating_tool = true;
    goto go_on;
case xge_2DWIN_SHEAR_TOOL:
    _2Dwin->shear_tool = true;
    goto go_on;
case xge_2DWIN_SELECTING_TOOL:
    _2Dwin->selecting_mode = true;
    goto go_on;
case xge_2DWIN_PANNING_TOOL:
    _2Dwin->panning = true;
go_on:
    _2Dwin->current_tool = tool;
    break;
default:
    _2Dwin->current_tool = xge_2DWIN_NO_TOOL;
    break;
  }
} /*xge_2DwinfEnableGeomWidget*/

void xge_2DwinfDrawGeomWidgets ( xge_widget *er )
{
  xge_2Dwinf  *_2Dwin;
  CameraRecf  *CPos;
  short       xc, yc, xp, yp, r;
  float       sx, sy;
  point2f     p, q;
  static char dl[2] = {2,3};

  if ( er->state == xgestate_2DWIN_MOVINGPOINT )
    return;
  _2Dwin = er->data0;
  CPos = &_2Dwin->CPos;
  xgeSetForeground ( xgec_Cyan );
  switch ( _2Dwin->current_tool ) {
case xge_2DWIN_MOVING_TOOL:
    xc = (short)(er->x + er->w/2);
    yc = (short)(er->y + er->h/2);
    xgeFillRectangle ( 5, 3, xc-2, yc-1 );
    xgeFillRectangle ( 3, 5, xc-1, yc-2 );
    break;

case xge_2DWIN_SCALING_TOOL:
    if ( (er->state == xgestate_2DWIN_USING_GEOM_WIDGET ||
          er->state == xgestate_2DWIN_ALTUSING_GEOM_WIDGET) ) {
      sx = _2Dwin->scaling_factors.x;
      sy = _2Dwin->scaling_factors.y;
    }
    else
      sx = sy = 1.0;
    CameraProjectPoint2f ( CPos, &_2Dwin->scaling_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    r = (short)(min ( er->w, er->h ));
    _2Dwin->scaling_size = r = (short)((3*r)/10);
    xgeFillRectangle ( 5, 3, xc-2, yc-1 );
    xgeFillRectangle ( 3, 5, xc-1, yc-2 );
    xgeFillRectangle ( 3, 3, xc-1, yc-r*sy-1 );
    xgeFillRectangle ( 3, 3, xc-1, yc+r*sy-1 );
    xgeFillRectangle ( 3, 3, xc-r*sx-1, yc-r*sy-1 );
    xgeFillRectangle ( 3, 3, xc-r*sx-1, yc-1 );
    xgeFillRectangle ( 3, 3, xc-r*sx-1, yc+r*sy-1 );
    xgeFillRectangle ( 3, 3, xc+r*sx-1, yc-r*sy-1 );
    xgeFillRectangle ( 3, 3, xc+r*sx-1, yc-1 );
    xgeFillRectangle ( 3, 3, xc+r*sx-1, yc+r*sy-1 );
    break;

case xge_2DWIN_ROTATING_TOOL:
    CameraProjectPoint2f ( CPos, &_2Dwin->rotating_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    r = (short)(min ( er->w, er->h ));
    _2Dwin->rotating_radius = r = (short)((3*r)/10);
    xgeFillRectangle ( 5, 3, xc-2, yc-1 );
    xgeFillRectangle ( 3, 5, xc-1, yc-2 );
    xgeSetDashes ( 2, dl, 0 );
    xgeSetLineAttributes ( 1, LineOnOffDash, CapRound, JoinRound );
    xgeDrawArc ( 2*r+1, 2*r+1, xc-r, yc-r, 0, 360*64 );
    xgeSetLineAttributes ( 1, LineSolid, CapRound, JoinRound );
    break;

case xge_2DWIN_SHEAR_TOOL:
    xgeSetDashes ( 2, dl, 0 );
    xgeSetLineAttributes ( 1, LineOnOffDash, CapRound, JoinRound );
    CameraProjectPoint2f ( CPos, &_2Dwin->shear_centre, &p );
    xc = (short)(p.x+0.5);   
    yc = (short)(p.y+0.5);
    p.x += (float)(min ( er->w, er->h ));
    CameraUnProjectPoint2f ( CPos, &p, &p );
    _2Dwin->shear_radius = 0.3*(p.x-_2Dwin->shear_centre.x);
    xgeFillRectangle ( 5, 3, xc-2, yc-1 );
    xgeFillRectangle ( 3, 5, xc-1, yc-2 );
    AddVector2Mf ( &_2Dwin->shear_centre, &_2Dwin->shear_axis[1],
                   _2Dwin->shear_radius, &p );
    CameraProjectPoint2f ( CPos, &p, &p );
    xp = (short)(p.x+0.5);
    yp = (short)(p.y+0.5);
    AddVector2Mf ( &_2Dwin->shear_centre, &_2Dwin->shear_axis[1],
                   -_2Dwin->shear_radius, &q );
    CameraProjectPoint2f ( CPos, &q, &q );
    xc = (short)(q.x+0.5);
    yc = (short)(q.y+0.5);
    xgeFillRectangle ( 3, 3, xp-1, yp-1 );
    xgeFillRectangle ( 3, 3, xc-1, yc-1 );
    xgeDrawLine ( xc, yc, xp, yp );
    AddVector2Mf ( &_2Dwin->shear_centre, &_2Dwin->shear_axis[0],
                   _2Dwin->shear_radius, &p );
    CameraProjectPoint2f ( CPos, &p, &p );
    xp = (short)(p.x+0.5);
    yp = (short)(p.y+0.5);
    AddVector2Mf ( &_2Dwin->shear_centre, &_2Dwin->shear_axis[0],
                   -_2Dwin->shear_radius, &q );
    CameraProjectPoint2f ( CPos, &q, &q );
    xc = (short)(q.x+0.5);
    yc = (short)(q.y+0.5);
    xgeFillRectangle ( 3, 3, xp-1, yp-1 );
    xgeFillRectangle ( 3, 3, xc-1, yc-1 );
    xgeDrawLine ( xc, yc, xp, yp );
    xgeSetLineAttributes ( 1, LineSolid, CapRound, JoinRound );
    break;

default:
    return;
  }
} /*xge_2DwinfDrawGeomWidgets*/

char xge_2DwinfIsItAGeomWidget ( xge_2Dwinf *_2Dwin, int key, short x, short y )
{
/* return value 0 - no widget has been pointed */
/* 1 - the widget should be applied */
/* 2 - the widget should be applied in the alternative mode */
/* 3 - the widget should be moved */
  xge_widget *er;
  short      xc, yc, r;
  CameraRecf *CPos;
  point2f    p;
  char       result;

  result = 0;
  if ( !(((key & xgemouse_LBUTTON_DOWN) && (key & xgemouse_LBUTTON_CHANGE)) ||
         ((key & xgemouse_RBUTTON_DOWN) && (key & xgemouse_RBUTTON_CHANGE))) )
    return result;
  er = _2Dwin->er;
  CPos = &_2Dwin->CPos;
  switch ( _2Dwin->current_tool ) {
case xge_2DWIN_MOVING_TOOL:
    xc = (short)(er->x + er->w/2);
    yc = (short)(er->y + er->h/2);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST )
      result = 1;
    break;

case xge_2DWIN_SCALING_TOOL:
    CameraProjectPoint2f ( CPos, &_2Dwin->scaling_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _2Dwin->saved_centre = _2Dwin->scaling_centre;
      result = 3;
    }
    else {
      if ( (key & xgemouse_LBUTTON_DOWN) &&
           (key & xgemouse_LBUTTON_CHANGE) )
        result = 1;
      else
        result = 2;
      r = _2Dwin->scaling_size;
      if ( abs(x-(xc-r))+abs(y-(yc-r)) < xge_MINDIST )
        _2Dwin->tool_mode = 1;
      else if ( abs(x-xc)+abs(y-(yc-r)) < xge_MINDIST )
        _2Dwin->tool_mode = 2;
      else if ( abs(x-(xc+r))+abs(y-(yc-r)) < xge_MINDIST )
        _2Dwin->tool_mode = 3;
      else if ( abs(x-(xc-r))+abs(y-yc) < xge_MINDIST )
        _2Dwin->tool_mode = 4;
      else if ( abs(x-(xc+r))+abs(y-yc) < xge_MINDIST )
        _2Dwin->tool_mode = 5;
      else if ( abs(x-(xc-r))+abs(y-(yc+r)) < xge_MINDIST )
        _2Dwin->tool_mode = 6;
      else if ( abs(x-xc)+abs(y-(yc+r)) < xge_MINDIST )
        _2Dwin->tool_mode = 7;
      else if ( abs(x-(xc+r))+abs(y-(yc+r)) < xge_MINDIST )
        _2Dwin->tool_mode = 8;
      else result = 0;
    }
    break;

case xge_2DWIN_ROTATING_TOOL:
    CameraProjectPoint2f ( CPos, &_2Dwin->rotating_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _2Dwin->saved_centre = _2Dwin->rotating_centre;
      result = 3;
    }
    else {
      r = (short)sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) );
      if ( abs ( r-_2Dwin->rotating_radius ) < xge_MINDIST )
        result = 1;
    }
    break;

case xge_2DWIN_SHEAR_TOOL:
    CameraProjectPoint2f ( CPos, &_2Dwin->shear_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _2Dwin->saved_centre = _2Dwin->shear_centre;
      _2Dwin->tool_mode = 0;
      result = 3;
      break;
    }
    if ( (key & xgemouse_LBUTTON_DOWN) || (key & xgemouse_LBUTTON_CHANGE) )
      result = 1;
    else if ( (key & xgemouse_RBUTTON_DOWN) || (key & xgemouse_RBUTTON_CHANGE) )
      result = 3;
    else {
      result = 0;
      break;
    }
    AddVector2Mf ( &_2Dwin->shear_centre, &_2Dwin->shear_axis[0],
                  _2Dwin->shear_radius, &p );
    CameraProjectPoint2f ( CPos, &p, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _2Dwin->saved_centre = _2Dwin->shear_axis[0];
      _2Dwin->tool_mode = 1;
      break;
    }
    AddVector2Mf ( &_2Dwin->shear_centre, &_2Dwin->shear_axis[0],
                  -_2Dwin->shear_radius, &p );
    CameraProjectPoint2f ( CPos, &p, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _2Dwin->saved_centre = _2Dwin->shear_axis[0];
      _2Dwin->tool_mode = 2;
      break;
    }
    AddVector2Mf ( &_2Dwin->shear_centre, &_2Dwin->shear_axis[1],
                  _2Dwin->shear_radius, &p );
    CameraProjectPoint2f ( CPos, &p, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _2Dwin->saved_centre = _2Dwin->shear_axis[1];
      _2Dwin->tool_mode = 3;
      break;
    }
    AddVector2Mf ( &_2Dwin->shear_centre, &_2Dwin->shear_axis[1],
                  -_2Dwin->shear_radius, &p );
    CameraProjectPoint2f ( CPos, &p, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _2Dwin->saved_centre = _2Dwin->shear_axis[1];
      _2Dwin->tool_mode = 4;
      break;
    }
    result = 0;
    break;

default:
    break;
  }
  return result;
} /*xge_2DwinfIsItAGeomWidget*/

void xge_2DwinfMoveGeomWidget ( xge_2Dwinf *_2Dwin, short x, short y )
{
  CameraRecf *CPos;
  point2f    p, q, r;
  vector2f   v;

  CPos = &_2Dwin->CPos;
  SetPoint2f ( &p, (float)_2Dwin->xx, (float)_2Dwin->yy );
  SetPoint2f ( &q, (float)x, (float)y );
  CameraUnProjectPoint2f ( CPos, &q, &r );
  CameraUnProjectPoint2f ( CPos, &p, &q );
  SubtractPoints2f ( &r, &q, &v );
  switch ( _2Dwin->current_tool ) {
case xge_2DWIN_SCALING_TOOL:
    AddVector2f ( &_2Dwin->saved_centre, &v, &_2Dwin->scaling_centre );
    break;
case xge_2DWIN_ROTATING_TOOL:
    AddVector2f ( &_2Dwin->saved_centre, &v, &_2Dwin->rotating_centre );
    break;
case xge_2DWIN_SHEAR_TOOL:
    switch ( _2Dwin->tool_mode ) {
  case 0:  /* move the centre */
      AddVector2f ( &_2Dwin->saved_centre, &v, &_2Dwin->shear_centre );
      break;
  case 1:  /* move the positive halfaxis x */
      AddVector2Mf ( &_2Dwin->saved_centre, &v,
                     (float)(1.0/_2Dwin->shear_radius), &_2Dwin->shear_axis[0] );
      goto cont_shear1;
  case 2:  /* move the negative halfaxis x */
      AddVector2Mf ( &_2Dwin->saved_centre, &v,
                     (float)(-1.0/_2Dwin->shear_radius), &_2Dwin->shear_axis[0] );
cont_shear1:
      SetVector2f ( &_2Dwin->shear_axis[1],
                    -_2Dwin->shear_axis[0].y, _2Dwin->shear_axis[0].x );
      NormalizeVector2f ( &_2Dwin->shear_axis[1] );
      break;
  case 3:
      AddVector2Mf ( &_2Dwin->saved_centre, &v,
                     (float)(1.0/_2Dwin->shear_radius), &_2Dwin->shear_axis[1] );
      goto cont_shear2;
  case 4:
      AddVector2Mf ( &_2Dwin->saved_centre, &v,
                     (float)(-1.0/_2Dwin->shear_radius), &_2Dwin->shear_axis[1] );
cont_shear2:
      SetVector2f ( &_2Dwin->shear_axis[0],
                    _2Dwin->shear_axis[1].y, -_2Dwin->shear_axis[1].x );
      NormalizeVector2f ( &_2Dwin->shear_axis[0] );
      break;
  default:
      break;
    }
    break;
default:
    break;
  }
  xge_callback ( _2Dwin->er,
                 xgemsg_2DWIN_CHANGE_TRANS, _2Dwin->current_tool, 0, 0 );
} /*xge_2DwinfMoveGeomWidget*/

boolean xge_2DwinfApplyGeomWidget ( xge_2Dwinf *_2Dwin, short x, short y,
                                    boolean alt )
{
  CameraRecf *CPos;
  point2f    p, q, r;
  vector2f   v, *sf;
  float      alpha, beta, x0, y0, x1, y1;
  int        xx, yy, xc, yc, s;
  trans2f    tr1, tr2;

  CPos = &_2Dwin->CPos;
  IdentTrans2f ( &_2Dwin->gwtrans );
  switch ( _2Dwin->current_tool ) {

case xge_2DWIN_MOVING_TOOL:
    SetPoint2f ( &p, (float)_2Dwin->xx, (float)_2Dwin->yy );
    SetPoint2f ( &q, (float)x, (float)y );
    CameraUnProjectPoint2f ( CPos, &q, &r );
    CameraUnProjectPoint2f ( CPos, &p, &q );
    SubtractPoints2f ( &r, &q, &v );
    ShiftTrans2f ( &_2Dwin->gwtrans, v.x, v.y );
    _2Dwin->trans_params = v;
    break;

case xge_2DWIN_SCALING_TOOL:
    sf = &_2Dwin->scaling_factors;
    CameraProjectPoint2f ( CPos, &_2Dwin->scaling_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    s = _2Dwin->scaling_size;
    switch ( _2Dwin->tool_mode ) {
  case 1: xx = (short)(xc - s);  yy = (short)(yc - s);  break;
  case 2: xx = (short)xc;        yy = (short)(yc - s);  break;
  case 3: xx = (short)(xc + s);  yy = (short)(yc - s);  break;
  case 4: xx = (short)(xc - s);  yy = (short)yc;        break;
  case 5: xx = (short)(xc + s);  yy = (short)yc;        break;
  case 6: xx = (short)(xc - s);  yy = (short)(yc + s);  break;
  case 7: xx = (short)xc;        yy = (short)(yc + s);  break;
  case 8: xx = (short)(xc + s);  yy = (short)(yc + s);  break;
 default: return false;
    }
    if ( xx != xc )
      sf->x = (float)(x-xc)/(float)(xx-xc);
    else
      sf->x = 1.0;
    if ( yy != yc )
      sf->y = (float)(y-yc)/(float)(yy-yc);
    else
      sf->y = 1.0;
    if ( alt )
      sf->x = sf->y = (float)(0.5*(sf->x+sf->y));
    ShiftTrans2f ( &_2Dwin->gwtrans,
                   -_2Dwin->scaling_centre.x, -_2Dwin->scaling_centre.y );
    ScaleTrans2f ( &_2Dwin->gwtrans, sf->x, sf->y );
    ShiftTrans2f ( &_2Dwin->gwtrans,
                   _2Dwin->scaling_centre.x, _2Dwin->scaling_centre.y );
    _2Dwin->trans_params = *sf;
    break;

case xge_2DWIN_ROTATING_TOOL:
    CameraProjectPoint2f ( CPos, &_2Dwin->rotating_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( ((_2Dwin->xx != xc) || (_2Dwin->yy != yc)) && ((x != xc) || (y != yc)) ) {
      beta = (float)atan2 ( _2Dwin->yy-yc, _2Dwin->xx-xc );
      alpha = (float)(beta - atan2( y-yc, x-xc ));
      ShiftTrans2f ( &_2Dwin->gwtrans,
                     -_2Dwin->rotating_centre.x, -_2Dwin->rotating_centre.y );
      RotTrans2f ( &_2Dwin->gwtrans, alpha );
      ShiftTrans2f ( &_2Dwin->gwtrans,
                     _2Dwin->rotating_centre.x, _2Dwin->rotating_centre.y );
      _2Dwin->trans_params.x = alpha;
    }
    else
      return false;
    break;

case xge_2DWIN_SHEAR_TOOL:
    SetPoint2f ( &p, (double)_2Dwin->xx, (double)_2Dwin->yy );
    SetPoint2f ( &q, (double)x, (double)y );
    CameraUnProjectPoint2f ( CPos, &q, &r );
    CameraUnProjectPoint2f ( CPos, &p, &q );
    SubtractPoints2f ( &r, &q, &v );
    switch ( _2Dwin->tool_mode ) {
  case 1:
      AddVector2Mf ( &_2Dwin->saved_centre, &v,
                     (float)(1.0/_2Dwin->shear_radius), &_2Dwin->shear_axis[0]);
      goto apply1;
  case 2:
      AddVector2Mf ( &_2Dwin->saved_centre, &v,
                     (float)(-1.0/_2Dwin->shear_radius), &_2Dwin->shear_axis[0]);
apply1:
      x0 = _2Dwin->shear_axis[0].x;  y0 = _2Dwin->shear_axis[0].y;
      x1 = _2Dwin->shear_axis[1].x;  y1 = _2Dwin->shear_axis[1].y;
      tr1.U0.a11 = x0*y1+x1*x1;  tr1.U0.a12 = (y1-x0)*x1;
      tr1.U0.a21 = (y0+x1)*y1;   tr1.U0.a22 = y1*y1-x1*y0;
      break;
  case 3:
      AddVector2Mf ( &_2Dwin->saved_centre, &v,
                     (float)(1.0/_2Dwin->shear_radius), &_2Dwin->shear_axis[1]);
      goto apply2;
  case 4:
      AddVector2Mf ( &_2Dwin->saved_centre, &v,
                     (float)(-1.0/_2Dwin->shear_radius), &_2Dwin->shear_axis[1]);
apply2:
      x0 = _2Dwin->shear_axis[0].x;  y0 = _2Dwin->shear_axis[0].y;
      x1 = _2Dwin->shear_axis[1].x;  y1 = _2Dwin->shear_axis[1].y;
      tr1.U0.a11 = x0*x0-x1*y0;  tr1.U0.a12 = x0*(x1+y0);
      tr1.U0.a21 = (x0-y1)*y0;   tr1.U0.a22 = y0*y0+x0*y1;
      break;
  default:
      return false;
    }
    tr1.U0.a13 = 0.0;  tr1.U0.a23 = 0.0;
    x0 = (tr1.U0.a11*tr1.U0.a22 - tr1.U0.a12*tr1.U0.a21);
    if ( x0 > 0.0 )
      tr1.U1.detsgn = 1;
    else if ( x0 < 0.0 )
      tr1.U1.detsgn = -1;
    else
      tr1.U1.detsgn = 0;
    ShiftTrans2f ( &_2Dwin->gwtrans,
                   -_2Dwin->shear_centre.x, -_2Dwin->shear_centre.y );
    CompTrans2f ( &tr2, &tr1, &_2Dwin->gwtrans );
    ShiftTrans2f ( &tr2, _2Dwin->shear_centre.x, _2Dwin->shear_centre.y );
    _2Dwin->gwtrans = tr2;
    break;

default:
    return false;
  }
  return true;
} /*xge_2DwinfApplyGeomWidget*/

void xge_2DwinfExitWidgetMode ( xge_2Dwinf *_2Dwin )
{
  if ( _2Dwin->current_tool == xge_2DWIN_SHEAR_TOOL )
    switch ( _2Dwin->tool_mode ) {
  case 1:  case 2:
      NormalizeVector2f ( &_2Dwin->shear_axis[1] );
      SetVector2f ( &_2Dwin->shear_axis[0],
                    _2Dwin->shear_axis[1].y, -_2Dwin->shear_axis[1].x );
      break;
  case 3:  case 4:
      NormalizeVector2f ( &_2Dwin->shear_axis[0] );
      SetVector2f ( &_2Dwin->shear_axis[1],
                    -_2Dwin->shear_axis[0].y, _2Dwin->shear_axis[0].x );
      break;
  default:
      break;
    }
} /*xge_2DwindExitWidgetMode*/

void xge_2DwinfResetGeomWidget ( xge_2Dwinf *_2Dwin )
{
  switch ( _2Dwin->current_tool ) {
case xge_2DWIN_SCALING_TOOL:
    SetPoint2f ( &_2Dwin->scaling_centre, 0.0, 0.0 );
    break;
case xge_2DWIN_ROTATING_TOOL:
    SetPoint2f ( &_2Dwin->rotating_centre, 0.0, 0.0 );
    break;
case xge_2DWIN_SHEAR_TOOL:
    if ( _2Dwin->tool_mode == 0 )
      SetPoint2f ( &_2Dwin->shear_centre, 0.0, 0.0 );
    else {
      SetVector2f ( &_2Dwin->shear_axis[0], 1.0, 0.0 );
      SetVector2f ( &_2Dwin->shear_axis[1], 0.0, 1.0 );
    }
    break;
default:
    return;
  }
  _2Dwin->er->state = xgestate_NOTHING;
  xge_ReleaseFocus ( _2Dwin->er );
  xge_SetClipping ( _2Dwin->er );
  _2Dwin->er->redraw ( _2Dwin->er, true );
} /*xge_2DwinfResetGeomWidget*/

void xge_2DwinfDrawCursorPos ( xge_2Dwinf *_2Dwin, short x, short y )
{
  xge_widget *er;
  char       s[10];
  point2f    p, q;
  int        xx;

  if ( _2Dwin->display_coord ) {
    er = _2Dwin->er;
    SetPoint2f ( &p, x, y );
    CameraUnProjectPoint2f ( &_2Dwin->CPos, &p, &q );
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
} /*xge_2DwinfDrawCursorPos*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *xge_New2Dwinf ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            xge_2Dwinf *_2Dwin,
                            void (*redraw)(xge_widget*, boolean) )
{
  xge_widget *er;

  er = xge_NewWidget ( window_num, prev, id, w, h, x, y,
                       _2Dwin, &_2Dwin->CPos, xge_2DwinfMsg, redraw );
  if ( er ) {
    _2Dwin->er = er;
    _2Dwin->panning = _2Dwin->selecting_mode = false;
    _2Dwin->display_coord = _2Dwin->inside = false;
    _2Dwin->current_tool = xge_2DWIN_NO_TOOL;
    _2Dwin->moving_tool = _2Dwin->scaling_tool = _2Dwin->rotating_tool =
    _2Dwin->shear_tool = false;
    xge_2DwinfSetDefBBox ( _2Dwin, -1.0, 1.0, -1.0, 1.0 );
    xge_2DwinfInitProjection ( _2Dwin, -1.0, 1.0, -1.0, 1.0 );
    xge_2DwinfResetGeomWidgetPos ( _2Dwin );
  }
  return er;
} /*xge_New2Dwinf*/

