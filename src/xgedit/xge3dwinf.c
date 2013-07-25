
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

/*#define SHOW_SIZE */

/* ///////////////////////////////////////////////////////////////////////// */
/* the four 3D window widget is a four window widget with three windows      */
/* showing parallel projections of a 3D object on three planes and one       */
/* perspective projection; it is intended to allow a~user to manipulate with */
/* this object via its control points in a unified way for all applications  */
/* ///////////////////////////////////////////////////////////////////////// */
/* the message processing procedures are part of the xgedit library, while   */
/* the drawing procedures must be supplied by the application.               */

/* ///////////////////////////////////////////////////////////////////////// */
static boolean xge_3DwinfParaMsg ( xge_widget *er,
                                   int msg, int key, short x, short y )
{
  int        id, i;
  xge_3Dwinf *_3Dwin;
  char       c;

  id = er->id % 0x03;
  _3Dwin = er->data0;
  switch ( msg ) {
case xgemsg_NULL:
    return true;

case xgemsg_ENTERING:
    if ( _3Dwin->panning )
      xge_SetCurrentWindowCursor ( xgeCURSOR_FLEUR );
    else if ( _3Dwin->selecting_mode || _3Dwin->special_selecting_mode )
      xge_SetCurrentWindowCursor ( xgeCURSOR_ARROW );
    else
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    _3Dwin->CoordWin = (char)id;
    if ( _3Dwin->display_coord )
      goto update_3pictures;
    return true;

case xgemsg_EXITING:
    if ( _3Dwin->panning || _3Dwin->selecting_mode ||
         _3Dwin->special_selecting_mode )
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    _3Dwin->CoordWin = -1;
    if ( _3Dwin->display_coord )
      goto update_3pictures;
    return true;

case xgemsg_MMOVE:
    xge_xx = x;
    xge_yy = y;
    break;

case xgemsg_RESIZE:
    if ( id == 0 ) {
      xge_3DwinfUpdatePerspProj ( _3Dwin );
      xge_3DwinfSetupParProj ( _3Dwin, &_3Dwin->RefBBox );
      xge_callback ( er, xgemsg_3DWIN_RESIZE, 0, x, y );
      if ( key )
        goto update_4pictures;
    }
    return true;

case xgemsg_SPECIAL_KEY:
    return false;

default:
    break;
  }

  switch ( er->state ) {
case xgestate_3DWIN_MOVINGPOINT:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      xge_callback ( er, xgemsg_3DWIN_MOVE_POINT, 0, x, y );
      goto update_4pictures;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_PARPANNING:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      break;
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      xge_3DwinfPanParWindows ( er, x, y );
      goto update_3pictures;
  default:
      break;
    }
    break;

case xgestate_3DWIN_PARZOOMING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( key & xgemouse_RBUTTON_DOWN ) {
        if ( y != _3Dwin->yy ) {
          xge_3DwinfZoomParWindows ( er, y );
          _3Dwin->yy = y;
          goto update_3pictures;
        }
      }
      else
        goto exit_mode;
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN ) )
        goto exit_mode;
      break;
  default:
      break;
    }
    break;

case xgestate_RESIZING_X:
case xgestate_RESIZING_Y:
case xgestate_RESIZING_XY:
    break;

case xgestate_3DWIN_SELECTING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      _3Dwin->selection_rect.x1 = x;
      _3Dwin->selection_rect.y1 = y;
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        xge_OrderSelectionRect ( &_3Dwin->selection_rect );
        xge_callback ( er, xgemsg_3DWIN_SELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        _3Dwin->selection_rect.x1 = x;
        _3Dwin->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &_3Dwin->selection_rect );
        xge_callback ( er, xgemsg_3DWIN_SELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_UNSELECTING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      _3Dwin->selection_rect.x1 = x;
      _3Dwin->selection_rect.y1 = y;
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
        xge_OrderSelectionRect ( &_3Dwin->selection_rect );
        xge_callback ( er, xgemsg_3DWIN_UNSELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
        _3Dwin->selection_rect.x1 = x;
        _3Dwin->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &_3Dwin->selection_rect );
        xge_callback ( er, xgemsg_3DWIN_UNSELECT_POINTS, 0, x, y );
exit_select_mode:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( _3Dwin->fww.er );
        xge_SetClipping ( _3Dwin->fww.er );
        _3Dwin->fww.er->redraw ( _3Dwin->fww.er, true );
      }
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_MOVING_GEOM_WIDGET:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( !(key & (xgemouse_LBUTTON_DOWN | xgemouse_RBUTTON_DOWN)) )
        goto exit_widget_mode;
      break;
  case xgemsg_MMOVE:
      if ( !(key & (xgemouse_LBUTTON_DOWN | xgemouse_RBUTTON_DOWN)) )
        goto exit_widget_mode;
      xge_3DwinfMoveGeomWidget ( _3Dwin, id, x, y );
      goto update_3pictures;
  case xgemsg_KEY:
      switch ( key ) {
    case 'r':  case 'R':
        xge_3DwinfResetGeomWidget ( _3Dwin );
        break;
    default:
        break;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_USING_GEOM_WIDGET:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_widget_mode;
      break;
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_widget_mode;
      if ( xge_3DwinfApplyGeomWidget ( _3Dwin, id, x, y, false ) ) {
        xge_callback ( _3Dwin->cwin[id], xgemsg_3DWIN_TRANSFORM_POINTS, 0, x, y );
        goto update_4pictures;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_ALTUSING_GEOM_WIDGET:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_widget_mode;
      break;
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
exit_widget_mode:
        xge_3DwinfExitWidgetMode ( _3Dwin );
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( _3Dwin->fww.er );
        goto update_3pictures;
      }
      if ( xge_3DwinfApplyGeomWidget ( _3Dwin, id, x, y, true ) ) {
        xge_callback ( _3Dwin->cwin[id], xgemsg_3DWIN_TRANSFORM_POINTS, 0, x, y );
        goto update_4pictures;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_USING_SPECIAL_WIDGET:
    break;

case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( _3Dwin->display_coord )
        goto update_3pictures;
      break;
  case xgemsg_MCLICK:
      if ( _3Dwin->panning ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_3DWIN_PARPANNING;
          xge_GrabFocus ( _3Dwin->fww.er, true );
          xge_GrabFocus ( er, true );
          _3Dwin->xx = x;
          _3Dwin->yy = y;
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          er->state = xgestate_3DWIN_PARZOOMING;
          xge_GrabFocus ( _3Dwin->fww.er, true );
          xge_GrabFocus ( er, true );
          _3Dwin->xx = _3Dwin->yy = y;
        }
      }
      else if ( _3Dwin->selecting_mode ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_3DWIN_SELECTING;
          goto continue_selecting_mode;
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) && 
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          er->state = xgestate_3DWIN_UNSELECTING;
continue_selecting_mode:
          xge_GrabFocus ( _3Dwin->fww.er, true );
          xge_GrabFocus ( er, true );
          _3Dwin->xx = _3Dwin->selection_rect.x0 = _3Dwin->selection_rect.x1 = x;
          _3Dwin->yy = _3Dwin->selection_rect.y0 = _3Dwin->selection_rect.y1 = y;
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
      }
      else if ( _3Dwin->special_selecting_mode ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          if ( xge_callback ( er, xgemsg_3DWIN_SPECIAL_SELECT, id, x, y ) ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          if ( xge_callback ( er, xgemsg_3DWIN_SPECIAL_UNSELECT, id, x, y ) ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
      }
      else if ( _3Dwin->current_tool != xge_3DWIN_NO_TOOL ) {
        c = xge_3DwinfIsItAGeomWidget ( _3Dwin, id, key, x, y );
        switch ( c ) {
      case 1:
      case 2:
          if ( (key & xgemouse_LBUTTON_DOWN) &&
               (key & xgemouse_LBUTTON_CHANGE) )
            er->state = xgestate_3DWIN_USING_GEOM_WIDGET;
          else if ( (key & xgemouse_RBUTTON_DOWN) &&
                    (key & xgemouse_RBUTTON_CHANGE) )
            er->state = xgestate_3DWIN_ALTUSING_GEOM_WIDGET;
          else
            break;
          xge_callback ( er, xgemsg_3DWIN_SAVE_POINTS, 0, x, y );
          _3Dwin->xx = x;
          _3Dwin->yy = y;
          xge_GrabFocus ( _3Dwin->fww.er, true );
          xge_GrabFocus ( er, true );
          break;
      case 3:
          if ( ((key & xgemouse_LBUTTON_DOWN) &&
                (key & xgemouse_LBUTTON_CHANGE)) ||
               ((key & xgemouse_RBUTTON_DOWN) &&
                (key & xgemouse_RBUTTON_CHANGE)) ) {
            er->state = xgestate_3DWIN_MOVING_GEOM_WIDGET;
            xge_GrabFocus ( _3Dwin->fww.er, true );
            xge_GrabFocus ( er, true );
            _3Dwin->xx = x;
            _3Dwin->yy = y;
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
          break;
      default:
          break;
        }
      }
      else {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
            /* find out whether the user has picked a control point */
          if ( xge_callback ( er, xgemsg_3DWIN_PICK_POINT, 0, x, y ) ) {
              /* yes */
            er->state = xgestate_3DWIN_MOVINGPOINT;
            xge_GrabFocus ( _3Dwin->fww.er, true );
            xge_GrabFocus ( er, true );
            _3Dwin->xx = x;
            _3Dwin->yy = y;
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
      }
      break;
  case xgemsg_KEY:
      switch ( key ) {
    case 'r':  case 'R':  /* reset projection */
        _3Dwin->RefBBox = _3Dwin->DefBBox;
        xge_3DwinfSetupParProj ( _3Dwin, &_3Dwin->RefBBox );
        xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, key, x, y );
        goto update_3pictures;
    case 'f':  case 'F':  /* find bounding box */
        if ( xge_callback ( er, xgemsg_3DWIN_FIND_REFBBOX, key, x, y ) ) {
          xge_3DwinfSetupParProj ( _3Dwin, &_3Dwin->RefBBox );
          xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, key, x, y );
          goto update_3pictures;
        }
        else
          break;
    case 'u':  case 'U':  /* undo, if possible */
        if ( xge_callback ( er, xgemsg_3DWIN_UNDO, id, x, y ) ) {
          xge_SetClipping ( _3Dwin->fww.er );
          _3Dwin->fww.er->redraw ( _3Dwin->fww.er, true );
        }
        break;
    default:  /* any other key to be consumed by the application */
        return xge_callback ( er, xgemsg_3DWIN_KEY, key, x, y );
      }
      break;

  default:
      break;
    }
    break;

default:
    switch ( msg ) {
  case xgemsg_MCLICK:
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
exit_mode:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( _3Dwin->fww.er );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;

  default:
      break;
    }
    break;
  }
  return true;

update_3pictures:
  i = _3Dwin->fww.zoomwin;
  if ( i >= 0 && i < 3 ) {
    xge_SetClipping ( _3Dwin->cwin[i] );
    _3Dwin->cwin[i]->redraw ( _3Dwin->cwin[i], true );
  }
  else {
    for ( i = 0; i < 3; i++ ) {
      xge_SetClipping ( _3Dwin->cwin[i] );
      _3Dwin->cwin[i]->redraw ( _3Dwin->cwin[i], true );
    }
  }
  return true;

update_4pictures:
  i = _3Dwin->fww.zoomwin;
  if ( i >= 0 && i < 4 ) {
    xge_SetClipping ( _3Dwin->cwin[i] );
    _3Dwin->cwin[i]->redraw ( _3Dwin->cwin[i], true );
  }
  else {
    for ( i = 0; i < 4; i++ ) {
      xge_SetClipping ( _3Dwin->cwin[i] );
      _3Dwin->cwin[i]->redraw ( _3Dwin->cwin[i], true );
    }
  }
  return true;
} /*xge_3DwinfParaMsg*/

void xge_3DwinfPanParWindows ( xge_widget *er, short x, short y )
{
  xge_3Dwinf *_3Dwin;
  int        id;
  point3f    a, b;
  float      w, h, d;
  vector3f   shift;

  _3Dwin = er->data0;
  if ( x != _3Dwin->xx || y != _3Dwin->yy ) {
    id = er->id & 0x03;
    SetPoint3f ( &a, (float)_3Dwin->xx, (float)_3Dwin->yy, 0.0 );
    CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &a, &a );
    SetPoint3f ( &b, (float)x, (float)y, 0.0 );
    CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &b, &b );
    SubtractPoints3f ( &a, &b, &shift );
/* bound the centre --- to be added */
    w = (float)(0.5*(_3Dwin->RefBBox.x1-_3Dwin->RefBBox.x0));
    h = (float)(0.5*(_3Dwin->RefBBox.y1-_3Dwin->RefBBox.y0));
    d = (float)(0.5*(_3Dwin->RefBBox.z1-_3Dwin->RefBBox.z0));
    SetPoint3f ( &a, (float)(0.5*(_3Dwin->RefBBox.x0+_3Dwin->RefBBox.x1)),
                     (float)(0.5*(_3Dwin->RefBBox.y0+_3Dwin->RefBBox.y1)),
                     (float)(0.5*(_3Dwin->RefBBox.z0+_3Dwin->RefBBox.z1)) );
    AddVector3f ( &a, &shift, &a );
    _3Dwin->RefBBox.x0 = a.x-w;  _3Dwin->RefBBox.x1 = a.x+w;
    _3Dwin->RefBBox.y0 = a.y-h;  _3Dwin->RefBBox.y1 = a.y+h;
    _3Dwin->RefBBox.z0 = a.z-d;  _3Dwin->RefBBox.z1 = a.z+d;
    xge_3DwinfSetupParProj ( _3Dwin, &_3Dwin->RefBBox );
    xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, id, x, y );
    _3Dwin->xx = x;
    _3Dwin->yy = y;
  }
} /*xge_3DwinfPanParWindows*/

void xge_3DwinfZoomParWindows ( xge_widget *er, short y )
{
  xge_3Dwinf *_3Dwin;
  float      f, w, h, d;
  point3f    centre;

  _3Dwin = er->data0;
  f = (float)(y-_3Dwin->xx)/(float)er->h;
  f = (float)exp(f);
  _3Dwin->xx = y;
  SetPoint3f ( &centre,
               (float)(0.5*(_3Dwin->RefBBox.x0+_3Dwin->RefBBox.x1)),
               (float)(0.5*(_3Dwin->RefBBox.y0+_3Dwin->RefBBox.y1)),
               (float)(0.5*(_3Dwin->RefBBox.z0+_3Dwin->RefBBox.z1)) );
  w = (float)(0.5*(_3Dwin->RefBBox.x1-_3Dwin->RefBBox.x0));
  h = (float)(0.5*(_3Dwin->RefBBox.y1-_3Dwin->RefBBox.y0));
  d = (float)(0.5*(_3Dwin->RefBBox.z1-_3Dwin->RefBBox.z0));
  _3Dwin->RefBBox.x0 = centre.x-f*w;  _3Dwin->RefBBox.x1 = centre.x+f*w;
  _3Dwin->RefBBox.y0 = centre.y-f*h;  _3Dwin->RefBBox.y1 = centre.y+f*h;
  _3Dwin->RefBBox.z0 = centre.z-f*d;  _3Dwin->RefBBox.z1 = centre.z+f*d;
  xge_3DwinfSetupParProj ( _3Dwin, &_3Dwin->RefBBox );
  xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, er->id & 0x03, 0, y );
  _3Dwin->yy = y;
} /*xge_3DwinfZoomParWindows*/

void xge_3DwinfPanPerspWindow ( xge_widget *er, short x, short y )
{
  xge_3Dwinf *_3Dwin;
  point3f    a, b, c;
  float      w, d, h;
  vector3f   shift;
  int        id;

  id = er->id & 0x03;
  if ( id != 3 )
    return;
  _3Dwin = er->data0;
  if ( x != _3Dwin->xx || y != _3Dwin->yy ) {
    SetPoint3f ( &c, (float)(0.5*(_3Dwin->PerspBBox.x1+_3Dwin->PerspBBox.x0)),
                     (float)(0.5*(_3Dwin->PerspBBox.y1+_3Dwin->PerspBBox.y0)),
                     (float)(0.5*(_3Dwin->PerspBBox.z1+_3Dwin->PerspBBox.z0)) );
    CameraProjectPoint3f ( &_3Dwin->CPos[3], &c, &b );
    SetPoint3f ( &a, (float)_3Dwin->xx, (float)_3Dwin->yy, b.z );
    CameraUnProjectPoint3f ( &_3Dwin->CPos[3], &a, &a );
    SetPoint3f ( &b, (float)x, (float)y, b.z );
    CameraUnProjectPoint3f ( &_3Dwin->CPos[3], &b, &b );
    SubtractPoints3f ( &a, &b, &shift );
/* bound the centre --- to be added */
    w = (float)(0.5*(_3Dwin->PerspBBox.x1-_3Dwin->PerspBBox.x0));
    h = (float)(0.5*(_3Dwin->PerspBBox.y1-_3Dwin->PerspBBox.y0));
    d = (float)(0.5*(_3Dwin->PerspBBox.z1-_3Dwin->PerspBBox.z0));
    AddVector3f ( &c, &shift, &c );
    _3Dwin->PerspBBox.x0 = c.x-w;  _3Dwin->PerspBBox.x1 = c.x+w;
    _3Dwin->PerspBBox.y0 = c.y-h;  _3Dwin->PerspBBox.x1 = c.y+h;
    _3Dwin->PerspBBox.z0 = c.z-d;  _3Dwin->PerspBBox.x1 = c.z+d;
    CameraMoveGf ( &_3Dwin->CPos[3], &shift );
    xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, 3, x, y );
    _3Dwin->xx = x;
    _3Dwin->yy = y;
  }
} /*xge_3DwinfPanPerspWindow*/

static void xge_3DwinfDrawXYZ ( CameraRecf *CPos, char cc, char left,
                                int xi, int eta, float xyz )
{
  char s[30];

  xgeSetForeground ( xgec_Green );
  switch ( cc ) {
case 0: sprintf ( s, "x=%5.3f", xyz );  break;
case 1: sprintf ( s, "y=%5.3f", xyz );  break;
case 2: sprintf ( s, "z=%5.3f", xyz );  break;
  }
  if ( left ) {
    xi = CPos->xmin+2;
    eta = max ( eta, CPos->ymin+11 );
  }
  else {
    xi = min ( xi+2, CPos->xmin+CPos->width-6*strlen(s) );
    eta = CPos->ymin+CPos->height-4;
  }
  xgeDrawString ( s, xi, eta );
} /*xge_3DwinfDrawXYZ*/

void xge_3DwinfDrawCursorPos ( xge_3Dwinf *_3Dwin,
                               int id, short x, short y )
{
  point3f p, q;
  char    zoomwin;

  id &= 0x03;
  if ( !_3Dwin->display_coord ||
       _3Dwin->CoordWin < 0 || _3Dwin->CoordWin > 2 ||
       id < 0 || id > 2 )
    return;
  zoomwin = _3Dwin->fww.zoomwin;
  xgeSetForeground ( xgec_Green5 );
  switch ( _3Dwin->CoordWin ) {
case 0:
    switch ( id ) {
  case 0:
      xgeDrawLine ( x, _3Dwin->CPos[id].ymin,
                    x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      xgeDrawLine ( _3Dwin->CPos[id].xmin, y,
                    _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y ); 
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xge_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 2, 1, x, y, p.z );
      if ( zoomwin == 0 )
        xge_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 0, 0, x, y, p.x );
      break;
  case 1:   
      xgeDrawLine ( _3Dwin->CPos[id].xmin, y,
                    _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y );
      break;
  case 2:   
      xgeDrawLine ( x, _3Dwin->CPos[id].ymin,
                    x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xge_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 0, 0, x, y, p.x );
      break;
    }
    break;
case 1:   
    switch ( id ) {
  case 0:
      xgeDrawLine ( _3Dwin->CPos[id].xmin, y,
                    _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y );
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xge_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 2, 1, x, y, p.z );
      break;
  case 1:   
      xgeDrawLine ( x, _3Dwin->CPos[id].ymin,
                    x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      xgeDrawLine ( _3Dwin->CPos[id].xmin, y,
                    _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y ); 
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xge_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 1, 0, x, y, p.y );
      if ( zoomwin == 1 )
        xge_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 2, 1, x, y, p.z );
      break;
  case 2:   
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[_3Dwin->CoordWin], &q, &p );
      CameraProjectPoint3f ( &_3Dwin->CPos[id], &p, &q );   
      xgeDrawLine ( _3Dwin->CPos[id].xmin, (int)q.y,
                    _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, (int)q.y );
      break;
    }
    break;
case 2:   
    switch ( id ) {
  case 0:
      xgeDrawLine ( x, _3Dwin->CPos[id].ymin,
                    x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      break;
  case 1:   
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[_3Dwin->CoordWin], &q, &p );
      CameraProjectPoint3f ( &_3Dwin->CPos[id], &p, &q );   
      xgeDrawLine ( (int)q.x, _3Dwin->CPos[id].ymin,
                    (int)q.x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      break;
  case 2:   
      xgeDrawLine ( x, _3Dwin->CPos[id].ymin,
                    x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      xgeDrawLine ( _3Dwin->CPos[id].xmin, y,
                    _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y ); 
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xge_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 1, 1, x, y, p.y );
      xge_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 0, 0, x, y, p.x );
      break;
    }
    break;
  }
} /*xge_3DwinfDrawCursorPos*/

void xge_3DwinfEnableGeomWidget ( xge_3Dwinf *_3Dwin, char tool )
{
  _3Dwin->selecting_mode = _3Dwin->special_selecting_mode =
  _3Dwin->panning =
  _3Dwin->moving_tool = _3Dwin->scaling_tool =
  _3Dwin->rotating_tool = _3Dwin->shear_tool =
  _3Dwin->special_trans_tool = false;
  switch ( tool ) {
case xge_3DWIN_MOVING_TOOL:
    _3Dwin->moving_tool = true;
    goto go_on;
case xge_3DWIN_SCALING_TOOL:
    _3Dwin->scaling_tool = true;
    goto go_on;
case xge_3DWIN_ROTATING_TOOL:
    _3Dwin->rotating_tool = true;
    goto go_on;
case xge_3DWIN_SHEAR_TOOL:
    _3Dwin->shear_tool = true;
    goto go_on;
case xge_3DWIN_SPECIAL_TRANS_TOOL:
    _3Dwin->special_trans_tool = true;
    goto go_on;
case xge_3DWIN_SELECTING_TOOL:
    _3Dwin->selecting_mode = true;
    goto go_on;
case xge_3DWIN_SPECIAL_SELECTING_TOOL:
    _3Dwin->special_selecting_mode = true;
    goto go_on;
case xge_2DWIN_PANNING_TOOL:
    _3Dwin->panning = true;
go_on:
    _3Dwin->current_tool = tool;
    break;
default:
    _3Dwin->current_tool = xge_3DWIN_NO_TOOL;
    break;
  }
} /*xge_3DwinfEnableGeomWidget*/

void xge_3DwinfResetGeomWidgets ( xge_3Dwinf *_3Dwin )
{
  SetVector3f ( &_3Dwin->scaling_factors, 1.0, 1.0, 1.0 );
  _3Dwin->scaling_size = 1;
  _3Dwin->tool_mode = 0;
  _3Dwin->rotating_radius = 0;
  SetVector3f ( &_3Dwin->shear_axis[0], 1.0, 0.0, 0.0 );
  SetVector3f ( &_3Dwin->shear_axis[1], 0.0, 1.0, 0.0 );
  SetVector3f ( &_3Dwin->shear_axis[2], 0.0, 0.0, 1.0 );
} /*xge_3DwinfResetGeomWidgets*/

void xge_3DwinfResetGeomWidgetPos ( xge_3Dwinf *_3Dwin )
{
  SetPoint3f ( &_3Dwin->scaling_centre, 0.0, 0.0, 0.0 );
  SetPoint3f ( &_3Dwin->rotating_centre, 0.0, 0.0, 0.0 );
  SetPoint3f ( &_3Dwin->shear_centre, 0.0, 0.0, 0.0 );
  xge_3DwinfResetGeomWidgets ( _3Dwin );
} /*xge_3DwinfResetGeomWidgetPos*/

static short _xge_GetWinfRadius ( xge_3Dwinf *_3Dwin )
{
  int   i;
  short r;
  CameraRecf *CPos;

  CPos = _3Dwin->CPos;
  i = _3Dwin->fww.zoomwin;
  if ( i >= 0 && i < 4 )
    r = (short)(min ( CPos[i].width, CPos[i].height ));
  else {
    r = (short)(min ( CPos[0].width, CPos[0].height ));
    for ( i = 1; i < 3; i++ ) {
      r = (short)(min ( r, CPos[i].width ));
      r = (short)(min ( r, CPos[i].height ));
    }
  }
  return r;
} /*_xge_GetWinfRadius*/

void xge_3DwinfDrawGeomWidgets ( xge_widget *er )
{
  xge_3Dwinf  *_3Dwin;
  CameraRecf  *CPos;
  int         id, wcase, i, s, t;
  short       xc, yc, r, xp, yp;
  static char dl[2] = {2,3};
  float       sx, sy;
  point3f     p, q;

  if ( er->state == xgestate_3DWIN_MOVINGPOINT )
    return;
  if ( (id = er->id & 0x03) < 3 ) {
    _3Dwin = er->data0;
    CPos = &_3Dwin->CPos[id];
    xgeSetForeground ( xgec_Cyan );
    switch ( _3Dwin->current_tool ) {
  case xge_3DWIN_MOVING_TOOL:
      xc = (short)(CPos->xmin + (float)(CPos->width)/2.0 + 0.5);
      yc = (short)(CPos->ymin + (float)(CPos->height)/2.0 + 0.5);
      xgeFillRectangle ( 5, 3, xc-2, yc-1 );
      xgeFillRectangle ( 3, 5, xc-1, yc-2 );
      break;

  case xge_3DWIN_SCALING_TOOL:
      if ( er->state == xgestate_3DWIN_USING_GEOM_WIDGET ||
           er->state == xgestate_3DWIN_ALTUSING_GEOM_WIDGET ) {
        wcase = 3*_3Dwin->CoordWin + id;
        switch ( wcase ) {
      case 0:  case 4:  case 8:
          sx = _3Dwin->scaling_factors.x;  sy = _3Dwin->scaling_factors.y;  break;
      case 1:  case 3:
          sx = _3Dwin->scaling_factors.z;  sy = _3Dwin->scaling_factors.y;  break;
      case 2:  case 6:
          sx = _3Dwin->scaling_factors.x;  sy = _3Dwin->scaling_factors.z;  break;
      case 5:
          sx = _3Dwin->scaling_factors.z;  sy = _3Dwin->scaling_factors.x;  break;
      case 7:
          sx = _3Dwin->scaling_factors.y;  sy = _3Dwin->scaling_factors.z;  break;
     default:
          return;
        }
      }
      else
        sx = sy = 1.0;
      CameraProjectPoint3f ( CPos, &_3Dwin->scaling_centre, &p );
      xc = (short)(p.x+0.5);
      yc = (short)(p.y+0.5);
      r = _xge_GetWinfRadius ( _3Dwin );
      _3Dwin->scaling_size = r = (short)((3*r)/10);
      xgeFillRectangle ( 5, 3, xc-2, yc-1 );
      xgeFillRectangle ( 3, 5, xc-1, yc-2 );
      xgeFillRectangle ( 3, 3, xc-1, yc-r*sy-1 );
      xgeFillRectangle ( 3, 3, xc-1, yc+r*sy-1 );
      xgeFillRectangle ( 3, 3, xc-r*sx-1, yc-1 );
      xgeFillRectangle ( 3, 3, xc+r*sx-1, yc-1 );
      xgeFillRectangle ( 3, 3, xc-r*sx-1, yc-r*sy-1 );
      xgeFillRectangle ( 3, 3, xc-r*sx-1, yc+r*sy-1 );
      xgeFillRectangle ( 3, 3, xc+r*sx-1, yc-r*sy-1 );
      xgeFillRectangle ( 3, 3, xc+r*sx-1, yc+r*sy-1 );
      break;

  case xge_3DWIN_ROTATING_TOOL:
      CameraProjectPoint3f ( CPos, &_3Dwin->rotating_centre, &p );
      xc = (short)(p.x+0.5);
      yc = (short)(p.y+0.5);
      r = _xge_GetWinfRadius ( _3Dwin );
      _3Dwin->rotating_radius = r = (short)((3*r)/10);
      xgeFillRectangle ( 5, 3, xc-2, yc-1 );
      xgeFillRectangle ( 3, 5, xc-1, yc-2 );
      xgeSetDashes ( 2, dl, 0 );
      xgeSetLineAttributes ( 1, LineOnOffDash, CapRound, JoinRound );
      xgeDrawArc ( 2*r+1, 2*r+1, xc-r, yc-r, 0, 360*64 );
      xgeSetLineAttributes ( 1, LineSolid, CapRound, JoinRound );
      break;

  case xge_3DWIN_SHEAR_TOOL:
      CameraProjectPoint3f ( &CPos[0], &_3Dwin->shear_centre, &p );
      xc = (short)(p.x+0.5);
      yc = (short)(p.y+0.5);
      r = _xge_GetWinfRadius ( _3Dwin );
      s = 0;
      for ( i = 0; i < 3; i++ ) {
        q = _3Dwin->shear_centre;
        switch ( i ) {
      case 0: q.x += 1.0;  break;
      case 1: q.y += 1.0;  break;
      case 2: q.z += 1.0;  break;
        }
        CameraProjectPoint3f ( &CPos[0], &q, &q );
        xp = (short)(q.x+0.5) - xc;
        yp = (short)(q.y+0.5) - yc;
        t = xp*xp + yp*yp;
        s = max ( s, t ); 
      }
      s = (int)(sqrt(s)+0.5);
      _3Dwin->shear_radius = (float)(0.3*(float)r/(float)s);
      CameraProjectPoint3f ( CPos, &_3Dwin->shear_centre, &p );
      xc = (short)(p.x+0.5);
      yc = (short)(p.y+0.5);
      xgeFillRectangle ( 5, 3, xc-2, yc-1 );
      xgeFillRectangle ( 3, 5, xc-1, yc-2 );
      for ( i = 0; i < 3; i++ ) {
        AddVector3Mf ( &_3Dwin->shear_centre, &_3Dwin->shear_axis[i],
                       _3Dwin->shear_radius, &q );
        CameraProjectPoint3f ( CPos, &q, &q );
        xc = (short)(q.x+0.5);
        yc = (short)(q.y+0.5);
        xgeFillRectangle ( 3, 3, xc-1, yc-1 );
        AddVector3Mf ( &_3Dwin->shear_centre, &_3Dwin->shear_axis[i],
                       -_3Dwin->shear_radius, &q );
        CameraProjectPoint3f ( CPos, &q, &q );
        xp = (short)(q.x+0.5);
        yp = (short)(q.y+0.5);
        xgeFillRectangle ( 3, 3, xp-1, yp-1 );
        xgeSetDashes ( 2, dl, 0 );
        xgeSetLineAttributes ( 1, LineOnOffDash, CapRound, JoinRound );
        xgeDrawLine ( xc, yc, xp, yp );
        xgeSetLineAttributes ( 1, LineSolid, CapRound, JoinRound );
      }
      break;

  default:
      break;
    }
  }
} /*xge_3DwinfDrawGeomWidgets*/

char xge_3DwinfIsItAGeomWidget ( xge_3Dwinf *_3Dwin, int id,
                                 int key, short x, short y )
{
/* return value 0 - no widget has been pointed */
/* 1 - the widget should be applied */
/* 2 - the widget should be applied in the alternative mode */
/* 3 - the widget should be moved */
  short      xc, yc, r;
  CameraRecf *CPos;
  point3f    p;
  char       result;
  int        i, d, e;

  result = 0;
  id &= 0x03;
  CPos = &_3Dwin->CPos[id];
  switch ( _3Dwin->current_tool ) {
case xge_3DWIN_MOVING_TOOL:
    xc = (short)(CPos->xmin + (float)(CPos->width)/2.0 + 0.5);
    yc = (short)(CPos->ymin + (float)(CPos->height)/2.0 + 0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST )
      result = 1;
    break;

case xge_3DWIN_SCALING_TOOL:
    CameraProjectPoint3f ( CPos, &_3Dwin->scaling_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _3Dwin->saved_centre = _3Dwin->scaling_centre;
      result = 3;
    }
    else {
      if ( key & xgemouse_LBUTTON_DOWN )
        result = 1;
      else
        result = 2;
      r = _3Dwin->scaling_size;
      if ( abs(x-(xc-r))+abs(y-(yc-r)) < xge_MINDIST )
        _3Dwin->tool_mode = 1;
      else if ( abs(x-xc)+abs(y-(yc-r)) < xge_MINDIST )
        _3Dwin->tool_mode = 2;
      else if ( abs(x-(xc+r))+abs(y-(yc-r)) < xge_MINDIST )
        _3Dwin->tool_mode = 3;
      else if ( abs(x-(xc-r))+abs(y-yc) < xge_MINDIST )
        _3Dwin->tool_mode = 4;
      else if ( abs(x-(xc+r))+abs(y-yc) < xge_MINDIST )
        _3Dwin->tool_mode = 5;
      else if ( abs(x-(xc-r))+abs(y-(yc+r)) < xge_MINDIST )
        _3Dwin->tool_mode = 6;
      else if ( abs(x-xc)+abs(y-(yc+r)) < xge_MINDIST )
        _3Dwin->tool_mode = 7;
      else if ( abs(x-(xc+r))+abs(y-(yc+r)) < xge_MINDIST )
        _3Dwin->tool_mode = 8;
      else result = 0;
    }
    break;

case xge_3DWIN_ROTATING_TOOL:
    CameraProjectPoint3f ( CPos, &_3Dwin->rotating_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _3Dwin->saved_centre = _3Dwin->rotating_centre;
      result = 3;
    }
    else {
      r = (short)sqrt ( (x-xc)*(x-xc) + (y-yc)*(y-yc) );
      if ( abs ( r-_3Dwin->rotating_radius ) < xge_MINDIST )
        result = 1;
    }
    break;

 case xge_3DWIN_SHEAR_TOOL:
    CameraProjectPoint3f ( CPos, &_3Dwin->shear_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( abs(x-xc)+abs(y-yc) < xge_MINDIST ) {
      _3Dwin->saved_centre = _3Dwin->shear_centre;
      _3Dwin->tool_mode = 0;
      result = 3;
    }
    else {
      e = xge_MINDIST+1;
      for ( i = 0; i < 3; i++ ) {
        AddVector3Mf ( &_3Dwin->shear_centre, &_3Dwin->shear_axis[i],
                       _3Dwin->shear_radius, &p );
        CameraProjectPoint3f ( CPos, &p, &p );
        xc = (short)(p.x+0.5);
        yc = (short)(p.y+0.5);
        d = abs(xc-x) + abs(yc-y);
        if ( d < e ) {
          e = d;
          _3Dwin->tool_mode = 2*i+1;
          _3Dwin->saved_centre = _3Dwin->shear_axis[i];
        }
        AddVector3Mf ( &_3Dwin->shear_centre, &_3Dwin->shear_axis[i],
                       -_3Dwin->shear_radius, &p );
        CameraProjectPoint3f ( CPos, &p, &p );
        xc = (short)(p.x+0.5);
        yc = (short)(p.y+0.5);
        d = abs(xc-x) + abs(yc-y);
        if ( d < e ) {
          e = d;
          _3Dwin->tool_mode = 2*i+2;
          _3Dwin->saved_centre = _3Dwin->shear_axis[i];
        }
      }
      if ( e < xge_MINDIST ) {
        if ( key & xgemouse_LBUTTON_DOWN )
          result = 1;
        else if ( key & xgemouse_RBUTTON_DOWN )
          result = 3;
      }
    }
    break;

 default:
    break;
  }
  return result;
} /*xge_3DwinfIsItAGeomWidget*/

static void _xge_3DwinfMoveShearFrame ( vector3f *sa1,
                                        vector3f *a1, vector3f *a2, vector3f *a3,
                                        vector3f *da1, boolean neg )
{
#define EPSANG 1.0e-6
  vector3f na1, v, aa1;
  float    d, c, s, phi;
  trans3f  tr;

  na1 = *sa1;
  if ( neg )
    { da1->x = -da1->x;  da1->y = -da1->y;  da1->z = -da1->z; }
  AddVector3f ( sa1, da1, &aa1 );
  d = sqrt ( DotProduct3f ( &aa1, &aa1 ) );
  if ( d < EPSANG )
    return;
  CrossProduct3f ( &na1, &aa1, &v );
  s = sqrt ( DotProduct3f ( &v, &v ) );
  if ( s > EPSANG ) {
    c = DotProduct3f ( &na1, sa1 );
    phi = atan2 ( s, c );
    NormalizeVector3f ( &v );
    IdentTrans3f ( &tr );
    RotVTrans3f ( &tr, &v, phi );
    TransVector3f ( &tr, a2, a2 );
    TransVector3f ( &tr, a3, a3 );
    *a1 = aa1;
    OrtVector3f ( a1, a2, a2 );
    OrtVector3f ( a1, a3, a3 );
    OrtVector3f ( a2, a3, a3 );
    NormalizeVector3f ( a2 );
    NormalizeVector3f ( a3 );
  }
#undef EPSANG
} /*_xge_3DwinfMoveShearFrame*/

void xge_3DwinfMoveGeomWidget ( xge_3Dwinf *_3Dwin, int id, short x, short y )
{
  CameraRecf *CPos;
  point3f    p, q, r;
  vector3f   v;

  id &= 0x03;
  CPos = &_3Dwin->CPos[id];
  SetPoint3f ( &p, _3Dwin->xx, _3Dwin->yy, 0.0 );
  SetPoint3f ( &q, x, y, 0.0 );
  CameraUnProjectPoint3f ( CPos, &q, &r );
  CameraUnProjectPoint3f ( CPos, &p, &q );
  SubtractPoints3f ( &r, &q, &v );
  switch ( _3Dwin->current_tool ) {
case xge_3DWIN_SCALING_TOOL:
    AddVector3f ( &_3Dwin->saved_centre, &v, &_3Dwin->scaling_centre );
    break;
case xge_3DWIN_ROTATING_TOOL:
    AddVector3f ( &_3Dwin->saved_centre, &v, &_3Dwin->rotating_centre );
    break;
case xge_3DWIN_SHEAR_TOOL:
    if ( _3Dwin->tool_mode == 0 )
      AddVector3f ( &_3Dwin->saved_centre, &v, &_3Dwin->shear_centre );
    else {
      MultVector3f ( 1.0/_3Dwin->shear_radius, &v, &v );
      switch ( _3Dwin->tool_mode ) {
    case 1:
        _xge_3DwinfMoveShearFrame ( &_3Dwin->saved_centre, &_3Dwin->shear_axis[0],
            &_3Dwin->shear_axis[1], &_3Dwin->shear_axis[2], &v, false );
        break;
    case 2:
        _xge_3DwinfMoveShearFrame ( &_3Dwin->saved_centre, &_3Dwin->shear_axis[0],
            &_3Dwin->shear_axis[1], &_3Dwin->shear_axis[2], &v, true );
        break;
    case 3:
        _xge_3DwinfMoveShearFrame ( &_3Dwin->saved_centre, &_3Dwin->shear_axis[1],
            &_3Dwin->shear_axis[2], &_3Dwin->shear_axis[0], &v, false );
        break;
    case 4:
        _xge_3DwinfMoveShearFrame ( &_3Dwin->saved_centre, &_3Dwin->shear_axis[1],
            &_3Dwin->shear_axis[2], &_3Dwin->shear_axis[0], &v, true );
        break;
    case 5:
        _xge_3DwinfMoveShearFrame ( &_3Dwin->saved_centre, &_3Dwin->shear_axis[2],
            &_3Dwin->shear_axis[0], &_3Dwin->shear_axis[1], &v, false );
        break;
    case 6:
        _xge_3DwinfMoveShearFrame ( &_3Dwin->saved_centre, &_3Dwin->shear_axis[2],
            &_3Dwin->shear_axis[0], &_3Dwin->shear_axis[1], &v, true );
        break;
    default:
        break;
      }
    }
    break;
default:
    break;
  }
  xge_callback ( _3Dwin->fww.er,
                 xgemsg_3DWIN_CHANGE_TRANS, _3Dwin->current_tool, 0, 0 );
} /*xge_3DwinfMoveGeomWidget*/

boolean xge_3DwinfApplyGeomWidget ( xge_3Dwinf *_3Dwin, int id, short x, short y,
                                    boolean alt )
{
  CameraRecf *CPos;
  point3f    p, q, r;
  vector3f   v, *sf;
  float      alpha, beta;
  short      xx, yy, xc, yc, s;
  trans3f    tr1, tr2, tr3;

  id &= 0x03;
  CPos = &_3Dwin->CPos[id];
  IdentTrans3f ( &_3Dwin->gwtrans );
  switch ( _3Dwin->current_tool ) {

case xge_3DWIN_MOVING_TOOL:
    SetPoint3f ( &p, _3Dwin->xx, _3Dwin->yy, 0.0 );
    SetPoint3f ( &q, x, y, 0.0 );
    CameraUnProjectPoint3f ( CPos, &q, &r );
    CameraUnProjectPoint3f ( CPos, &p, &q );
    SubtractPoints3f ( &r, &q, &v );
    ShiftTrans3f ( &_3Dwin->gwtrans, v.x, v.y, v.z );
    _3Dwin->trans_params = v;
    break;

case xge_3DWIN_SCALING_TOOL:
    sf = &_3Dwin->scaling_factors;
    CameraProjectPoint3f ( CPos, &_3Dwin->scaling_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    s = _3Dwin->scaling_size;
    switch ( _3Dwin->tool_mode ) {
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
    sf->z = 1.0;

    if ( alt ) {
      switch ( _3Dwin->tool_mode ) {
  case 1:  case 3:  case 6:  case 8:
        sf->x = sf->y = (float)(0.5*(sf->x+sf->y));
        break;
  case 2:  case 7:
        sf->x = sf->z = sf->y;
        break;
  case 4:  case 5:
        sf->y = sf->z = sf->x;
        break;
      }
    }

    ShiftTrans3f ( &_3Dwin->gwtrans, -_3Dwin->scaling_centre.x,
                   -_3Dwin->scaling_centre.y, -_3Dwin->scaling_centre.z );
    EulerRotTrans3f ( &_3Dwin->gwtrans, -CPos->psi, -CPos->theta, -CPos->phi );
    ScaleTrans3f ( &_3Dwin->gwtrans, sf->x, sf->y, sf->z );
    EulerRotTrans3f ( &_3Dwin->gwtrans, CPos->phi, CPos->theta, CPos->psi );
    ShiftTrans3f ( &_3Dwin->gwtrans, _3Dwin->scaling_centre.x,
                   _3Dwin->scaling_centre.y, _3Dwin->scaling_centre.z );
    _3Dwin->trans_params = *sf;
    break;

case xge_3DWIN_ROTATING_TOOL:
    SetPoint3f ( &p, 0.0, 0.0, 0.0 );
    SetPoint3f ( &q, 0.0, 0.0, 1.0 );
    CameraUnProjectPoint3f ( CPos, &q, &r );
    CameraUnProjectPoint3f ( CPos, &p, &q );
    SubtractPoints3f ( &r, &q, &v );
    NormalizeVector3f ( &v );
    CameraProjectPoint3f ( CPos, &_3Dwin->rotating_centre, &p );
    xc = (short)(p.x+0.5);
    yc = (short)(p.y+0.5);
    if ( ((_3Dwin->xx != xc) || (_3Dwin->yy != yc)) && ((x != xc) || (y != yc))) {
      beta = (float)atan2 ( _3Dwin->yy-yc, _3Dwin->xx-xc );
      alpha = (float)(atan2 ( y-yc, x-xc ) - beta);
      ShiftTrans3f ( &_3Dwin->gwtrans, -_3Dwin->rotating_centre.x,
                     -_3Dwin->rotating_centre.y, -_3Dwin->rotating_centre.z );
      RotVTrans3f ( &_3Dwin->gwtrans, &v, alpha );
      ShiftTrans3f ( &_3Dwin->gwtrans, _3Dwin->rotating_centre.x,
                     _3Dwin->rotating_centre.y, _3Dwin->rotating_centre.z );
      _3Dwin->trans_params.x = alpha;
    }
    else
      return false;
    break;

case xge_3DWIN_SHEAR_TOOL:  
    SetPoint3f ( &p, _3Dwin->xx, _3Dwin->yy, 0.0 );
    SetPoint3f ( &q, x, y, 0.0 );
    CameraUnProjectPoint3f ( CPos, &p, &p );
    CameraUnProjectPoint3f ( CPos, &q, &q );
    SubtractPoints3f ( &q, &p, &v );
    MultVector3f ( 1.0/_3Dwin->shear_radius, &v, &v );
    switch ( _3Dwin->tool_mode ) {
  case 2:
      v.x = -v.x;  v.y = -v.y;  v.z = -v.z;
  case 1:
      AddVector3f ( &_3Dwin->saved_centre, &v, &_3Dwin->shear_axis[0] );
      break;
  case 4:
      v.x = -v.x;  v.y = -v.y;  v.z = -v.z;
  case 3:
      AddVector3f ( &_3Dwin->saved_centre, &v, &_3Dwin->shear_axis[1] );
      break;
  case 6:
      v.x = -v.x;  v.y = -v.y;  v.z = -v.z;
  case 5:
      AddVector3f ( &_3Dwin->saved_centre, &v, &_3Dwin->shear_axis[2] );
      break;
  default:
      return false;
    }
    tr1.U0.a11 = tr2.U0.a11 = _3Dwin->shear_axis[0].x;
    tr1.U0.a21 = tr2.U0.a12 = _3Dwin->shear_axis[0].y;
    tr1.U0.a31 = tr2.U0.a13 = _3Dwin->shear_axis[0].z;
    tr1.U0.a12 = tr2.U0.a21 = _3Dwin->shear_axis[1].x;
    tr1.U0.a22 = tr2.U0.a22 = _3Dwin->shear_axis[1].y;
    tr1.U0.a32 = tr2.U0.a23 = _3Dwin->shear_axis[1].z;
    tr1.U0.a13 = tr2.U0.a31 = _3Dwin->shear_axis[2].x;
    tr1.U0.a23 = tr2.U0.a32 = _3Dwin->shear_axis[2].y;
    tr1.U0.a33 = tr2.U0.a33 = _3Dwin->shear_axis[2].z;
    tr1.U0.a14 = tr1.U0.a24 = tr1.U0.a34 =
    tr2.U0.a14 = tr2.U0.a24 = tr2.U0.a34 = 0.0;
    switch ( _3Dwin->tool_mode ) {
  case 1:  case 2:
      CrossProduct3f ( &_3Dwin->shear_axis[1], &_3Dwin->shear_axis[2], &v );
      tr2.U0.a11 = v.x;  tr2.U0.a12 = v.y;  tr2.U0.a13 = v.z;
      break;
  case 3:  case 4:
      CrossProduct3f ( &_3Dwin->shear_axis[2], &_3Dwin->shear_axis[0], &v );
      tr2.U0.a21 = v.x;  tr2.U0.a22 = v.y;  tr2.U0.a23 = v.z;
      break;
  case 5:  case 6:
      CrossProduct3f ( &_3Dwin->shear_axis[0], &_3Dwin->shear_axis[1], &v );
      tr2.U0.a31 = v.x;  tr2.U0.a32 = v.y;  tr2.U0.a33 = v.z;
      break;
    }
    ShiftTrans3f ( &_3Dwin->gwtrans, -_3Dwin->shear_centre.x,
                   -_3Dwin->shear_centre.y, -_3Dwin->shear_centre.z );
    CompTrans3f ( &tr3, &tr2, &_3Dwin->gwtrans );
    CompTrans3f ( &_3Dwin->gwtrans, &tr1, &tr3 );
    ShiftTrans3f ( &_3Dwin->gwtrans, _3Dwin->shear_centre.x,
                   _3Dwin->shear_centre.y, _3Dwin->shear_centre.z );
    break;

default:
    return false;
  }
  return true;
} /*xge_3DwinfApplyGeomWidget*/

void xge_3DwinfExitWidgetMode ( xge_3Dwinf *_3Dwin )
{
  if ( _3Dwin->current_tool == xge_3DWIN_SHEAR_TOOL ) {
    switch ( _3Dwin->tool_mode ) {
  case 1:  case 2:
      CrossProduct3f ( &_3Dwin->shear_axis[1], &_3Dwin->shear_axis[2],  
                       &_3Dwin->shear_axis[0] );
      NormalizeVector3f ( &_3Dwin->shear_axis[0] );
      break;
  case 3:  case 4:
      CrossProduct3f ( &_3Dwin->shear_axis[2], &_3Dwin->shear_axis[0],
                       &_3Dwin->shear_axis[1] );
      NormalizeVector3f ( &_3Dwin->shear_axis[1] );
      break;
  case 5:  case 6:
      CrossProduct3f ( &_3Dwin->shear_axis[0], &_3Dwin->shear_axis[1],
                       &_3Dwin->shear_axis[2] );
      NormalizeVector3f ( &_3Dwin->shear_axis[2] );
      break;
  default:
      break;
    }
  }
} /*xge_3DwindExitWidgetMode*/

void xge_3DwinfResetGeomWidget ( xge_3Dwinf *_3Dwin )
{
  xge_widget *er;

  switch ( _3Dwin->current_tool ) {
case xge_3DWIN_SCALING_TOOL:
    SetPoint3f ( &_3Dwin->scaling_centre, 0.0, 0.0, 0.0 );
    break;
case xge_3DWIN_ROTATING_TOOL:
    SetPoint3f ( &_3Dwin->rotating_centre, 0.0, 0.0, 0.0 );
    break;
case xge_3DWIN_SHEAR_TOOL:
    if ( _3Dwin->tool_mode == 0 )
      SetPoint3f ( &_3Dwin->shear_centre, 0.0, 0.0, 0.0 );
    else {
      SetVector3f ( &_3Dwin->shear_axis[0], 1.0, 0.0, 0.0 );
      SetVector3f ( &_3Dwin->shear_axis[1], 0.0, 1.0, 0.0 );
      SetVector3f ( &_3Dwin->shear_axis[2], 0.0, 0.0, 1.0 );
    }
    break;
default:
    return;
  }
  er = xge_GetFocusWidget ( _3Dwin->fww.er->window_num );
  if ( er )
    er->state = xgestate_NOTHING;
  xge_ReleaseFocus ( _3Dwin->fww.er );
  xge_SetClipping ( _3Dwin->fww.er );
  _3Dwin->fww.er->redraw ( _3Dwin->fww.er, true );
} /*xge_3DwindResetGeomWidget*/

void xge_3DwinfSavePerspCamera ( xge_3Dwinf *_3Dwin )
{
  _3Dwin->CPos[4] = _3Dwin->CPos[3];
} /*xge_3DwinfSavePerspCamera*/

void xge_3DwinfSwapPerspCameras ( xge_3Dwinf *_3Dwin )
{
  CameraRecf *CPos;

  CPos = _3Dwin->CPos;
  pkv_Exchange ( &CPos[3], &CPos[4], sizeof(CameraRecf) );
  CameraInitFramef ( &CPos[3], false, false, CPos[4].width, CPos[4].height,
                     CPos[4].xmin, CPos[4].ymin,
                     CPos[4].aspect, CPos[4].ncplanes );
  CameraSetMappingf ( &CPos[3] );
} /*xge_3DwinfSwapPerspCameras*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean xge_3DwinfPerspMsg ( xge_widget *er,
                                    int msg, int key, short x, short y )
{
  xge_3Dwinf *_3Dwin;
  CameraRecf *CPos;
  int        id;
  vector3f   v;
  float      a, angle;

  _3Dwin = er->data0;
  id = er->id & 0x03;            /* id == 3 anyway */
  CPos = &_3Dwin->CPos[id];
  switch ( msg ) {
case xgemsg_ENTERING:
    if ( _3Dwin->panning )
      xge_SetCurrentWindowCursor ( xgeCURSOR_FLEUR );
    else if ( _3Dwin->selecting_mode || _3Dwin->special_selecting_mode )
      xge_SetCurrentWindowCursor ( xgeCURSOR_ARROW );
    else
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    _3Dwin->CoordWin = (char)id;
    return true;

case xgemsg_EXITING:
    if ( _3Dwin->panning || _3Dwin->selecting_mode ||
         _3Dwin->special_selecting_mode )
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    _3Dwin->CoordWin = -1;
    return true;

case xgemsg_SPECIAL_KEY:
    return false;

default:
    break;
  }

  switch ( er->state ) {
case xgestate_3DWIN_TURNING_VIEWER:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( x != _3Dwin->xx || y != _3Dwin->yy ) {
          SetVector3f ( &v, (y-_3Dwin->yy)/CPos->yscale,
                        (_3Dwin->xx-x)/CPos->xscale, 0.0 );
          a = (float)(1000.0*sqrt ( ((float)v.x*v.x + v.y*v.y)/
                    ((float)CPos->width*CPos->width +
                     (float)CPos->height*CPos->height) ));
          angle = (float)(-(2.0*a+1.0)*a);
          CameraRotVCf ( CPos, &v, angle );
          xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, id, x, y );
          xge_SetClipping ( er ); 
          er->redraw ( er, true );
          _3Dwin->xx = x;
          _3Dwin->yy = y;
        }
      }
      else
        goto exit_mode;
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_ZOOMING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( key & xgemouse_RBUTTON_DOWN ) {
        if ( y != _3Dwin->yy ) {
          a = (float)(_3Dwin->xx-y)/(float)CPos->height;
          a = (float)(_3Dwin->perspzoom*exp(a));
          a = (float)max ( a, xge_3DWIN_MIN_ZOOM );
          a = (float)min ( a, xge_3DWIN_MAX_ZOOM );
          CameraSetFf ( CPos, a );
          xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, id, x, y );
          xge_SetClipping ( er );
          er->redraw ( er, true );
          _3Dwin->yy = y;
        }
      }
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) )
        goto exit_mode;
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_PANNING:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      break;
  case xgemsg_MMOVE:
      if ( !(key & xgemouse_LBUTTON_DOWN) )
        goto exit_mode;
      xge_3DwinfPanPerspWindow ( er, x, y );
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  default:
      break;
    }
    break;

case xgestate_RESIZING_X:
case xgestate_RESIZING_Y:
case xgestate_RESIZING_XY:
#ifdef SHOW_SIZE
    if ( msg == xgemsg_RESIZE )
      printf ( "w=%d, h=%d\n", er->w, er->h );
#endif
    break;

case xgestate_3DWIN_SELECTING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      _3Dwin->selection_rect.x1 = x;
      _3Dwin->selection_rect.y1 = y;
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        xge_OrderSelectionRect ( &_3Dwin->selection_rect );
        xge_callback ( er, xgemsg_3DWIN_SELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        _3Dwin->selection_rect.x1 = x;
        _3Dwin->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &_3Dwin->selection_rect );
        xge_callback ( er, xgemsg_3DWIN_SELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      break;
  default:
      break;
    }
    break;

case xgestate_3DWIN_UNSELECTING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      _3Dwin->selection_rect.x1 = x;
      _3Dwin->selection_rect.y1 = y;
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
        xge_OrderSelectionRect ( &_3Dwin->selection_rect );
        xge_callback ( er, xgemsg_3DWIN_UNSELECT_POINTS, 0, x, y );
        goto exit_select_mode;
      }
      xge_SetClipping ( er );
      er->redraw ( er, true );
      break;
  case xgemsg_MCLICK:
      if ( !(key & xgemouse_RBUTTON_DOWN) ) {
        _3Dwin->selection_rect.x1 = x;
        _3Dwin->selection_rect.y1 = y;
        xge_OrderSelectionRect ( &_3Dwin->selection_rect );
        xge_callback ( er, xgemsg_3DWIN_UNSELECT_POINTS, 0, x, y );
exit_select_mode:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( _3Dwin->fww.er );
        xge_SetClipping ( _3Dwin->fww.er );
        _3Dwin->fww.er->redraw ( _3Dwin->fww.er, true );
      }
      break;
  default:
      break;
    }
    break;

case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_MMOVE:
      break;
  case xgemsg_MCLICK:
      if ( _3Dwin->panning ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_3DWIN_PANNING;
          xge_GrabFocus ( _3Dwin->fww.er, true );
          xge_GrabFocus ( er, true );
          _3Dwin->xx = x;
          _3Dwin->yy = y;
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) && 
                  (key & xgemouse_RBUTTON_CHANGE) )
          goto enter_zooming_mode;
      }
      else if ( _3Dwin->selecting_mode ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_3DWIN_SELECTING;
          goto continue_selecting_mode;
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) && 
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          er->state = xgestate_3DWIN_UNSELECTING;
continue_selecting_mode:
          xge_GrabFocus ( _3Dwin->fww.er, true );
          xge_GrabFocus ( er, true );
          _3Dwin->xx = _3Dwin->selection_rect.x0 = _3Dwin->selection_rect.x1 = x;
          _3Dwin->yy = _3Dwin->selection_rect.y0 = _3Dwin->selection_rect.y1 = y;
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
      }
      else if ( _3Dwin->special_selecting_mode ) {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          if ( xge_callback ( er, xgemsg_3DWIN_SPECIAL_SELECT, id, x, y ) ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
                  (key & xgemouse_RBUTTON_CHANGE) ) {
          if ( xge_callback ( er, xgemsg_3DWIN_SPECIAL_UNSELECT, id, x, y ) ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
      }
      else {
        if ( (key & xgemouse_LBUTTON_DOWN) &&
             (key & xgemouse_LBUTTON_CHANGE) ) {
          er->state = xgestate_3DWIN_TURNING_VIEWER;
          xge_GrabFocus ( _3Dwin->fww.er, true );
          xge_GrabFocus ( er, true );
          _3Dwin->xx = x;
          _3Dwin->yy = y;
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
        else if ( (key & xgemouse_RBUTTON_DOWN) &&
                  (key & xgemouse_RBUTTON_CHANGE) ) {
enter_zooming_mode:
          er->state = xgestate_3DWIN_ZOOMING;
          xge_GrabFocus ( _3Dwin->fww.er, true );
          xge_GrabFocus ( er, true );
          _3Dwin->xx = _3Dwin->yy = y;
          _3Dwin->perspzoom = CPos->vd.persp.f;
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
      }
      break;
  case xgemsg_KEY:
      switch ( key ) {
    case 'C':
        xge_3DwinfSavePerspCamera ( _3Dwin );
        break;
    case 'c':
        xge_3DwinfSwapPerspCameras ( _3Dwin );
        xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, key, x, y );
        xge_SetClipping ( er );
        er->redraw ( er, true );
        break;
    case 'r':  case 'R':  /* reset projection */
        _3Dwin->PerspBBox = _3Dwin->DefBBox;
        xge_3DwinfSetupPerspProj ( _3Dwin, true );
        xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, key, x, y );
        xge_SetClipping ( er );
        er->redraw ( er, true );
        break;
    case 'f':  case 'F':  /* find bounding box */
        _3Dwin->PerspBBox = _3Dwin->RefBBox;
        if ( xge_callback ( er, xgemsg_3DWIN_FIND_REFBBOX, key, x, y ) ) {
          pkv_Exchange ( &_3Dwin->RefBBox, &_3Dwin->PerspBBox, sizeof(Box3f) );
          xge_3DwinfSetupPerspProj ( _3Dwin, false );
          xge_callback ( er, xgemsg_3DWIN_PROJCHANGE, key, x, y );
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
        break;
    case 'u':  case 'U':  /* undo, if possible */
        if ( xge_callback ( er, xgemsg_3DWIN_UNDO, id, x, y ) ) {
          xge_SetClipping ( _3Dwin->fww.er );
          _3Dwin->fww.er->redraw ( _3Dwin->fww.er, true );
        }
        break;
    default:  /* any other key to be consumed by the application */
        return xge_callback ( er, xgemsg_3DWIN_KEY, key, x, y );
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

exit_mode:
  er->state = xgestate_NOTHING;
  xge_ReleaseFocus ( _3Dwin->fww.er );
  xge_SetClipping ( er );
  er->redraw ( er, true );
  return true;
} /*xge_3DwinfPerspMsg*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_3DwinfSetDefBBox ( xge_3Dwinf *_3Dwin, float x0, float x1,
                            float y0, float y1, float z0, float z1 )
{
  _3Dwin->DefBBox.x0 = x0;  _3Dwin->DefBBox.x1 = x1;
  _3Dwin->DefBBox.y0 = y0;  _3Dwin->DefBBox.y1 = y1;
  _3Dwin->DefBBox.z0 = z0;  _3Dwin->DefBBox.z1 = z1;
} /*xge_3DwinfSetDefBBox*/

void xge_3DwinfSetupParProj ( xge_3Dwinf *_3Dwin, Box3f *bbox )
{
  point3f centre, p, q, r;
  float   dx, dy, dz, bx, by, bz, s;
  int     i;

  bx = bbox->x1-bbox->x0;
  by = bbox->y1-bbox->y0;
  bz = bbox->z1-bbox->z0;
  SetPoint3f ( &centre, (float)(0.5*(bbox->x0+bbox->x1)),
                        (float)(0.5*(bbox->y0+bbox->y1)),
                        (float)(0.5*(bbox->z0+bbox->z1)) );
  switch ( _3Dwin->fww.zoomwin ) {
case 0:
    dx = (float)_3Dwin->cwin[0]->w;
    dy = (float)_3Dwin->cwin[0]->h/xge_aspect;
    s = min ( (dx/bx), (dy/bz) );
    CameraInitFramef ( &_3Dwin->CPos[0], true, false,
                       _3Dwin->cwin[0]->w, _3Dwin->cwin[0]->h,
                       _3Dwin->cwin[0]->x, _3Dwin->cwin[0]->y, xge_aspect, 0 );
    _3Dwin->CPos[0].position = centre;
    _3Dwin->CPos[0].vd.para.wdt = (float)_3Dwin->cwin[0]->w/s;
    _3Dwin->CPos[0].vd.para.dim_case = 1;
    CameraTurnGf ( &_3Dwin->CPos[0], 0.0, (float)(-0.5*PI), 0.0 );
    break;

case 1:
    dx = (float)_3Dwin->cwin[1]->w;
    dy = (float)_3Dwin->cwin[1]->h/xge_aspect;
    s = min ( (dx/by), (dy/bz) );
    CameraInitFramef ( &_3Dwin->CPos[1], true, false,
                       _3Dwin->cwin[1]->w, _3Dwin->cwin[1]->h,
                       _3Dwin->cwin[1]->x, _3Dwin->cwin[1]->y, xge_aspect, 0 );
    _3Dwin->CPos[1].position = centre;
    _3Dwin->CPos[1].vd.para.wdt = (float)_3Dwin->cwin[1]->w/s;
    _3Dwin->CPos[1].vd.para.dim_case = 1;
    CameraTurnGf ( &_3Dwin->CPos[1], 0.0, (float)(-0.5*PI), (float)(-0.5*PI) );
    break;

case 2:
    dx = (float)_3Dwin->cwin[2]->w;
    dy = (float)_3Dwin->cwin[2]->h/xge_aspect;
    s = min ( (dx/bx), (dy/by) );
    CameraInitFramef ( &_3Dwin->CPos[2], true, false,
                       _3Dwin->cwin[2]->w, _3Dwin->cwin[2]->h,
                       _3Dwin->cwin[2]->x, _3Dwin->cwin[2]->y, xge_aspect, 0 );
    _3Dwin->CPos[2].position = centre;
    _3Dwin->CPos[2].vd.para.wdt = (float)_3Dwin->cwin[2]->w/s;
    _3Dwin->CPos[2].vd.para.dim_case = 1;
    CameraTurnGf ( &_3Dwin->CPos[2], 0.0, PI, 0.0 );
    break;

default:
    dx = (float)min ( _3Dwin->cwin[0]->w, _3Dwin->cwin[2]->w );
    dy = min ( _3Dwin->cwin[1]->w, _3Dwin->cwin[2]->h/xge_aspect );
    dz = min ( _3Dwin->cwin[0]->h, _3Dwin->cwin[1]->h )/xge_aspect;
    s = min ( (dx/bx), (dy/by) );
    s = min ( s, (dz/bz) );
    for ( i = 0; i < 3; i++ ) {
      CameraInitFramef ( &_3Dwin->CPos[i], true, false,
                         _3Dwin->cwin[i]->w, _3Dwin->cwin[i]->h,
                         _3Dwin->cwin[i]->x, _3Dwin->cwin[i]->y, xge_aspect, 0 );
      _3Dwin->CPos[i].position = centre;
      _3Dwin->CPos[i].vd.para.wdt = _3Dwin->cwin[i]->w/s;
      _3Dwin->CPos[i].vd.para.dim_case = 1;
      switch ( i ) {
    case 0:
        CameraTurnGf ( &_3Dwin->CPos[0], 0.0, (float)(-0.5*PI), 0.0 );
        break;
    case 1:
        CameraTurnGf ( &_3Dwin->CPos[1], 0.0, (float)(-0.5*PI), (float)(-0.5*PI) );
        break;
    case 2:
        CameraTurnGf ( &_3Dwin->CPos[2], 0.0, PI, 0.0 );
        break;
      }
    }
    break;
  }
        /* find the largest box fitting in the parallel windows */
  CameraProjectPoint3f ( &_3Dwin->CPos[0], &centre, &p );
  SetPoint3f ( &q, _3Dwin->CPos[0].xmin, _3Dwin->CPos[0].ymin, p.z );
  CameraUnProjectPoint3f ( &_3Dwin->CPos[0], &q, &r );
  _3Dwin->WinBBox.x0 = r.x;
  _3Dwin->WinBBox.z1 = r.z;
  SetPoint3f ( &q, _3Dwin->CPos[0].xmin+_3Dwin->CPos[0].width-1,
                   _3Dwin->CPos[0].ymin+_3Dwin->CPos[0].height-1, p.z );
  CameraUnProjectPoint3f ( &_3Dwin->CPos[0], &q, &r );
  _3Dwin->WinBBox.x1 = r.x;
  _3Dwin->WinBBox.z0 = r.z; 
  CameraProjectPoint3f ( &_3Dwin->CPos[1], &centre, &p );
  SetPoint3f ( &q, _3Dwin->CPos[1].xmin, _3Dwin->CPos[1].ymin, p.z );
  CameraUnProjectPoint3f ( &_3Dwin->CPos[1], &q, &r );
  _3Dwin->WinBBox.y1 = r.y;
  _3Dwin->WinBBox.z1 = min ( _3Dwin->WinBBox.z1, r.z );
  SetPoint3f ( &q, _3Dwin->CPos[1].xmin+_3Dwin->CPos[1].width-1,
                   _3Dwin->CPos[1].ymin+_3Dwin->CPos[1].height-1, p.z );
  CameraUnProjectPoint3f ( &_3Dwin->CPos[1], &q, &r );
  _3Dwin->WinBBox.y0 = r.y;
  _3Dwin->WinBBox.z0 = max ( _3Dwin->WinBBox.z0, r.z );
  CameraProjectPoint3f ( &_3Dwin->CPos[2], &centre, &p );
  SetPoint3f ( &q, _3Dwin->CPos[2].xmin, _3Dwin->CPos[2].ymin, p.z );
  CameraUnProjectPoint3f ( &_3Dwin->CPos[2], &q, &r );
  _3Dwin->WinBBox.x0 = max ( _3Dwin->WinBBox.x0, r.x );
  _3Dwin->WinBBox.y1 = min ( _3Dwin->WinBBox.y1, r.y );
  SetPoint3f ( &q, _3Dwin->CPos[2].xmin+_3Dwin->CPos[2].width-1,
                   _3Dwin->CPos[2].ymin+_3Dwin->CPos[2].height-1, p.z );
  CameraUnProjectPoint3f ( &_3Dwin->CPos[2], &q, &r );
  _3Dwin->WinBBox.x1 = min ( _3Dwin->WinBBox.x1, r.x );
  _3Dwin->WinBBox.y0 = max ( _3Dwin->WinBBox.y0, r.y );
} /*xge_3DwinfSetupParProj*/

void xge_3DwinfSetupPerspProj ( xge_3Dwinf *_3Dwin, boolean resetpos )
{
  CameraRecf *CPos;
  float      bx, by, bz, r;
  point3f    p;
  vector3f   v;

    /* the initial parameters of the camera are to be reviewed */
  CPos = &_3Dwin->CPos[3];
  CameraInitFramef ( CPos, false, false,
                     _3Dwin->cwin[3]->w, _3Dwin->cwin[3]->h,
                     _3Dwin->cwin[3]->x, _3Dwin->cwin[3]->y, xge_aspect, 0 );
  SetPoint3f ( &p, (float)(0.5*(_3Dwin->PerspBBox.x0+_3Dwin->PerspBBox.x1)),
                   (float)(0.5*(_3Dwin->PerspBBox.y0+_3Dwin->PerspBBox.y1)),
                   (float)(0.5*(_3Dwin->PerspBBox.z0+_3Dwin->PerspBBox.z1)) );
  bx = _3Dwin->PerspBBox.x1-_3Dwin->PerspBBox.x0;
  by = _3Dwin->PerspBBox.y1-_3Dwin->PerspBBox.y0;
  bz = _3Dwin->PerspBBox.z1-_3Dwin->PerspBBox.z0;
  r = (float)sqrt ( bx*bx+by*by+bz*bz );
  if ( resetpos ) {
    CameraInitPosf ( CPos );
    CPos->position = CPos->c_centre = CPos->g_centre = p;
    CPos->c_fixed = false;
    SetVector3f ( &v, 0.0, 0.0, (float)(-10.0*r) );
    CameraMoveGf ( CPos, &v );
    CameraRotXCf ( CPos, PI );
  }
  else {
    CPos->position = CPos->c_centre = CPos->g_centre = p;
    CPos->c_fixed = false;
    SetVector3f ( &v, 0.0, 0.0, (float)(-10.0*r) );
    CameraMoveCf ( CPos, &v );
  }
  CPos->c_fixed = true;
  CameraSetDepthRangef ( CPos, 8.0*r, 12.0*r );
  CameraSetFf ( CPos, 5.0 );
} /*xge_3DwinfSetupPerspProj*/

void xge_3DwinfUpdatePerspProj ( xge_3Dwinf *_3Dwin )
{
  CameraInitFramef ( &_3Dwin->CPos[3], false, false,
                     _3Dwin->cwin[3]->w, _3Dwin->cwin[3]->h,
                     _3Dwin->cwin[3]->x, _3Dwin->cwin[3]->y, xge_aspect, 0 );
  CameraSetMappingf ( &_3Dwin->CPos[3] );
} /*xge_3DwinfUpdatePerspProj*/

void xge_3DwinfInitProjections ( xge_3Dwinf *_3Dwin,
                   float x0, float x1, float y0, float y1, float z0, float z1 )
{
  int i;

  _3Dwin->RefBBox.x0 = x0;  _3Dwin->RefBBox.x1 = x1;
  _3Dwin->RefBBox.y0 = y0;  _3Dwin->RefBBox.y1 = y1;
  _3Dwin->RefBBox.z0 = z0;  _3Dwin->RefBBox.z1 = z1;
  _3Dwin->PerspBBox = _3Dwin->RefBBox;
  _3Dwin->perspzoom = 1.0;
  for ( i = 0; i < 3; i++ )
    CameraSetDepthRangef ( &_3Dwin->CPos[i], -10.0, 10.0 );
  CameraSetDepthRangef ( &_3Dwin->CPos[3], 1.0, 100.0 );
  xge_3DwinfSetupParProj ( _3Dwin, &_3Dwin->RefBBox );
  xge_3DwinfSetupPerspProj ( _3Dwin, true );
  _3Dwin->CPos[4] = _3Dwin->CPos[3];
} /*xge_3DwinfInitProjections*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *xge_New3Dwinf ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            xge_3Dwinf *_3Dwin,
                            void (*pararedraw)(xge_widget*, boolean),
                            void (*perspredraw)(xge_widget*, boolean) )
{
  int        i;
  xge_widget *fww;

  _3Dwin->cwin[3] = xge_NewWidget ( window_num, NULL, id+3, 100, 100, 0, 0,
                    _3Dwin, &_3Dwin->CPos[3], xge_3DwinfPerspMsg, perspredraw );
  for ( i = 2; i >= 0; i-- )
    _3Dwin->cwin[i] = xge_NewWidget ( window_num, _3Dwin->cwin[i+1], id+i,
                      100, 100, 0, 0,
                      _3Dwin, &_3Dwin->CPos[i], xge_3DwinfParaMsg, pararedraw );

  fww = xge_NewFourWW ( window_num, prev, id,
                        w, h, x, y, _3Dwin->cwin[0], &_3Dwin->fww );
  if ( fww ) {
    _3Dwin->panning = _3Dwin->selecting_mode = false;
    _3Dwin->display_coord = false;
    _3Dwin->current_tool = xge_3DWIN_NO_TOOL;
    _3Dwin->moving_tool = _3Dwin->scaling_tool = _3Dwin->rotating_tool = false;
    _3Dwin->CoordWin = -1;
    xge_3DwinfSetDefBBox ( _3Dwin, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0 );
    xge_3DwinfInitProjections ( _3Dwin, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0 );
    xge_3DwinfResetGeomWidgets ( _3Dwin );
  }
  return fww;
} /*xge_New3Dwinf*/

