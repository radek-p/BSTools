
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

void RysujOkno ( xge_widget *er, boolean onscreen )
{
  int id;
  xge_3Dwind *sw;

  sw = er->data0;
  id = er->id & 0x03;  /* ought to be 0 .. 2 */
  xge_DrawGeomWinBackground ( er );
  xge_3DwindDrawCursorPos ( sw, id, xge_xx, xge_yy );

  if ( view_surf_spatches )
    DrawGHSurfPatches ( id );
  if ( view_surf_1 && options1.surf_valid )
    DrawGHSurfFillingPatches ( id, 1 );
  if ( view_surf_2 && options2.surf_valid )
    DrawGHSurfFillingPatches ( id, 2 );
  switch ( swind_ed_switch ) {
  case SWIN_SELECT_VIEW:
      if ( view_surf_cp )
        DrawGHControlNet ( id );
      if ( view_constraints_frame ) {
        if ( constraints1 )
          DrawConstraintFrame ( id, 1 );
        else if ( constraints2 )
          DrawConstraintFrame ( id, 2 );
      }
      break;
case SWIN_EDITING_SURFACE:
    DrawGHControlNet ( id );
    break;
case SWIN_EDITING_CONSTRAINTS:
    if ( view_surf_cp )
      DrawGHControlNet ( id );
    if ( constraints1 )
      DrawConstraintFrame ( id, 1 );
    else if ( constraints2 )
      DrawConstraintFrame ( id, 2 );
    break;
case SWIN_EDITING_LIGHT:
    if ( view_surf_cp )
      DrawGHControlNet ( id );
    DrawLightPoints ( id );
    break;
default:
    break;
  }
  if ( view_surf_numbers )
    DrawGHSurfNumbers  ( id );

  xge_DrawGeomWinSelectionRect ( er, &sw->selection_rect );
  xge_3DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujOkno*/

void RysujPOkno ( xge_widget *er, boolean onscreen )
{
  int id;
  xge_3Dwind *sw;

  sw = er->data0;
  id = er->id & 0x03;  /* ought to be 3 */
  if ( swind_picture ) {
    XPutImage ( xgedisplay, xgepixmap, xgegc, rendimage,
                er->x, er->y, er->x, er->y, er->w, er->h );
    if ( onscreen )
      xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
  }
  else {
    xge_DrawGeomWinBackground ( er );

    if ( view_surf_spatches )
      DrawGHSurfPatches ( id );
    if ( view_surf_1 && options1.surf_valid )
      DrawGHSurfFillingPatches ( id, 1 );
    if ( view_surf_2 && options2.surf_valid )
      DrawGHSurfFillingPatches ( id, 2 );
    switch ( swind_ed_switch ) {
  case SWIN_SELECT_VIEW:
      if ( view_surf_cp )
        DrawGHControlNet ( id );
      if ( view_constraints_frame ) {
        if ( constraints1 )
          DrawConstraintFrame ( id, 1 );
        else if ( constraints2 )
          DrawConstraintFrame ( id, 2 );
      }
      break;
  case SWIN_EDITING_SURFACE:
      DrawGHControlNet ( id );
      break;
  case SWIN_EDITING_CONSTRAINTS:
      if ( view_surf_cp )
        DrawGHControlNet ( id );
      if ( constraints1 )
        DrawConstraintFrame ( id, 1 );
      else if ( constraints2 )
        DrawConstraintFrame ( id, 2 );
      break;
  case SWIN_EDITING_LIGHT:
      if ( view_surf_cp )
        DrawGHControlNet ( id );
      DrawLightPoints ( id );
      break;
  default:
      break;
    }
    if ( view_surf_numbers )
      DrawGHSurfNumbers  ( id );

    xge_DrawGeomWinSelectionRect ( er, &sw->selection_rect );
  }
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujPOkno*/

