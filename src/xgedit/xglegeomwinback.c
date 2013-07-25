
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2010                            */
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

#include <GL/gl.h>
#include <GL/glx.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"

#include "xgedit.h"
#include "xgledit.h"
#include "xgeprivate.h"


/* ///////////////////////////////////////////////////////////////////////// */
void xgle_DrawGeomWinBackground ( xge_widget *er, GLbitfield mask )
{
  glDisable ( GL_LIGHTING );
  switch ( er->state ) {
case xgestate_2DWIN_MOVINGPOINT:
case xgestate_2DWIN_ZOOMING:
case xgestate_2DWIN_PANNING:
case xgestate_2DWIN_SELECTING:
case xgestate_2DWIN_UNSELECTING:
case xgestate_2DWIN_MOVING_GEOM_WIDGET:
case xgestate_2DWIN_USING_GEOM_WIDGET:
case xgestate_2DWIN_ALTUSING_GEOM_WIDGET:
case xgestate_3DWIN_MOVINGPOINT:
case xgestate_3DWIN_ZOOMING:
case xgestate_3DWIN_PARPANNING:
case xgestate_3DWIN_PARZOOMING:
case xgestate_3DWIN_TURNING_VIEWER:
case xgestate_3DWIN_SELECTING:
case xgestate_3DWIN_UNSELECTING:
case xgestate_3DWIN_MOVING_GEOM_WIDGET:
case xgestate_3DWIN_USING_GEOM_WIDGET:
case xgestate_3DWIN_ALTUSING_GEOM_WIDGET:
case xgestate_KNOTWIN_MOVINGKNOT:
case xgestate_KNOTWIN_PANNING:
case xgestate_KNOTWIN_ZOOMING:
case xgestate_T2KNOTWIN_MOVINGKNOT_U:
case xgestate_T2KNOTWIN_MOVINGKNOT_V:
case xgestate_T2KNOTWIN_MOVING_POINT:
case xgestate_T2KNOTWIN_PANNING:
case xgestate_T2KNOTWIN_ZOOMING:
case xgestate_T2KNOTWIN_SELECTING:
case xgestate_T2KNOTWIN_UNSELECTING:
    xgleClearColor3fv ( xglec_Blue7 );
    break;
default:
    xgleClearColor3fv ( xglec_Blue6 );
    break;
  }
  glClear ( GL_COLOR_BUFFER_BIT | mask );
} /*xgle_DrawGeomWinBackground*/

