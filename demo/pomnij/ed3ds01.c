
/* //////////////////////////////////////////////////// */ 
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
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
#include "xgedit.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"

/* ///////////////////////////////////////////////////////////////////////// */
void RysujSPROkno ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *_2Dwin;

  _2Dwin = er->data0;
  xge_DrawGeomWinBackground ( er );
  DisplayEqMerAxes ();
  if ( _2Dwin->inside && _2Dwin->display_coord )
    xge_2DwindDrawCursorPos ( _2Dwin, xge_xx, xge_yy );

  if ( display_eqmer_Bezier_polygons )
    DisplayEqMerBezierPolygons ();
  if ( display_eqmer_control_polygon )
    DisplayEqMerControlPolygon ();
  if ( display_eqmer_curve )
    DisplayEqMerCurve ();
  if ( display_eqmer_control_polygon )
    DisplayEqMerControlPoints ();

  xge_DrawGeomWinSelectionRect ( er, &_2Dwin->selection_rect );
  xge_2DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujSPROkno*/

