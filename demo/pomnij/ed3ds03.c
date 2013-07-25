
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
void RysujPOkno ( xge_widget *er, boolean onscreen )
{
  int id;

  id = er->id & 0x03;  /* ought to be 3 */
  if ( swind_picture ) {
    XPutImage ( xgedisplay, xgepixmap, xgegc, rendimage,
                er->x, er->y, er->x, er->y, er->w, er->h );
    if ( onscreen )
      xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
  }
  else {
    xge_DrawGeomWinBackground ( er );
    if ( display_surface )
      DisplaySurface ( id );
    if ( display_Bezier_nets )
      DisplayBezierNets ( id );
    if ( sw_blending_constraints )
      DisplayConstraintCurves ( id );
    if ( display_control_net ) {
      DisplayControlNet ( id );
      DisplayControlPoints ( id );
    }
    if ( win1_contents == WIN1_BLENDING ) {
      if ( display_constr_poly )
        DisplayG2BlendingConstrCP ( id );
      if ( display_trans_net )
        DisplayPreTransControlNet ( id );
    }
    xge_DrawGeomWinSelectionRect ( er, &swind.selection_rect );
    xge_3DwindDrawGeomWidgets ( er );
    xge_DrawGeomWinFrame ( er, onscreen );
  }
} /*RysujPOkno*/

