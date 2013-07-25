
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
#include "xgedit.h"

#include "spltemplate.h"
#include "edtempwidgets.h"
#include "edtemplate.h"

void RysujDomOkno ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *dw;

  dw = er->data0;
  xge_DrawGeomWinBackground ( er );
  if ( dw->inside && dw->display_coord )
    xge_2DwindDrawCursorPos ( dw, xge_xx, xge_yy );

  xge_DrawGeomWinSelectionRect ( er, &dw->selection_rect );
  xge_2DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujDomOkno*/

void RysujKnotOkno ( xge_widget *er, boolean onscreen )
{
  xge_DrawGeomWinBackground ( er );

  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujKnotOkno*/

boolean KnotOknoMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_RESIZE:
    er->w = x;
    er->h = y;
    if ( key ) {
      xge_SetClipping ( er );
      er->redraw ( er, true );
    }
    return 1;

default:
    break;
  }
  return false;
} /*KnotOknoMsg*/

