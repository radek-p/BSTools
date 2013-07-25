
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

void RysujOkno ( xge_widget *er, boolean onscreen )
{
  int id;
  xge_3Dwind *sw;

  sw = er->data0;
  id = er->id & 0x03;
  xge_DrawGeomWinBackground ( er );
  xge_3DwindDrawCursorPos ( sw, id, xge_xx, xge_yy );


  xge_DrawGeomWinSelectionRect ( er, &sw->selection_rect );
  xge_3DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujOkno*/

void RysujPOkno ( xge_widget *er, boolean onscreen )
{
  xge_3Dwind *sw;

  sw = er->data0;
  xge_DrawGeomWinBackground ( er );


  xge_DrawGeomWinSelectionRect ( er, &sw->selection_rect );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujPOkno*/

