
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

void RysujDomOkno ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *dw;

  dw = er->data0;
  xge_DrawGeomWinBackground ( er );
  if ( dw->inside && dw->display_coord )
    xge_2DwindDrawCursorPos ( dw, xge_xx, xge_yy );

  if ( view_dom_spatches )
    DrawDomainSurrPatches ();
  if ( view_dom_patches1 )
    DrawDomainPatches ( 1 );
  if ( view_dom_patches2 )
    DrawDomainPatches ( 2 );
  if ( view_dom_cp )
    DrawGHDomainControlNet ();
  if ( view_dom_numbers )
    DrawDomainNumbers ();

  xge_DrawGeomWinSelectionRect ( er, &dw->selection_rect );
  xge_2DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujDomOkno*/

