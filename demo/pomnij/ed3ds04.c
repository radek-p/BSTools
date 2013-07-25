
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
void RysujDOkno ( xge_widget *er, boolean onscreen )
{
  xge_T2KnotWind *T2win;

  T2win = er->data0;
  xge_DrawGeomWinBackground ( er );

  if ( T2win->display_coord && T2win->inside )
    xge_T2KnotWindDrawCursorPos ( T2win, xge_xx, xge_yy );
  DisplayDomain ();
  if ( display_domain_net )
    DisplayDomainNet ();
  xge_DrawGeomWinSelectionRect ( er, &T2win->selection_rect );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujDOkno*/

