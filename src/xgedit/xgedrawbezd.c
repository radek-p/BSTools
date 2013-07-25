
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2009                            */
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

#include "pkgeom.h"
#include "multibs.h"

#include "xgedit.h"
#include "xgeprivate.h"


void xge_DrawBC2d ( int n, const point2d *cp )
{
  mbs_RasterizeBC2d ( n, cp, xge_OutPixels, true );
} /*xge_DrawBC2d*/

void xge_DrawBC2Rd ( int n, const point3d *cp )
{
  mbs_RasterizeBC2Rd ( n, cp, xge_OutPixels, true );
} /*xge_DrawBC2Rd*/

