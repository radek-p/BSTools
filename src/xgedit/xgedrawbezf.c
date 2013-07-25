
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


void xge_DrawBC2f ( int n, const point2f *cp )
{
  mbs_RasterizeBC2f ( n, cp, xge_OutPixels, true );
} /*xge_DrawBC2f*/

void xge_DrawBC2Rf ( int n, const point3f *cp )
{
  mbs_RasterizeBC2Rf ( n, cp, xge_OutPixels, true );
} /*xge_DrawBC2Rf*/

