
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2009                            */
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

void xge_RepositionWidgets ( short w, short h, short x, short y,
                             xge_widget *edr )
{
/* depending on its rpos, position each widget relative to the */
/* 0 - upper left, 1 - upper right, 2 - lower left, 3 - lower right corner */
  while ( edr ) {
    if ( edr->rpos >= 0 ) {
      edr->x = (short)(x + edr->xofs);
      if ( edr->rpos & 0x01 )
        edr->x = (short)(edr->x+w);
      edr->y = (short)(y + edr->yofs);
      if ( edr->rpos & 0x02 )
        edr->y = (short)(edr->y+h);
    }
    edr = edr->prev;
  }
} /*xge_RepositionWidgets*/

