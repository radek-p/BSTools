
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>   
#include <stdio.h>
#include <math.h>  
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>  
#include <GL/glu.h> 
#include <GL/glx.h>

#include "pkvaria.h" 
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "arrowbtn.h"

void DrawButtonArrowUp ( xge_widget *er, boolean onscreen )
{
  short x, y1, y2;

  xge_DrawButton ( er, false );
  x = er->x + (er->w-1)/2;
  y1 = er->y + er->h/2;
  y2 = y1 + 5;
  y1 -= 5;
  xgeDrawLine ( x, y1, x, y2 );
  xgeDrawLine ( x, y1, x-3, y1+3 );
  xgeDrawLine ( x, y1, x+3, y1+3 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*DrawButtonArrowUp*/

void DrawButtonArrowDown ( xge_widget *er, boolean onscreen )
{
  short x, y1, y2;

  xge_DrawButton ( er, false );
  x = er->x + (er->w-1)/2;
  y1 = er->y + er->h/2;
  y2 = y1 + 5;
  y1 -= 5;
  xgeDrawLine ( x, y1, x, y2 );
  xgeDrawLine ( x, y2, x-3, y2-3 );
  xgeDrawLine ( x, y2, x+3, y2-3 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*DrawButtonArrowDown*/

