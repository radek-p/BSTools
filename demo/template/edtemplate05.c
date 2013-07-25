
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

boolean FilenameCorrect ( const char *fn )
{
  return (boolean)(fn[0] != 0);
} /*FilenameCorrect*/

void OpenFile ( const char *fn )
{
  xge_DisplayErrorMessage ( ErrMsgNoAction, -1 );
} /*OpenFile*/

void SaveFile ( const char *fn )
{
  if ( !FilenameCorrect ( fn ) ) {
    xge_DisplayErrorMessage ( ErrMsgBadFilename, -1 );
    return;
  }
  xge_DisplayErrorMessage ( ErrMsgNoAction, -1 );
} /*SaveFile*/

