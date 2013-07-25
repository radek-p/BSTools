
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

char txtNULL[]        = "";
char txtFile[]        = "File";
char txtAbout[]       = "About";
char txtMark_unmark[] = "Mark/Unmark";
char txtMove[]        = "Move";
char txtScale[]       = "Scale";
char txtRotate[]      = "Rotate";
char txtPan_zoom[]    = "Pan/Zoom";
char txtCoordinates[] = "Coordinates";
char txtNew[]         = "New";
char txtOpen[]        = "Open";
char txtSave[]        = "Save";
char txtSave_as[]     = "Save as";
char txtCancel[]      = "Cancel";
char txtExit[]        = "Exit";
char txtReset[]       = "Reset";

char ErrMsgNoAction[] = "Error: no action is implemented for now.";
char ErrMsgBadFilename[] = "Error: incorrect filename.";

char *InfoMsg[6] =
  { "Template for a two-window geometry editor program",
    "using the libxgedit widgets",
    "*** an unfinished version ***",
    "This program is a part of the BSTools package.",
    "(C) Copyright by Przemyslaw Kiciak, 2008",
    NULL };
