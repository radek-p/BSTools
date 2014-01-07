
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

#include "pkgeom.h"
#include "cameraf.h"
#include "multibs.h"
#include "xgedit.h"

#include "spline2d.h"
#include "ed2dwidgets.h"

char txtFile[]    = "File";
char txtOpen[]    = "Open";
char txtSave[]    = "Save";
char txtSaveAs[]  = "Save as";
char txtExport[]  = "Export";
char txtExit[]    = "Exit";
char txtCancel[]  = "Cancel";

char b1text[]     = "Edit";
char b2text[]     = "View";
char b3text[]     = "About";
char b4text[]     = "Refine";
char sw19text[]   = "uniform";
char b5text[]     = "Reset";
char sw0text[]    = "select/unselect";
char sw1text[]    = "move";
char sw2text[]    = "scale";
char sw3text[]    = "rotate";
char txtShear[]   = "shear";
char intw2text[]  = "degree";
char sw4text[]    = "curve";
char sw5text[]    = "control polygon";
char sw6text[]    = "Bezier polygons";
char sw7text[]    = "ticks";
char sw8text[]    = "convex hulls";
char sw9text[]    = "NURBS";
char sw10text[]   = "closed";
char sw11text[]   = "function";
char sw12text[]   = "basis functions";
char sw13text[]   = "diagonal forms";
char sw14text[]   = "curvature graph";
char intw15text[] = "graph density";
char sw16text[]   = "pan & zoom";
char sw17text[]   = "coordinates";
char sw18text[]   = "move many knots";

char *InfoMsg[5] =
  { "Planar B-spline & NURBS curve demonstration program",
    "",
    "This program is a part of the BSTools package.",
    "(C) Copyright by Przemyslaw Kiciak, 2005, 2014",
    NULL };

char MsgReconsider[] = "Please, reconsider it and then proceed.";

char ErrMsgRaiseDeg[]      = "Error: Cannot raise degree above the limit.";
char ErrMsgReduceDeg[]     = "Error: Cannot reduce degree below 1.";
char ErrMsgCannotInsert[]  = "Error: Cannot insert a knot at this point.";
char ErrMsgToManyKnots[]   = "Error: Too many knots.";
char ErrMsgCannotRemove[]  = "Error: Cannot remove this knot.";
char ErrMsgCannotClose[]   = "Error: Not enough knots to close the curve.";
char ErrMsgCannotSave[]    = "Error: Cannot write the file.";
char ErrMsgCannotRead[]    = "Error: Cannot read this file.";
char ErrMsgCannotRefine[]  = "Error: Can refine only curves with uniform knots.";
