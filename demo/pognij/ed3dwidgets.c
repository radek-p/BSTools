
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <string.h>

#include "ed3dwidgets.h"

char txtFile[]           = "File";
char txtOpen[]           = "Open";
char txtSave[]           = "Save";
char txtSaveAs[]         = "Save as";
char txtExit[]           = "Exit";
char txtCancel[]         = "Cancel";

char txtEdit[]           = "Edit";
char txtView[]           = "View";
char txtAbout[]          = "About";
char txtUniform[]        = "uniform";
char txtRefine[]         = "Refine";
char txtReset[]          = "Reset";

char txtDegree[]         = "degree";
char txtSelectUnselect[] = "select/unselect";
char txtMove[]           = "move";
char txtScale[]          = "scale";
char txtRotate[]         = "rotate";
char txtShear[]          = "shear";
char txtControlPolygon[] = "control polygon";
char txtCurve[]          = "curve";
char txtTicks[]          = "ticks";
char txtBezierPolygons[] = "Bezier polygons";
char txtConvexHulls[]    = "convex hulls";
char txtNURBS[]          = "NURBS";
char txtClosed[]         = "closed";
char txtDiagonalForms[]  = "diagonal forms";
char txtCurvatureGraph[] = "curvature graph";
char txtTorsionGraph[]   = "torsion graph";
char txtGraphDensity[]   = "graph density";
char txtPanZoom[]        = "pan & zoom";
char txtCoordinates[]    = "coordinates";
char txtMoveManyKnots[]  = "move many knots";

char *InfoMsg[5] =
  { "3D B-spline & NURBS curve demonstration program",
    "",
    "This program is a part of the BSTools package.",
    "(C) Copyright by Przemyslaw Kiciak, 2005, 2011",
    NULL };

char MsgReconsider[] = "Please, reconsider it and then proceed.";

char ErrMsgCannotSave[] = "";
char ErrMsgCannotRead[] = "";
char ErrMsgCannotRefine[] = "";

