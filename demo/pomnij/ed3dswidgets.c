
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "xgedit.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"
#include "ed3dswidgets.h"

char txtNULL[]              = "";
char txtReset[]             = "Reset";
char txtInit[]              = "Init";
char txtEdit[]              = "Edit";
char txtView[]              = "View";
char txtPicture[]           = "Picture";
char txtAbout[]             = "About";
char txtDegree[]            = "degree";
char txtDegree_u[]          = "degree u";
char txtDegree_v[]          = "degree v";
char txtEquidist_u[]        = "equidistant u";
char txtEquidist_v[]        = "equidistant v";
char txtFlip[]              = "Flip";
char txtPan_zoom[]          = "pan & zoom";
char txtCoordinates[]       = "coordinates";
char txtControl_net[]       = "control net ";
char txtSurface[]           = "surface";
char txt_Surface[]          = "Surface";
char txt_Light[]            = "Light";
char txtMove[]              = "move";
char txtScale[]             = "scale";
char txtRotate[]            = "rotate";
char txtShear[]             = "shear";
char txtU_density[]         = "u density";
char txtV_density[]         = "v density";
char txtMark_unmark[]       = "mark/unmark";
char txtBezier_nets[]       = "Bezier nets";
char txtConstrPoly[]        = "constraint p.";
char txtTransfNet[]         = "transformed net";
char txtMove_many_knots[]   = "move many knots";
char txtClosed[]            = "closed";
char txtClosed_u[]          = "closed u";
char txtClosed_v[]          = "closed v";
char txtDomain_net[]        = "domain net";
char txtFile[]              = "File";
char txtOpen[]              = "Open";
char txtSave[]              = "Save";
char txtSave_as[]           = "Save as";
char txtExport[]            = "Export";
char txtExit[]              = "Exit";
char txtOK[]                = "OK";
char txtCancel[]            = "Cancel";
char txtSurface_Type[]      = "Surface type";
char txtGeneral[]           = "General";
char txtGeneral_BSpline[]   = "General B-spline";
char txtSpherical[]         = "Spherical";
char txtSpherical_product[] = "Spherical product";
char txtSwept[]             = "Swept";
char txtBlending[]          = "Blending";
char txtEquator[]           = "Equator";
char txtMeridian[]          = "Meridian";
char txtBind[]              = "Bind";
char txtNURBS[]             = "NURBS";
char txtCurve[]             = "curve";
char txtControl_polygon[]   = "control polygon";
char txtBezier_polygons[]   = "Bezier polygons";
char txtTicks[]             = "ticks";
char txtRender[]            = "Render";
char txtC[]                 = "c.";
char txtD_shape_func[]      = "d. shape func.";
char txtGaussian[]          = "Gaussian";
char txtMean[]              = "mean";
char txtIsophotes[]         = "isophotes";
char txtReflection[]        = "reflection";
char txtHighlight[]         = "highlight";
char txtSections[]          = "sections";
char txtParamQual[]         = "param. qual.";
char txtShadows[]           = "shadows";
char txtAntialias[]         = "antialias";
char txtLight_0[]           = "light 0";
char txtLight_1[]           = "light 1";
char txtLight_2[]           = "light 2";
char txtLight_3[]           = "light 3";
char txtAmbient[]           = "ambient";
char txtReflectionFrame[]   = "refl. frame";
char txtHighlightFrame[]    = "highlight frame";
char txtSectionsFrame[]     = "sections frame";
char txtArcs[]              = "Arcs";
char txtQuarter[]           = "Quarter";
char txtHalf[]              = "Half";
char txtFull_circle[]       = "Full circle";
char txtArc_angle[]         = "Arc angle";
char txtSwept_Surface[]     = "Swept surface";
char txtBlending_Surface[]  = "Blending surface";
char txtClamped[]           = "clamped";
char txtBiharmonic[]        = "biharmonic";
char txtTriharmonic[]       = "triharmonic";
char txtRefine[]            = "Refine";
char txtConstraints[]       = "constraints";
char txtTransform[]         = "Transform";
char txtOptions[]           = "Options";
char txtShowSteps[]         = "show steps";
char txtInfo[]              = "Info";
char txtOptimize[]          = "Optimize";
char txtInterrupt[]         = "Interrupt";
char txtIterationLimit[]    = "iteration limit";
char txtQuadratureKnots[]   = "quadrature knots";
char txtDumpData[]          = "dump data";
char txtIdentity[]          = "Identity";
char txtG1[]                = "G1";
char txtG2[]                = "G2";

/* information message */
char *InfoMsg[5] =
  { "B-spline & NURBS surface demonstration program",
    "*** an unfinished version ***",
    "This program is a part of the BSTools package.",
    "(C) Copyright by Przemyslaw Kiciak, 2007, 2010",
    NULL };

char MsgReconsider[] = "Please, reconsider it and then proceed.";

/* error messages */
char ErrorMsgCannotOpen[]          = "Error: This file cannot be read.";
char ErrorMsgCannotSave[]          = "Error: Cannot save anything yet.";
char ErrorMsgCannotRaiseDeg[]      = "Error: Cannot raise degree above the limit.";
char ErrorMsgCannotReduceDeg[]     = "Error: Cannot reduce degree below 1.";
char ErrorMsgCannotInsertKnot[]    = "Error: Cannot insert a knot at this point.";
char ErrorMsgTooManyKnots[]        = "Error: Too many knots.";
char ErrorMsgCannotRemoveKnot[]    = "Error: Cannot remove this knot.";
char ErrorMsgNotEnoughKnots[]      = "Error: Not enough knots.";
char ErrorMsgNotEnoughUKnots[]     = "Error: Not enough knots in the \"u\" sequence.";
char ErrorMsgNotEnoughVKnots[]     = "Error: Not enough knots in the \"v\" sequence.";
char ErrorMsgFileWritingError[]    = "Error: Could not write the file.";
char ErrorMsgRefineBlending[]      = "Error: Cannot refine the surface.";
char ErrorMsgBindBlending[]        = "Error: Blending surface must be bound.";
char ErrorMsgNonlinBlending[]      = "Error: Nonlinear optimization failed.";
char ErrorMsgCannotAddConstraint[] = "Error: Cannot add a constraint.";
char ErrorMsgIncorrectNumber[]     = "Error: Incorrect number.";
char ErrorMsgBiquadratic[]         = "Error: G1 blending patch must be biquadratic.";
char ErrorMsgBicubic[]             = "Error: G1 blending patch must be bicubic.";

char ErrorChildProcessNotActive[]  = "Error: Child process is not active.";
char ErrorChildProcessTerminated[] = "Error: Child process terminated.";

char ErrorMsgNotImplemented[]      =
    "Error: This part of the program is not implemented yet.";

char WarningMsgPovExport[]         =
    "Warning: PovRay accepts only polynomial patches of degree (3,3).";

