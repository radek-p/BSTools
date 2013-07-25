
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
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "splhole.h"  

char txtNULL[]            = "";
char txtFile[]            = "File";
char txtData[]            = "Data";
char txtEdit[]            = "Edit";
char txtConstraints[]     = "Constraints";
char txtLight[]           = "Light";
char txtView[]            = "View"; 
char txtPicture[]         = "Picture";
char txtAbout[]           = "About";
char txtSides[]           = "Number of sides";
char txtMark_unmark[]     = "Mark/Unmark";
char txtMove[]            = "Move";
char txtScale[]           = "Scale";
char txtRotate[]          = "Rotate";
char txtShear[]           = "Shear";
char txtPan_zoom[]        = "Pan/Zoom";
char txtCoordinates[]     = "Coordinates";
char txtNew[]             = "New";
char txtOpen[]            = "Open";
char txtSave[]            = "Save";
char txtSave_as[]         = "Save as";
char txtCancel[]          = "Cancel";
char txtExit[]            = "Exit";
char txtReset[]           = "Reset";
char txtOptions[]         = "Options";
char txtInfo[]            = "Info";
char txtControl_net[]     = "control net";
char txtSurr_patches[]    = "surr. patches";
char txtDomain_patches1[] = "dom. patches 1";
char txtDomain_patches2[] = "dom. patches 2";
char txtNumbers[]         = "numbers";
char txtFirst[]           = "1st";
char txtSecond[]          = "2nd";
char txtConstrFrame[]     = "constr. frame";
char txtGCOrder[]         = "GC order";
char txtRestricted[]      = "Restricted";
char txtCoons[]           = "Coons";
char txtBezier[]          = "Bezier";
char txtSpline[]          = "Spline";
char txtNknots[]          = "number of knots";
char txtMult1[]           = "multiplicity 1";
char txtMult2[]           = "multiplicity 2";
char txtSurface[]         = "Surface";
char txt_surface[]        = "surface";
char txtLinear[]          = "linear";
char txtQuasiG2[]         = "quasi G2";
char txtAltCentre[]       = "alt. centre";
char txtGaussLegendre[]   = "Gauss-Legendre";
char txtRender[]          = "Render";
char txtStop[]            = "Stop";
char txtC[]               = "c.";
char txtD_shape_func[]    = "d. shape func.";
char txtGaussian[]        = "Gaussian";
char txtMean[]            = "mean";
char txtIsophotes[]       = "isophotes";
char txtReflection[]      = "reflection";
char txtHighlight[]       = "highlight";
char txtSections[]        = "sections";
char txtShadows[]         = "shadows";
char txtAntialias[]       = "antialias";
char txtLight_0[]         = "light 0";
char txtLight_1[]         = "light 1";
char txtLight_2[]         = "light 2";
char txtLight_3[]         = "light 3";
char txtAmbient[]         = "ambient";
char txtReflectionFrame[] = "refl. frame";
char txtHighlightFrame[]  = "highlight frame";
char txtSectionsFrame[]   = "sections frame";
char txtFirstType[]       = "1st type";
char txtSecondType[]      = "2nd type";
char txtZero[]            = "zero";
char txtCurrent[]         = "Current";
char txtCentralPoint[]    = "central point";

char *InfoMsg[6] =
  { "Demonstration program for filling holes in piecewise",
    "bicubic spline surface, with G^1 and G^2 continuity",
    "*** an almost finished version ***",
    "This program is a part of the BSTools package.",
    "(C) Copyright by Przemyslaw Kiciak, 2008, 2009",
    NULL };

char MsgReconsider[] = "Please, reconsider it and then proceed.";

char ErrMsgNoAction[]       = "Error: no action is implemented for now.";
char ErrMsgBadFilename[]    = "Error: incorrect filename.";
char ErrMagBadConstrFrame[] = "Error: cannot get the constraint frame.";
char ErrMsgFileReading[]    = "Error: cannot read this file.";
char ErrMsgFileWriting[]    = "Error: cannot write the file.";


