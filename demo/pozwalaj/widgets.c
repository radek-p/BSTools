
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>

#include "widgets.h"


char txtNull[]              = "";
char txtFile[]              = "File";
char txtObjects[]           = "Objects";
char txtEdit[]              = "Edit";
char txtTransform[]         = "Transform";
char txtPicture[]           = "Picture";
char txtAbout[]             = "About";
char txtBit[]               = "bit";
char txtMarkUnmark[]        = "mark/unmark";
char txtMark[]              = "mark";
char txtUnmark[]            = "unmark";
char txtTranslate[]         = "translate";
char txtScale[]             = "scale";
char txtRotate[]            = "rotate";
char txtShear[]             = "shear";
char txtPanZoom[]           = "pan/zoom";
char txtCoordinates[]       = "coordinates";
char txtOK[]                = "OK";
char txtCancel[]            = "Cancel";
char txtExit[]              = "Exit";
char txtOpen[]              = "Open";
char txtSave[]              = "Save";
char txtSaveAs[]            = "Save as";
char txtOptions[]           = "Options";
char txtView[]              = "View";
char txtData[]              = "Data";
char txtName[]              = "Name:";
char txtAdd[]               = "Add";
char txtReplace[]           = "Replace";
char txtNew[]               = "New";
char txtCopy[]              = "Copy";
char txtDelete[]            = "Delete";
char txtPurge[]             = "Purge";
char txtActive[]            = "Active";
char txtBezierCurve[]       = "Bezier curve";
char txtBSplineCurve[]      = "B-spline curve";
char txtBezierPatch[]       = "Bezier patch";
char txtBSplinePatch[]      = "B-spline patch";
char txtBSplineMesh[]       = "B-spline mesh";
char txtBSplineHole[]       = "B-spline hole";
char txtDegree[]            = "degree";
char txtDegreeU[]           = "degree u";
char txtDegreeV[]           = "degree v";
char txt2D[]                = "2D";
char txt3D[]                = "3D";
char txtRational[]          = "rational";
char txtRefine[]            = "refine";
char txtDouble[]            = "double";
char txtAverage[]           = "average";
char txtTetrahedron[]       = "tetrahedron";
char txtCube[]              = "cube";
char txtDodecahedron[]      = "dodecahedron";
char txt_gon[]              = "-gon";
char txt_gonalPrism[]       = "-gonal prism";
char txtControlPolygon[]    = "control polygon";
char txtControlNet[]        = "control net";
char txtCurve[]             = "curve";
char txtSurface[]           = "surface";
char txtHoleFilling[]       = "hole filling";
char txtVertex[]            = "vertex";
char txtLine[]              = "line";
char txtEdge[]              = "edge";
char txtFacet[]             = "facet";
char txtRemove[]            = "remove";
char txtFilter[]            = "filter";
char txtDoubleEdges[]       = "double edges";
char txtGlueLoops[]         = "glue loops";
char txtShrink[]            = "shrink";
char txtContract[]          = "contract";
char txtDensity[]           = "density";
char txtDensityU[]          = "density u";
char txtDensityV[]          = "density v";
char txtLevel[]             = "level";
char txtSubdivision[]       = "subdivision";
char txtBlending[]          = "blending";
char txtSpecialNets[]       = "special nets";
char txtUniform[]           = "uniform";
char txtUniformU[]          = "uniform u";
char txtUniformV[]          = "uniform v";
char txtClosed[]            = "closed";
char txtClosedU[]           = "closed u";
char txtClosedV[]           = "closed v";
char txtFlip[]              = "flip";
char txtCoons[]             = "Coons";
char txtBezier[]            = "Bezier";
char txtQKnots[]            = "quad. knots";
char txtMaxIter[]           = "max iter.";
char txtStartFrom[]         = "start from";
char txtLevels[]            = "levels";
char txtBlocks[]            = "blocks";
char txtConstraints[]       = "constraints";
char txtShapeOnly[]         = "shape only";
char txtAltML[]             = "alt. multilevel";
char txtUseCoarseMesh[]     = "use coarse mesh";
char txtOptimize[]          = "optimize";
char txtPreTransf[]         = "pretransf.";
char txtPreTransformation[] = "Pre-transformation";
char txtG1[]                = "G1";
char txtG2[]                = "G2";
char txtG1quasiG2[]         = "G1 quasi G2";
char txtRender[]            = "render";
char txtCamera[]            = "camera";
char txtInterrupt[]         = "interrupt";
char txtLight[]             = "light";
char txtC[]                 = "c.";
char txtD_shape_func[]      = "d. shape func.";
char txtGaussian[]          = "Gaussian";
char txtMean[]              = "mean";
char txtIsophotes[]         = "isophotes";
char txtReflection[]        = "reflection";
char txtHighlight[]         = "highlight";
char txtSections[]          = "sections";
char txtReflectionFrame[]   = "refl. frame";
char txtHighlightFrame[]    = "highlight frame";
char txtSectionsNormal[]    = "sections normal";
char txtShadows[]           = "shadows";
char txtAntialias[]         = "antialias";
char txtShapeFunction[]     = "shape f.";
char txtLight0[]            = "light 0";
char txtLight1[]            = "light 1";
char txtLight2[]            = "light 2";
char txtLight3[]            = "light 3";
char txtAmbient[]           = "ambient";
char txtGeneral[]           = "general";
char txtBlendingG1[]        = "blending G1";
char txtBlendingG2[]        = "blending G2";
char txtSwept[]             = "swept";
char txtSpherical[]         = "spherical";
char txtLofted[]            = "lofted";
char txtBiharmonic[]        = "biharmonic";
char txtTriharmonic[]       = "triharmonic";
char txtClamped[]           = "clamped";
char txtColour[]            = "colour";
char txtSelectVertex[]      = "select vertex";
char txtSelectEdge[]        = "select edge";
char txtOnLine[]            = "on line";
char txtOnPlane[]           = "on plane";
char txtOnSphere[]          = "on sphere";
char txtOnCylinder[]        = "on cylinder";
char txtRefPoint[]          = "ref. point";
char txtRefVector[]         = "ref. vector";
char txtRadius[]            = "radius:";
char txtAngleDeg[]          = "angle (deg):";
char txtX[]                 = "x:";
char txtY[]                 = "y:";
char txtZ[]                 = "z:";
char txtProject[]           = "project";
char txtLogIt[]             = "log it";
char txtNPThreads[]         = "#pthreads";
char txtSpecials[]          = "specials";
char txtQuestionmark[]      = "?";
char txtReset[]             = "Reset";
char txtPosition[]          = "position";
char txtF[]                 = "f:";
char txtEquator[]           = "equator";
char txtMeridian[]          = "meridian";
char txtSphericalProduct[]  = "Spherical product";
char txtG1BlendingPatch[]   = "G1 blending patch";
char txtG2BlendingPatch[]   = "G2 blending patch";


char *InfoMsg[5] =
  { "pozwalaj - a simple curve and surface modeller",
    "     *** a program under development ***",
    "This program is a part of the BSTools package.",
    "(C) Copyright by Przemyslaw Kiciak, 2010, 2013.",
    NULL };


char MsgReconsider[] = "Please, reconsider it and then proceed.";

char ErrorMsgCannotOpen[]             = "Error: Cannot read this file.";
char ErrorMsgCannotSave[]             = "Error: Cannot write the file.";
char ErrorMsgIncorrectFilename[]      = "Error: Incorrect file name.";
char ErrorMsgCannotRemoveVertex[]     = "Error: Cannot remove this vertex.";
char ErrorMsgCannotContractEdge[]     = "Error: Cannot contract this edge.";
char ErrorMsgCannotGlueLoops[]        = "Error: Cannot glue boundary loops.";
char ErrorMsgCannotRemoveFacet[]      = "Error: Cannot remove this facet.";
char ErrorMsgCannotDoubleFacetEdges[] = "Error: Cannot double edges of this facet.";
char ErrorMsgMeshIntegrity[]          = "Error: Mesh integrity test failed.";
char ErrorMsgCannotClose[]            = "Error: Not enough knots to close.";
char ErrorMsgCannotRender[]           = "Error: Only 3D objects may be rendered.";
char ErrorMsgCannotLaunchAChild[]     = "Error: Cannot launch a child process.";
char ErrorMsgChildProcessBusy[]       = "Error: Child process is busy now.";
char ErrorMsgChildProcessFailure[]    = "Error: Child process failure.";
char ErrorMsgNotImplemented[]         = "Error: This is not implemented yet.";
char ErrorMsgCannotDoIt[]             = "Error: Cannot do it.";
char ErrorMsgOptimizationFailure[]    = "Error: Optimization procedure failed.";
char ErrorMsgNothingToMark[]          = "Error: Nothing to mark.";
char ErrorMsgMeshCannotBeOptimized[]  = "Error: This mesh cannot be optimized.";
char ErrorMessageCannotCangeDegree[]  = "Error: Cannot change degree.";
char ErrorMsgCannotFlip[]             = "Error: Cannot flip the parameters.";
char ErrorMsgCannotFindObject[]       = "Error: Cannot find this object.";

