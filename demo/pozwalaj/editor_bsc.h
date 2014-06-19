
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_BSC_KNOTS   1024

#define BSC_MIN_PIPE_DIAMETER  0.01
#define BSC_MAX_PIPE_DIAMETER  1.0

typedef struct {
    geom_object me;
    int         maxknots;
    int         degree, lastknot;
    double      *knots;
    double      *cpoints;
    double      *weightpoints;
    int         savedsize;
    double      *savedcpoints;
    byte        *mkcp;
    boolean     rational, closed, uniform;
    boolean     view_curve, view_cpoly, view_bpoly,
                view_curvature, view_torsion;
    double      curvature_scale, torsion_scale;
    int         graph_dens;
    double      pipe_diameter;  /* on ray traced pictures */
  } GO_BSplineCurve;

/* display list masks */
#define BSC_DLM_CPOLY     DLISTMASK_CPOINTS /*0x0001*/
#define BSC_DLM_CURVE     0x0002
#define BSC_DLM_BPOLY     0x0004
#define BSC_DLM_CURVATURE 0x0008
#define BSC_DLM_TORSION   0x0010

/* number of display lists */
#define BSC_NDL        5


/* B-spline curves processing procedures */
boolean GeomObjectInitBSplineCurve ( GO_BSplineCurve *obj,
                                     char spdimen, boolean rational );
geom_object *GeomObjectAddBSplineCurve ( const char *name,
                                         char spdimen, boolean rational );
geom_object *GeomObjectCopyBSplineCurve ( GO_BSplineCurve *obj );
void GeomObjectDeleteBSplineCurve ( GO_BSplineCurve *obj );
boolean GeomObjectBSplineCurveSetRational ( GO_BSplineCurve *obj );
boolean GeomObjectBSplineCurveSetNonRational ( GO_BSplineCurve *obj );
void GeomObjectDrawBSplineCurve ( GO_BSplineCurve *obj );
void GeomObjectDrawBSplineCPoly ( GO_BSplineCurve *obj );
void GeomObjectDrawBSplineBPoly ( GO_BSplineCurve *obj );
void GeomObjectDrawBSplineCurvature ( GO_BSplineCurve *obj );
void GeomObjectDrawBSplineTorsion ( GO_BSplineCurve *obj );
void GeomObjectDisplayBSplineCurve ( GO_BSplineCurve *obj );
boolean GeomObjectBSplineCurveSetDegree ( GO_BSplineCurve *obj, int deg );
void GeomObjectBSplineCurveFindBBox ( GO_BSplineCurve *obj, Box3d *box );
boolean GeomObjectBSplineCurveFindCPoint ( GO_BSplineCurve *obj,
                                           CameraRecd *CPos, short x, short y,
                                           int *dist );
void GeomObjectBSplineCurveSetCPoint ( GO_BSplineCurve *obj,
                                       CameraRecd *CPos, short x, short y );
void GeomObjectBSplineCurveMarkCPoints ( GO_BSplineCurve *obj,
                                         CameraRecd *CPos, Box2s *box,
                                         char mask, boolean clear );
void GeomObjectBSplineCurveMarkCPoint ( GO_BSplineCurve *obj,
                                        char mask, boolean clear );
boolean GeomObjectBSplineCurveSaveCPoints ( GO_BSplineCurve *obj );
void GeomObjectBSplineCurveUndoLastTransformation ( GO_BSplineCurve *obj );
void GeomObjectBSplineCurveTransformCPoints ( GO_BSplineCurve *obj,
                                              trans3d *tr, char mask );

boolean GeomObjectBSplineCurveGetPointCoord ( GO_BSplineCurve *obj, int p,
                          int *spdimen, int *cpdimen, double **pc );

void GeomObjectBSplineCurveMoveKnot ( GO_BSplineCurve *obj, int kni );
boolean GeomObjectBSplineCurveInsertKnot ( GO_BSplineCurve *obj,
                                           double knot );
boolean GeomObjectBSplineCurveRemoveKnot ( GO_BSplineCurve *obj,
                                           int kni );
boolean GeomObjectBSplineCurveSetUniformKnots ( GO_BSplineCurve *obj, boolean uniform );
boolean GeomObjectBSplineCurveRefine ( GO_BSplineCurve *obj );
boolean GeomObjectBSplineCurveSetClosed ( GO_BSplineCurve *obj, boolean closed );

void GeomObjectBSplineCurveSetCurvatureGraph ( GO_BSplineCurve *obj,
                                boolean view_curvature, double curvature_scale,
                                boolean view_torsion, double torsion_scale,
                                int graph_dens );

boolean GeomObjectWriteBSCAttributes ( GO_BSplineCurve *obj );
boolean GeomObjectWriteBSplineCurve ( GO_BSplineCurve *obj );
boolean GeomObjectBSCResolveDependencies ( GO_BSplineCurve *obj );
void GeomObjectReadBSplineCurve ( void *usrdata, const char *name, int ident,
                                  int degree, int lastknot,
                                  const double *knots, boolean closed,
                                  const point4d *cpoints, int spdimen,
                                  boolean rational );

void GeomObjectBSplineCurveOutputToRenderer ( GO_BSplineCurve *obj );
void GeomObjectBSplineCurveDisplayInfoText ( GO_BSplineCurve *obj );

boolean GeomObjectBSplineCurveProcessDep ( GO_BSplineCurve *obj, geom_object *go );
void GeomObjectBSplineCurveProcessDeletedDep ( GO_BSplineCurve *obj, geom_object *go );

GO_BSplineCurve *GeomObjectFindBSCurve2DByName ( char *name );

