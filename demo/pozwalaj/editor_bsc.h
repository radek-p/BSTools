
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_BSC_KNOTS   1024

#define BSC_MIN_PIPE_DIAMETER  0.01
#define BSC_MAX_PIPE_DIAMETER  1.0

#define BSC_MC_MIN_EXP    3.2
#define BSC_MC_MAX_EXP   32.0
#define BSC_MC_MIN_PPAR   1.0e-2
#define BSC_MC_MAX_PPAR   1.0e+12

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
    boolean     rational, closed, uniform, mengerc;
    boolean     view_curve, view_cpoly, view_bpoly,
                view_curvature, view_torsion;
    double      curvature_scale, torsion_scale;
    int         graph_dens;
    double      pipe_diameter;  /* on ray traced pictures */
        /* Menger curvature related stuff */
    double      mc_exponent, mc_pparam[MENGERC_NPPARAM];
    int         mc_nqkn;
    int         mc_ppopt;
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
void GeomObjectBSplineCurveInitMC ( GO_BSplineCurve *obj );
void GeomObjectBSplineCurveCopyMC ( GO_BSplineCurve *obj, GO_BSplineCurve *copy );
boolean GeomObjectInitBSplineCurve ( GO_BSplineCurve *obj,
                                     char spdimen, boolean rational );
geom_object *GeomObjectAddBSplineCurve ( const char *name,
                                         char spdimen, boolean rational );
geom_object *GeomObjectCopyBSplineCurve ( GO_BSplineCurve *obj );
void GeomObjectDeleteBSplineCurve ( GO_BSplineCurve *obj );
void GeomObjectAssignBSCurve ( GO_BSplineCurve *obj, int spdimen, boolean rational,
                               int degree, int lastknot, double *knots,
                               double *cpoints, byte *mkcp, boolean closed );
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
                                         char mask, int action );
void GeomObjectBSplineCurveMarkCPoint ( GO_BSplineCurve *obj,
                                        char mask, int action );
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
boolean GeomObjectBSplineCurveSetMengerc ( GO_BSplineCurve *obj, boolean mengerc );

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

