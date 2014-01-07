
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

typedef struct {
    geom_object me;
    int         degree;
    double      *cpoints;
    double      *weightpoints;
    int         savedsize;
    double      *savedcpoints;
    byte        *mkcp;
    boolean     rational;
    boolean     view_curve, view_cpoly;
  } GO_BezierCurve;

/* display list masks */
#define BEZC_DLM_CPOLY  DLISTMASK_CPOINTS /*0x0001*/
#define BEZC_DLM_CURVE  0x0002

/* number of display lists */
#define BEZC_NDL        2


/* Bezier curves processing procedures */
boolean GeomObjectInitBezierCurve ( GO_BezierCurve *obj,
                                    char spdimen, boolean rational );
geom_object *GeomObjectAddBezierCurve ( const char *name,
                                        char spdimen, boolean rational );
geom_object *GeomObjectCopyBezierCurve ( GO_BezierCurve *obj );
void GeomObjectDeleteBezierCurve ( GO_BezierCurve *obj );
void GeomObjectDrawBezierCurve ( GO_BezierCurve *obj );
void GeomObjectDrawBezierCPoly ( GO_BezierCurve *obj );
void GeomObjectDisplayBezierCurve ( GO_BezierCurve *obj );
boolean GeomObjectBezierCurveSetDegree ( GO_BezierCurve *obj, int deg );
void GeomObjectBezierCurveFindBBox ( GO_BezierCurve *obj, Box3d *box );
boolean GeomObjectBezierCurveFindCPoint ( GO_BezierCurve *obj,
                                          CameraRecd *CPos, short x, short y,
                                          int *dist );
void GeomObjectBezierCurveSetCPoint ( GO_BezierCurve *obj,
                                      CameraRecd *CPos, short x, short y );
void GeomObjectBezierCurveMarkCPoints ( GO_BezierCurve *obj,
                                        CameraRecd *CPos, Box2s *box,
                                        char mask, boolean clear );
void GeomObjectBezierCurveMarkCPoint ( GO_BezierCurve *obj,
                                       char mask, boolean clear );
boolean GeomObjectBezierCurveSaveCPoints ( GO_BezierCurve *obj );
void GeomObjectBezierCurveUndoLastTransformation ( GO_BezierCurve *obj );
void GeomObjectBezierCurveTransformCPoints ( GO_BezierCurve *obj,
                                             trans3d *tr, char mask );

boolean GeomObjectBezierCurveGetPointCoord ( GO_BezierCurve *obj, int p,
                          int *spdimen, int *cpdimen, double **pc );

boolean GeomObjectWriteBezierCurve ( GO_BezierCurve *obj );
void GeomObjectReadBezierCurve ( void *usrdata,
                                 const char *name, int degree,
                                 const point4d *cpoints, int spdimen,
                                 boolean rational );

void GeomObjectBezierCurveOutputToRenderer ( GO_BezierCurve *obj );
void GeomObjectBezierCurveDisplayInfoText ( GO_BezierCurve *obj );

boolean GeomObjectBezierCurveProcessDep ( GO_BezierCurve *obj, geom_object *go );
void GeomObjectBezierCurveProcessDeletedDep ( GO_BezierCurve *obj, geom_object *go );

