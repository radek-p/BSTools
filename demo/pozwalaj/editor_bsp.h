
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_BSP_KNOTS    256

/* various types of B-spline tensor product patches*/
#define BSP_TYPE_GENERAL     0
#define BSP_TYPE_SWEPT       1
#define BSP_TYPE_SPHERICAL   2
#define BSP_TYPE_LOFTED      3
#define BSP_TYPE_BLENDING_G1 4
#define BSP_TYPE_BLENDING_G2 5

/* blending patch optimization parameter range */
#define BSP_BL_MIN_CPARAM  0.0001
#define BSP_BL_MAX_CPARAM  0.1
#define BSP_BL_MIN_NKN     3
#define BSP_BL_MAX_NKN    10
#define BSP_BL_MAX_ITER   50

/* display list masks */
#define BSP_DLM_CNET  DLISTMASK_CPOINTS /*0x0001*/
#define BSP_DLM_PATCH 0x0002

/* number of display lists */
#define BSP_NDL       2

typedef struct {
    geom_object me;
    short int   bsp_type;
    int         maxknots_u, maxknots_v;
    int         degree_u, lastknot_u, degree_v, lastknot_v;
    double      *knots_u, *knots_v;
    double      *cpoints;
    int         savedsize;
    double      *savedcpoints;
    byte        *mkcp;
    boolean     rational, closed_u, closed_v, uniform_u, uniform_v,
                clamped, nharmonic;
    int         dens_u, dens_v;
    boolean     view_surf, view_cnet;
        /* equator and meridian name for the spherical product */
    char        eqname[MAX_NAME_LENGTH+1],
                mername[MAX_NAME_LENGTH+1];
        /* data for optimization of blending patches */
    int         blp_range[4];
    double      *blp_amat, **blp_arow;  /* discretized biharmonic */
    int         *blp_prof;              /* or triharmonic operator matrix */
    int         blp_lknu, blp_lknv, blp_n;  /* control net dimensions */
                                        /* and number of optimized control points */
    byte        blp_deg;
    boolean     blp_closed_u;
    double      blp_C;
    int         nkn1, nkn2, maxit;      /* nonlinear optimization parameters */
  } GO_BSplinePatch;


/* B-spline patches processing procedures */
boolean GeomObjectInitBSplinePatch ( GO_BSplinePatch *obj,
                                     char spdimen, boolean rational );
geom_object *GeomObjectAddBSplinePatch ( const char *name,
                                         char spdimen, boolean rational );
geom_object *GeomObjectCopyBSplinePatch ( GO_BSplinePatch *obj );
void GeomObjectDeleteBSplinePatch ( GO_BSplinePatch *obj );
void GeomObjectAssignBSPatch ( GO_BSplinePatch *obj, int spdimen, boolean rational,
                               int degree_u, int lastknot_u, double *knots_u,
                               int degree_v, int lastknot_v, double *knots_v,
                               double *cpoints, byte *mkcp,
                               boolean closed_u, boolean closed_v );
boolean GeomObjectBSplinePatchSetRational ( GO_BSplinePatch *obj, boolean rational );
void GeomObjectDrawBSplinePatch ( GO_BSplinePatch *obj );
void GeomObjectDrawBSplinePNet ( GO_BSplinePatch *obj );
void GeomObjectDisplayBSplinePatch ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchSetDensityU ( GO_BSplinePatch *obj, int densu );
boolean GeomObjectBSplinePatchSetDensityV ( GO_BSplinePatch *obj, int densv );
boolean GeomObjectBSplinePatchSetDegreeU ( GO_BSplinePatch *obj, int degu );
boolean GeomObjectBSplinePatchSetDegreeV ( GO_BSplinePatch *obj, int degv );
void GeomObjectBSplinePatchFindBBox ( GO_BSplinePatch *obj, Box3d *box );
boolean GeomObjectBSplinePatchFindCPoint ( GO_BSplinePatch *obj,
                                           CameraRecd *CPos, short x, short y,
                                           int *dist );
void GeomObjectBSplinePatchSetCPoint ( GO_BSplinePatch *obj,
                                       CameraRecd *CPos, short x, short y );
void GeomObjectBSplinePatchMarkCPoints ( GO_BSplinePatch *obj,
                                         CameraRecd *CPos, Box2s *box,
                                         byte mask, boolean clear );
void GeomObjectBSplinePatchMarkCPoint ( GO_BSplinePatch *obj,
                                        byte mask, boolean clear );
boolean GeomObjectBSplinePatchSaveCPoints ( GO_BSplinePatch *obj );
void GeomObjectBSplinePatchUndoLastTransformation ( GO_BSplinePatch *obj );
void GeomObjectBSplinePatchTransformCPoints ( GO_BSplinePatch *obj,
                                              trans3d *tr, byte mask );

boolean GeomObjectBSplinePatchGetPointCoord ( GO_BSplinePatch *obj, int p,
                          int *spdimen, int *cpdimen, double **pc );

void GeomObjectBSplinePatchMoveKnotU ( GO_BSplinePatch *obj, int kni );
boolean GeomObjectBSplinePatchInsertKnotU ( GO_BSplinePatch *obj,
                                            double knot );
boolean GeomObjectBSplinePatchRemoveKnotU ( GO_BSplinePatch *obj,
                                           int kni );
void GeomObjectBSplinePatchMoveKnotV ( GO_BSplinePatch *obj, int kni );
boolean GeomObjectBSplinePatchInsertKnotV ( GO_BSplinePatch *obj,
                                            double knot );
boolean GeomObjectBSplinePatchRemoveKnotV ( GO_BSplinePatch *obj,
                                           int kni );
boolean GeomObjectBSplinePatchSetUniformKnotsU ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchSetUniformKnotsV ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchSetClosedU ( GO_BSplinePatch *obj, boolean closed );
boolean GeomObjectBSplinePatchSetClosedV ( GO_BSplinePatch *obj, boolean closed );
boolean GeomObjectBSplinePatchFlipUV ( GO_BSplinePatch *obj );
boolean GeomObjectWriteBSplinePatch ( GO_BSplinePatch *obj );
void GeomObjectReadBSplinePatch ( void *usrdata,
                 const char *name,
                 int degreeu, int lastknotu, const double *knotsu,
                 int degreev, int lastknotv, const double *knotsv,
                 boolean closed_u, boolean closed_v,
                 int pitch, const point4d *cpoints,
                 int spdimen, boolean rational, byte *_mkcp );
void GeomObjectBSplinePatchOutputToRenderer3D ( GO_BSplinePatch *obj );

void GeomObjectBSplinePatchAdjustGeneral ( GO_BSplinePatch *obj );

boolean GeomObjectBSplinePatchGenSphericalProduct ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchAdjustSProduct ( GO_BSplinePatch *obj );

void GeomObjectBSplinePatchDestroyBlMat ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchIsClamped ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchInitBlG1 ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchAdjustBlG1 ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchInitBlG2 ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchAdjustBlG2 ( GO_BSplinePatch *obj );
void GeomObjectBSplinePatchSetEntireBlendingRange ( GO_BSplinePatch *obj );
void GeomObjectBSplinePatchFreeToClamped ( GO_BSplinePatch *obj );
void GeomObjectBSplinePatchClampedToFree ( GO_BSplinePatch *obj );
void GeomObjectBSPatchMarkCPGeneral ( GO_BSplinePatch *obj );
void GeomObjectBSplinePatchMarkBLRange ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchRefine ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchSetupNHBlMat ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchUpdNHarmonic ( GO_BSplinePatch *obj );
boolean GeomObjectBSplinePatchSetNHarmonic ( GO_BSplinePatch *obj,
                                             boolean nharmonic );

void GeomObjectBSplinePatchDisplayInfoText ( GO_BSplinePatch *obj );

boolean GeomObjectBSplinePatchProcessDep ( GO_BSplinePatch *obj, geom_object *go );
void GeomObjectBSplinePatchProcessDeletedDep ( GO_BSplinePatch *obj, geom_object *go );

