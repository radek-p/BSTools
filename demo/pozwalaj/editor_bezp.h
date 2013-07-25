
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

typedef struct {
    geom_object me;
    int         degree_u, degree_v;
    double      *cpoints;
    int         savedsize;
    double      *savedcpoints;
    byte        *mkcp;
    boolean     rational;
    int         dens_u, dens_v;
    boolean     view_surf, view_cnet;
  } GO_BezierPatch;

/* display list masks */
#define BEZP_DLM_CNET  DLISTMASK_CPOINTS /*0x0001*/
#define BEZP_DLM_PATCH 0x0002

/* number of display lists */
#define BEZP_NDL       2


/* Bezier patches processing procedures */
boolean GeomObjectInitBezierPatch ( GO_BezierPatch *obj,
                                    char spdimen, boolean rational );
geom_object *GeomObjectAddBezierPatch ( const char *name,
                                        char spdimen, boolean rational );
geom_object *GeomObjectCopyBezierPatch ( GO_BezierPatch *obj );
void GeomObjectDeleteBezierPatch ( GO_BezierPatch *obj );
void GeomObjectDrawBezierPatch ( GO_BezierPatch *obj );
void GeomObjectDrawBezierPNet ( GO_BezierPatch *obj );
void GeomObjectDisplayBezierPatch ( GO_BezierPatch *obj );
boolean GeomObjectBezierPatchSetDensityU ( GO_BezierPatch *obj, int densu );
boolean GeomObjectBezierPatchSetDensityV ( GO_BezierPatch *obj, int densv );
boolean GeomObjectBezierPatchSetDegreeU ( GO_BezierPatch *obj, int degu );
boolean GeomObjectBezierPatchSetDegreeV ( GO_BezierPatch *obj, int degv );
boolean GeomObjectBezierPatchFlipUV ( GO_BezierPatch *obj );
void GeomObjectBezierPatchFindBBox ( GO_BezierPatch *obj, Box3d *box );
boolean GeomObjectBezierPatchFindCPoint ( GO_BezierPatch *obj,
                                          CameraRecd *CPos, short x, short y,
                                          int *dist );
void GeomObjectBezierPatchSetCPoint ( GO_BezierPatch *obj,
                                      CameraRecd *CPos, short x, short y );
void GeomObjectBezierPatchMarkCPoints ( GO_BezierPatch *obj,
                                        CameraRecd *CPos, Box2s *box,
                                        char mask, boolean clear );
void GeomObjectBezierPatchMarkCPoint ( GO_BezierPatch *obj,
                                      char mask, boolean clear );
boolean GeomObjectBezierPatchSaveCPoints ( GO_BezierPatch *obj );
void GeomObjectBezierPatchUndoLastTransformation ( GO_BezierPatch *obj );
void GeomObjectBezierPatchTransformCPoints ( GO_BezierPatch *obj,
                                             trans3d *tr, char mask );

boolean GeomObjectBezierPatchGetPointCoord ( GO_BezierPatch *obj, int p,
                          int *spdimen, int *cpdimen, double **pc );

boolean GeomObjectWriteBezierPatch ( GO_BezierPatch *obj );
void GeomObjectReadBezierPatch ( void *usrdata,
                                 const char *name, int degreeu, int degreev,
                                 int pitch, const point4d *cpoints, int spdimen,
                                 boolean rational, byte *_mkcp );
void GeomObjectBezierPatchOutputToRenderer3D ( GO_BezierPatch *obj );

void GeomObjectBezierPatchDisplayInfoText ( GO_BezierPatch *obj );

boolean GeomObjectBezierPatchProcessDep ( GO_BezierPatch *obj, geom_object *go );
void GeomObjectBezierPatchProcessDeletedDep ( GO_BezierPatch *obj, geom_object *go );

