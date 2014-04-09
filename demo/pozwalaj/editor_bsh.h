
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

typedef struct {
    geom_object me;
    int         maxcp;
    int         hole_k;
    double      *knots;
    point2d     *domain_cp;
    double      *hole_cp;
    int         savedsize;
    double      *savedcpoints;
    byte        *mkcp;
    boolean     rational;
    boolean     view_surf, view_cnet;
  } GO_BSplineHole;


/* B-spline holes processing procedures */
boolean GeomObjectInitBSplineHole ( GO_BSplineHole *obj, boolean rational );
geom_object *GeomObjectAddBSplineHole ( const char *name, boolean rational );
geom_object *GeomObjectCopyBSplineHole ( GO_BSplineHole *obj );
void GeomObjectDeleteBSplineHole ( GO_BSplineHole *obj );
void GeomObjectDrawBSplineHole ( GO_BSplineHole *obj );
void GeomObjectDrawBSplineHNet ( GO_BSplineHole *obj );
void GeomObjectDisplayBSplineHole ( GO_BSplineHole *obj );
void GeomObjectBSplineHoleFindBBox ( GO_BSplineHole *obj, Box3d *box );
boolean GeomObjectBSplineHoleFindCPoint ( GO_BSplineHole *obj,
                                          CameraRecd *CPos, short x, short y,
                                          int *dist );
void GeomObjectBSplineHoleSetCPoint ( GO_BSplineHole *obj,
                                      CameraRecd *CPos, short x, short y );
void GeomObjectBSplineHoleMarkCPoints ( GO_BSplineHole *obj,
                                        CameraRecd *CPos, Box2s *box,
                                        byte mask, boolean clear );
void GeomObjectBSplineHoleMarkCPoint ( GO_BSplineHole *obj,
                                       byte mask, boolean clear );
boolean GeomObjectBSplineHoleSaveCPoints ( GO_BSplineHole *obj );
void GeomObjectBSplineHoleUndoLastTransformation ( GO_BSplineHole *obj );
void GeomObjectBSplineHoleTransformCPoints ( GO_BSplineHole *obj,
                                             trans3d *tr, byte mask );

boolean GeomObjectBSplineHoleGetPointCoord ( GO_BSplineHole *obj, int p,  
                          int *spdimen, int *cpdimen, double **pc );

boolean GeomObjectWriteBSHAttributes ( GO_BSplineHole *obj );
boolean GeomObjectWriteBSplineHole ( GO_BSplineHole *obj );
boolean GeomObjectBSHResolveDependencies ( GO_BSplineHole *obj );
void GeomObjectReadBSplineHole ( void *usrdata, const char *name, int ident,
                                 int hole_k, const double *knots,
                                 const point2d *domain_cp, const point4d *hole_cp,
                                 int spdimen, boolean rational );

void GeomObjectBSplineHoleDisplayInfoText ( GO_BSplineHole *obj );

boolean GeomObjectBSplineHoleProcessDep ( GO_BSplineHole *obj, geom_object *go );
void GeomObjectBSplineHoleProcessDeletedDep ( GO_BSplineHole *obj, geom_object *go );

