
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

/* some arbitrary limitations */
#define MAX_NAME_LENGTH   64

  /* maximal degree of polynomial curves and patches */
#define MAX_DEGREE        10
  /* maximal density of nets of constant parameter curves for */
  /* wireframe polynomial patches */
#define MAX_PNET_DENSITY  16
  /* maximal subdivision level for images of subdivision surfaces */
#define MAX_SUBDIV_LEVEL   5
  /* maximal number of dependencies */
#define MAX_DEPNUM       256

/* object type identifiers */
#define GO_BEZIER_CURVE    1
#define GO_BSPLINE_CURVE   2
#define GO_BEZIER_PATCH    3
#define GO_BSPLINE_PATCH   4
#define GO_BSPLINE_MESH    5
#define GO_BSPLINE_HOLE    6

/* tolerance for picking the points */
#define MAXPIXDIST        10

/* bit masks for control points of the objects */
#define MASK_CP_MOVEABLE  0x01
#define MASK_CP_BOUNDARY  0x02
#define MASK_CP_SPECIAL   0x04  /* same as in pozwalaj_proc.h */
#define MASK_CP_MARKED_0  0x08
#define MASK_CP_MARKED_1  0x10
#define MASK_CP_MARKED_2  0x20
#define MASK_CP_MARKED_3  0x40
#define MASK_CP_MARKED_4  0x80
#define MASK_MARKED (MASK_CP_MARKED_0 | MASK_CP_MARKED_1 | \
                     MASK_CP_MARKED_2 | MASK_CP_MARKED_3 | MASK_CP_MARKED_4)

#define MARK_SELECT    0
#define MARK_UNSELECT  1
#define MARK_TGSELECT  2

/* range of coefficients for G1 quasi G2 polygonal hole filling */
#define Q2COEFF_MIN       1.0
#define Q2COEFF_MAX     100.0

/* masks for display lists - common to all types of objects */
#define DLISTMASK_CPOINTS 0x0001

/* options for writing data in files */
#define GO_WRITE_CURRENT 0
#define GO_WRITE_ACTIVE  1
#define GO_WRITE_ALL     2

/* common object attributes */
typedef struct geom_object {
    struct  geom_object *prev, *next;
    int     ident;
    int     maxdn;
    struct  geom_object **dependencies;
    struct  geom_object *dep_list;
    char    obj_type;
    char    spdimen, cpdimen;
    boolean active;
    short   dlistmask;
    GLint   displaylist;
    double  colour[3];
    char    name[MAX_NAME_LENGTH+1];
    boolean bound_with_a_child;
    boolean display_pretrans;
    trans3d pretrans;
    char    tag;
        /* dependencies read in from a file, before resolving */
    int     filedepname;
    int     filedepnum;
    int     *filedepid;
  } geom_object;

/* auxiliary data structure for reading optional object */
/* attributes from data files; */
typedef struct {
    geom_object    *go_being_read;
    int            obj_type;
    int            nmk;
    byte           *mk;
    int            nhemk;
    byte           *hemk;
    double         colour[3];
        /* dependencies read in from a file */
    int            filedepname;
    int            filedepnum;
    int            *filedepid;
        /* elements of a trimmed domain */
    int            ntrd, trdl;
    mbs_polycurved *trdelem;
    boolean        trd_error;
  } rw_object_attributes;

extern geom_object *first_go, *last_go, *current_go, *currentp_go;
extern int         current_point_ind;
extern double      *current_point;
extern char        current_point_dim;
extern boolean     current_point_rational;
extern byte        marking_mask;


void SetStatusText ( char *text, boolean onscreen );

/* object list processing procedures */
void GeomObjectInitList ( void );
void GeomObjectPurgeList ( void );
void GeomObjectDeleteCurrent ( void );
boolean GeomObjectCopyCurrent ( void );
int GeomObjectNumber ( void );
int GeomObjectCurrentNumber ( void );
boolean GeomObjectSelect ( int objno );
boolean GeomObjectIsInTheList ( geom_object *go );
void GeomObjectDisplayInfoText ( geom_object *obj );
geom_object *GeomObjectFindByName ( geom_object *obj, const char *name );
boolean GeomObjectExchangeInTheList ( boolean up );

void GeomObjectSetupIniPoints ( char spdimen, boolean rational,
                                char *cpdimen,
                                int np, const double *pt4, double *pt );

/* drawing procedures */
void DrawAVertex ( int cpdimen, int spdimen, const double *cp );
void DrawAPolyline ( int cpdimen, int spdimen, int np, const double *cp );
void DrawLineSegments ( int cpdimen, int spdimen, int np,
                        const double *cp, const double *wp );
void DrawARectNet ( int cpdimen, int spdimen, int ncols, int nrows,
                    int pitch, const double *cp );
void DrawCPoints ( int cpdimen, int spdimen, int np, const double *cp,
                   const byte *mkcp, byte mask1, byte mask2 );
void DrawPoints ( int cpdimen, int spdimen, int np, const double *cp );
void DrawBezierCurve ( int cpdimen, int spdimen, int deg, double *cp, int dd );
void DrawBSplineCurve ( int cpdimen, int spdimen, int deg, int lkn, double *knots,
                        double *cp, int dd );
void DrawBezierPatchWF ( int cpdimen, int spdimen, int degu, int degv,
                         int pitch, const double *cpoints,
                         int densu, int densv, int denscu, int denscv,
                         boolean firstu, boolean lastu,
                         boolean firstv, boolean lastv );
void DrawBSplinePatchWF ( int cpdimen, int spdimen,
                          int degu, int lknu, const double *knu,
                          int degv, int lknv, const double *knv,
                          int pitch, const double *cpoints,
                          int densu, int densv, boolean firstu, boolean lastu,
                          boolean firstv, boolean lastv );
void GeomObjectDrawMeshEdges ( int cpdimen, int spdimen,
                               int nv, BSMvertex *mv, double *mvc, int *mvhei,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               byte edgemask );
void GeomObjectDrawMeshFacet ( int cpdimen, int spdimen,
                               int nv, BSMvertex *mv, double *mvc, int *mvhei,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               int facetnum );
void GeomObjectDrawMeshHalfedge ( int cpdimen, int spdimen,
                                  int nv, BSMvertex *mv, double *mvc, int *mvhei,
                                  int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  int henum );
void GeomObjectDrawMeshHalfedges ( int cpdimen, int spdimen,
                                   int nv, BSMvertex *mv, double *mvc, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   byte edgemask,
                                   byte *hemark, byte hemask, boolean inv,
                                   boolean half0, boolean half1 );
void GeomObjectDrawMeshVertex ( int cpdimen, int spdimen,
                                int nv, BSMvertex *mv, double *mvc, int *mvhei,
                                int nhe, BSMhalfedge *mhe,
                                int nfac, BSMfacet *mfac, int *mfhei,
                                int vertnum );

/* auxiliary editing procedures */
void GeomObjectSetCurrentPointPtr ( double *currpt, char dim, boolean rational );
void GeomObjectSetupWeightPoints ( char cpdimen, int np, double *cpoints,
                                   double *weightpoints );
void GeomObjectUnmarkPoints ( int np, byte *mkcp );
void GeomObjectFindBBox ( char cpdimen, boolean rational, int np, double *cp,
                          Box3d *box );
boolean GeomObjectFindNearestPoint ( char cpdimen, char spdimen,
                                     int np, double *cp, byte *mkcp,
                                     byte mask, CameraRecd *CPos,
                                     short x, short y, int *dist );
void GeomObjectSetPoint ( CameraRecd *CPos, short x, short y,
                          char cpdimen, char spdimen, double *pt );
void GeomObjectSetWeightPoint ( CameraRecd *CPos, short x, short y,
                                char cpdimen, double *pt, double *wpt );
boolean GeomObjectPointInBox ( char cpdimen, char spdimen, double *pt,
                               CameraRecd *CPos, Box2s *box );
void GeomObjectMarkPoints ( char cpdimen, char spdimen,
                            int np, byte *mkcp, double *cp,
                            CameraRecd *CPos, Box2s *box,
                            byte mask, int action );
void GeomObjectMarkPoint ( int np, byte *mkcp, byte mask, int action );
void GeomObjectTransformPoint ( char cpdimen, char spdimen,
                                double *scp, double *cp, trans3d *tr );
void GeomObjectTransformPoints ( char cpdimen, char spdimen,
                                 int np, byte *mkcp, byte mask,
                                 double *savedp, double *cp, trans3d *tr );
boolean GeomObjectSelectPoint ( char spdimen, CameraRecd *CPos, short x, short y );
boolean GeomObjectSelectEdge ( char spdimen, CameraRecd *CPos, short x, short y );
void GeomObjectUnselectPoint ( char spdimen );
void GeomObjectUnselectEdge ( char spdimen );

/* a dispatcher */
void GeomObjectDisplayActive ( char spdimen );
boolean GeomObjectFindBoundingBox ( char spdimen, Box3d *box );
boolean GeomObjectFindObjCPoint ( char spdimen, CameraRecd *CPos, short x, short y,
                                  geom_object *obj );
boolean GeomObjectFindCPoint ( char spdimen, CameraRecd *CPos, short x, short y );
void GeomObjectSetCPoint ( CameraRecd *CPos, short x, short y );
void GeomObjectMarkCPoints ( char spdimen, CameraRecd *CPos, Box2s *box,
                             boolean clear );
void GeomObjectMarkHalfedges ( char spdimen, CameraRecd *CPos, Box2s *box,
                               int action );
void GeomObjectSetMarkingMask ( boolean bits[4] );
boolean GeomObjectSaveCPoints ( char spdimen );
void GeomObjectUndoLastTransformation ( char spdimen );
void GeomObjectTransformCPoints2D ( trans2d *tr, byte mask );
void GeomObjectTransformCPoints3D ( trans3d *tr, byte mask );

boolean GeomObjectGetPointCoord ( geom_object *go, int p,
                                  int *spdimen, int *cpdimen, double **pc );

void GeomObjectSavePretransformation ( geom_object *go );
void GeomObjectCompPretransformation ( geom_object *go, trans3d *tr );
void GeomObjectSetPretransformation ( geom_object *go, trans3d *tr );

/* special transformations */
void GeomObjectSpecial3DTranslate ( vector3d *v );
void GeomObjectSpecial3DScale ( point3d *p, vector3d *s );
void GeomObjectSpecial3DRotate ( point3d *p, vector3d *v, double phi );
void GeomObjectProject3DPointsOnALine ( point3d *p, vector3d *v );
void GeomObjectProject3DPointsOnAPlane ( point3d *p, vector3d *v, char dir );
void GeomObjectProject3DPointsOnASphere ( point3d *p, double r, char dir );
void GeomObjectProject3DPointsOnACylinder ( point3d *p, vector3d *v, double r,
                                            char dir );

/* output to the renderer */
void GeomObjectOutputToRenderer3D ( boolean all );

/* disk I/O */
void GeomObjectBeginReading ( void *usrdata, int obj_type );
void GeomObjectEndReading ( void *usrdata, int obj_type, boolean success );
void GeomObjectReadColour ( void *usrdata, point3d *colour );
void GeomObjectReadCPMK ( void *usrdata, int ncp, unsigned int *mk );
void GeomObjectReadHEMK ( void *usrdata, int nhe, unsigned int *mk );
void GeomObjectReadDependency ( void *usrdata,
                                int depname, int ndep, int *dep );
boolean GeomObjectResolveDependencies ( void );
void GeomObjectReadTrimmedDomain ( void *usrdata, mbs_polycurved *elem );
boolean GeomObjectReadFile ( char *filename, bsf_Camera_fptr CameraReader );
boolean GeomObjectWriteAttributes ( void *usrdata );
boolean GeomObjectWriteObj ( geom_object *go );
boolean GeomObjectMakeIdentifiers ( char *filename, boolean append );
boolean GeomObjectWriteFile ( char *filename, char whattowrite,
                              boolean (*writeotherdata)( void *usrdata ),
                              void *usrdata,
                              boolean append );

/* debugging etc. */
boolean GeomObjectOpenDumpFile ( void );
void GeomObjectCloseDumpFile ( void );
void GeomObjectDumpBezierPatch ( int spdim, int cpdim, boolean rational,
                                 int udeg, int vdeg, const double *cp );

/* dependencies */
boolean GeomObjectSortDependencies ( void );
boolean GeomObjectAddDependency ( geom_object *go, int maxdn,
                                  int dn, geom_object *dgo );
void GeomObjectRemoveDependency ( geom_object *go, int dn );
void GeomObjectDeleteDependencies ( geom_object *go );
boolean GeomObjectProcessDependencies ( geom_object *go );
void GeomObjectProcessDeletedDep ( geom_object *dgo );

