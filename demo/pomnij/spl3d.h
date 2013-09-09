
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_DEGREE    9
#define MAX_KNOTS   128

#define MAX_BLENDING_CONSTRAINTS  20

#define SCRATCH_MEM_SIZE 268435456 /* 256MB */

#define KNOT_MARGIN 50

#define MAX_CPOINTS ((MAX_DEGREE+1)*(MAX_KNOTS-1)*(MAX_DEGREE+1)*(MAX_KNOTS-1))

#define KNOT_EPS 1.0e-4 /* the smallest distance between neighbouring knots */

#define WIN1_GENERAL   0
#define WIN1_SPHERICAL 1
#define WIN1_SWEPT     2
#define WIN1_BLENDING  3

#define NLBLENDING_CMIN  0.0001
#define NLBLENDING_CMAX  0.1

extern int win1_contents;

/* general B-spline patch representation */
extern xge_3Dwind     swind;
extern xge_T2KnotWind kwind;

/* the degrees and numbers of knots are stored in the */
/* widget data structures above */

extern int     degree_u, lastknot_u, degree_v, lastknot_v;
extern double  knots_u[], knots_v[];
extern point4d cpoints[];
extern point4d savedcpoints[MAX_CPOINTS];
extern point2d rpoints[4][MAX_CPOINTS];
extern boolean clpoints[4][MAX_CPOINTS];
extern byte    mkpoints[MAX_CPOINTS];

/* curves for spherical product representations */
extern xge_2Dwind     eq_cwind, mer_cwind;
extern xge_KnotWind   eq_ckwind, mer_ckwind;
extern xge_int_widget eqdeg;

extern boolean equator, meridian;
extern boolean bind_spr;
extern boolean eqmer_nurbs;
extern boolean eqmer_closed;
extern double  arc_angle;

extern double  meridian_knots[];
extern point3d meridian_cpoints[];
extern point2d meridian_rpoints[];
extern byte    meridian_mkpoints[];
extern boolean meridian_nurbs;

extern double  equator_knots[];
extern point3d equator_cpoints[];
extern point2d equator_rpoints[];
extern byte    equator_mkpoints[];
extern boolean equator_nurbs;

extern int     neqmerpoints;
extern point3d *eqmer_cpoints;
extern point2d *eqmer_rpoints;
extern byte    *eqmer_mkpoints;

extern boolean display_eqmer_curve, display_eqmer_control_polygon,
               display_eqmer_Bezier_polygons, display_eqmer_ticks;

/* the variables below are used in the construction */
/* of triharmonic blending surfaces */
extern boolean sw_bind_blending;
extern boolean blending_mat_valid;
extern int     blending_n, blending_lknu, blending_lknv;
extern double  *blending_Amat;
extern double  **blending_Arow;
extern int     *blending_prof;

extern boolean sw_blending_g1, sw_blending_g2;
extern boolean sw_clamped_blending;
extern boolean sw_triharmonic_blending, sw_nonlin_blending;
extern double  blending_factor;
extern int     blending_lmt_iter;
extern boolean sw_blending_constraints, sw_blending_opt_entire;
extern int     blending_opt_part[4];
extern trans3d blending_opt_transform;
extern boolean display_trans_net;

extern boolean sw_blending_opt_dump;

/* blending constraints data */
extern int max_blending_constraints, n_blending_constraints;
extern double  blending_constr_knots[MAX_BLENDING_CONSTRAINTS+1];
extern boolean blending_constr_poly_valid[MAX_BLENDING_CONSTRAINTS+1];
extern point4d blending_constr_cp[MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];
extern point2d blending_constr_rp[4][MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];
extern boolean clblending_constr[4][MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];
extern byte    mkblending_cp[MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];
extern point4d savedbl_constr_cp[MAX_BLENDING_CONSTRAINTS*(MAX_KNOTS)];

extern boolean display_surface, display_control_net, display_Bezier_nets;
extern boolean display_domain_net;
extern boolean display_constr_poly;

extern int display_bez_dens_u, display_bez_dens_v;


/* ////////////////////////////////////////////////////////////////////////// */
void ClearPointMarking ( int npoints, boolean *mkpoints );
boolean FindNearestPoint ( int npoints, const point2d *rpoints,
                           short x, short y, short *dist, int *k );
void MarkPoints ( int npoints, const point2d *rpoints, boolean *mkpoints,
                  const Box2s *sbox, boolean mk );

/* ////////////////////////////////////////////////////////////////////////// */
boolean InsertKnotU ( double newknot );
boolean InsertKnotV ( double newknot );
boolean RemoveKnotU ( int knotnum );
boolean RemoveKnotV ( int knotnum );

boolean SetClosedUSurface ( void );
boolean SetClosedVSurface ( void );

void SetEquidistantU ( void );
void SetEquidistantV ( void );

boolean FindNearestCPoint ( int id, short x, short y );
void SetCPoint ( int id, short x, short y );
void MkPoint ( int i, int j, boolean mk );
void MarkBlConstrPoint ( int i, boolean mark );
void SelectPoints ( int id, const Box2s *sbox, boolean mk );

boolean RzutujPunkt3 ( int id, point3d *p, point2d *q );
boolean RzutujPunkt ( int id, point4d *p, point2d *q );

void ProjectSurface ( int id );
void ResetObject ( void );
void ResizeObject ( void );
void FindBoundingBox ( Box3d *box );
void FlipPatch ( void );

boolean DegreeElevationU ( void );
boolean DegreeElevationV ( void );
boolean DegreeReductionU ( void );
boolean DegreeReductionV ( void );

void DisplayLine3 ( int id, point3d *p0, point3d *p1 );
void DisplayLine4 ( int id, point4d *p0, point4d *p1 );

void DisplayBezCurve4 ( CameraRecd *CPos, int degree, point4d *cp );
void DisplayBezPatch4 ( CameraRecd *CPos,
                        int degree_u, int degree_v, point4d *cp,
                        boolean u_edge, boolean v_edge );
void DisplaySurface ( int id );
void DisplayControlNet ( int id );
void DisplayBezierNets ( int id );
void DisplayControlPoints ( int id );
void DisplayConstraintCurves ( int id );
void DisplayPreTransControlNet ( int id );

void DisplayKnots ( boolean high );
void DisplayKnotLines ( void );
void DisplayBlendingKnots ( void );
double GrevilleAbscissa ( int degree, double *knots, int i );
void DisplayDomainNet ( void );
void DisplayDomain ( void );
void SelectDomPoints ( const Box2s *sbox, boolean mk );

void SaveControlPoints ( void );
void TransformMarkedControlPoints ( trans3d *tr );

void ErrorHandler ( int module, const char *file, int line,
                    int errcode, const char *errstr );
void DumpData ( void );

/* ////////////////////////////////////////////////////////////////////////// */
boolean SetClosedEqMer ( void );
boolean FindNearestEqMerCPoint ( int x, int y, int mindist );

void ProjectEqMerCurve ( void );

void DisplayEqMerControlPolygon ( void );
void DisplayEqMerCurve ( void );
void DisplayEqMerControlPoints ( void );
void DisplayEqMerAxes ( void );
void DisplayEqMerBezierPolygons ( void );

void InsertEqMerKnot ( void );
void RemoveEqMerKnot ( void );
boolean EqMerDegreeElevation ( void );
boolean EqMerDegreeReduction ( void );

void EqMerSetKnot ( void );
void EqMerSetCPoint ( short x, short y );

void SelectEqMerPoints ( const Box2s *sbox, boolean mk );
void ClearEqMerPointMarking ( void );
void SaveEqMerControlPoints ( void );
void TransformEqMerMarkedControlPoints ( trans2d *tr );

void ResetEquator ( void );
void ResetMeridian ( void );
void ResetEqMer ( void );
void FindEqMerRefBox ( Box2d *box );

void SelectEquator ( void );
void SelectMeridian ( void );
void SetupNURBSEqMer ( void );
void BindSphericalProduct ( void );

void EqMerSetQuarterCircleArc ( void );
void EqMerSetHalfCircleArc ( void );
void EqMerSetFullCircleArc ( void );
void EqMerSetCircleArc ( double angle );

/* ////////////////////////////////////////////////////////////////////////// */
void InitG1BlendingSurface ( void );
boolean RefineG1BlendingSurface ( void );
boolean SetupG1BlendingMatrix ( void );
boolean ConstructG1BlendingSurface ( void );
void G1FreeBoundaryToClamped ( void );
void G1ClampedBoundaryToFree ( void );
void G1CheckIfClamped ( void );

void InitG2BlendingSurface ( void );
boolean RefineG2BlendingSurface ( void );
boolean SetupG2BlendingMatrix ( void );
boolean ConstructTriharmG2BlendingSurface ( void );
boolean ConstructTriharmG2BlendingConstrSurface ( void );
boolean ConstructTriharmG2BlendingConstrClosedSurface ( void );
boolean ConstructG2BlendingSurface ( void );
void G2FreeBoundaryToClamped ( void );
void G2ClampedBoundaryToFree ( void );
void G2CheckIfClamped ( void );

void ValidateG1ConstrKnots ( void );
boolean InsertG1BlendingConstrKnot ( double newknot );
boolean RemoveG1BlendingConstrKnot ( int knotnum );
void ChangeG1BlendingConstrKnot ( int oldnum, int knotnum );
void ProjectG1BlendingConstrCP ( int id, int knotnum );
void ProjectG1BlendingConstrCPoly ( int id );
void FindG1BlendingConstrCP ( int knotnum );
void DisplayG1BlendingConstrCP ( int id );

void ValidateG2ConstrKnots ( void );
boolean InsertG2BlendingConstrKnot ( double newknot );
boolean RemoveG2BlendingConstrKnot ( int knotnum );
void ChangeG2BlendingConstrKnot ( int oldnum, int knotnum );
void ProjectG2BlendingConstrCP ( int id, int knotnum );
void ProjectG2BlendingConstrCPoly ( int id );
void FindG2BlendingConstrCP ( int knotnum );
void DisplayG2BlendingConstrCP ( int id );

void SendOptionsToChild ( void );
int PatchDataSize ( void );
void SendPatchToChild ( void );
int ConstraintDataSize ( void );
void SendConstraintsToChild ( void );
void GetOptionsFromChild ( void );
void GetPatchFromChild ( void );
void GetConstraintsFromChild ( void );

/* ////////////////////////////////////////////////////////////////////////// */
boolean PovRayExport ( void );
boolean CExport ( void );

