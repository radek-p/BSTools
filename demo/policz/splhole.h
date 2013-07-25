
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#define MAX_CPOINTS (1+12*GH_MAX_K)

#define SWIN_SELECT_VIEW          0
#define SWIN_EDITING_SURFACE      1
#define SWIN_EDITING_CONSTRAINTS  2
#define SWIN_EDITING_LIGHT        3

#define Q2COEFF_MIN   1.0
#define Q2COEFF_MAX 100.0

#define MAX_CONSTRAINTS     257 /* see also the file "edpolep.h" */
#define MAX_CONSTR_FRAME_CP  20 /* 7+1+NK_MAX*M1_MAX for G^2 */

#define MAX_SAVED_CPOINTS max ( MAX_CPOINTS, MAX_CONSTRAINTS )

typedef struct {
          int      order;                    /* 1 -- G1, 2 -- G2 */
          boolean  coons, bezier, spline;    /* basis type */
          boolean  restricted;
          boolean  lin, quasiG2;             /* functional to minimise */
          boolean  altcentre;                /* central point construction */
          boolean  gausslegendre;            /* quadrature type */
          int      nk, m1, m2;               /* spline basis parameters */
          double   slp, q2coeff;             /* coefficient for the */
                                             /* G1 quasi G2 construction */
          char     constr_type;              /* 0 - constraints off */
                                             /* 1 - first type */
                                             /* 2 - second type */
          char     nconstrsw;                /* number of constraint switches */
          boolean  constrsw[MAX_CONSTRAINTS];
          char     ncp;                      /* number of control points per curve */
          unsigned char mccp[MAX_CONSTRAINTS];
          point3d  constrcp[GH_MAX_K*MAX_CONSTR_FRAME_CP];
          vector3d constrnv;                 /* "normal" vector for the second */
                                             /* order constraints */
          short    constr_no;                /* number of constraints set */
          char     constr_spec[MAX_CONSTRAINTS][2];
          double   *bcmat;
          boolean  spl_basis_valid;
          boolean  constr_matrix_valid;
          boolean  fast;                     /* true if the construction is fast enough */
          boolean  surf_valid;               /* true if the surface is ready to display */
        } GHoptions;

typedef void (*OutputBezPatch3) ( int n, int m, const point3d *cp,
                                  void *usrptr );
typedef void (*OutputBSPatch3) ( int n, int lknu, const double *knu,
                                 int m, int lknv, const double *knv,
                                 const point3d *cp, void *usrptr );

extern xge_3Dwind swind;
extern xge_2Dwind domwind;
extern ghKnotWind knwind;
extern boolean swind_picture;  /* if nonzero, a ray-traced picture is on screen */
extern char swind_ed_switch;

extern GHoptions  options1, options2;
extern char    constr_surfno;

extern int     hole_k;
extern int     nctrlp;
extern double  knots[GH_MAX_K*11];
extern point2d domain_cp[1+12*GH_MAX_K];
extern point3d hole_cp[1+12*GH_MAX_K];
extern unsigned char mkdcp[MAX_CPOINTS];
extern unsigned char mkhcp[MAX_CPOINTS];
extern point3d saved_cp[MAX_SAVED_CPOINTS];

extern GHoleDomaind *domain1, *domain2;
extern double *acoeff1, *acoeff2;

extern point2d rpoints[4][MAX_CPOINTS];
extern point2d rdpoints[MAX_CPOINTS];  

/* final patches - Bezier or B-spline */
extern int final_np1, final_deg1, final_lkn1;
extern point3d *final_cp1;
extern double  *final_knots1;
extern int final_np2, final_deg2, final_lkn2;
extern point3d *final_cp2;
extern double  *final_knots2;

/* domain patches - always Bezier */
extern int domain_np1, domain_deg1;
extern point2d *domain_bcp1;
extern int domain_np2, domain_deg2;
extern point2d *domain_bcp2;

/* parameters for the built-in constructions of the domain and surface */
extern vector2d domcvect[GH_MAX_K];
extern double   domcparam[];
extern vector3d surfcvect[GH_MAX_K];
extern double   surfcparam[];

boolean Option1SetOrder ( int k );
boolean Option1SetNK ( int k );
boolean Option1SetM1 ( int k );
boolean Option1SetM2 ( int k );
boolean Option2SetOrder ( int k );
boolean Option2SetNK ( int k );
boolean Option2SetM1 ( int k );
boolean Option2SetM2 ( int k );

int MyG1OptionProc1 ( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata );
int MyG2OptionProc1 ( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata );
int MyG1OptionProc2 ( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata );
int MyG2OptionProc2 ( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata );

void OutputDomainPatch1 ( int n, int m, const point2d *cp );
void OutputDomainPatch2 ( int n, int m, const point2d *cp );
boolean UpdateDomain1 ( void );
boolean UpdateDomain2 ( void );
boolean UpdateDomains ( void );
boolean InitGHObject ( int k );

void ProjectSurfaceNet ( void );
void ProjectDomainNet ( void );
void FindBoundingBox ( Box3d *bbox );
void FindDomainBoundingBox ( Box2d *bbox );

void MyDrawLined ( point2d *p0, point2d *p1 );
void MyMarkPointd ( point2d *p );
void MyMarkPoint2s ( point2s *p );
void MyWriteNumber ( short x, short y, int n );

void DrawGHControlNet ( int id );
boolean GetSurfBezPatch ( int i, int j, point3d *bcp );
void DrawGHSurfPatch ( int n, int m, const point3d *cp );
void DrawGHSurfPatches ( int id );
void DrawGHSurfFillingPatches ( int id, int sn );
void DrawGHSurfNumbers ( int id );

void DrawGHDomainControlNet ( void );
boolean GetDomBezPatch ( int i, int j, point2d *bcp );
void DrawDomSurrndPatch ( int n, int m, const point2d *cp );
void DrawDomainSurrPatches ( void );
void DrawDomainPatches ( int dn );
void DrawDomainNumbers ( void );

boolean FindNearestPoint ( int npts, point2d *pts,
                           short x, short y, short mindist,
                           int *current_point );

boolean FindNearestSurfCPoint ( int id, short x, short y, short mindist );
void SetSurfCPoint ( int id, short x, short y );
void SelectCPoints ( Box2s *sel_rect, int npts, point2d *rpts,
                     boolean *mkpts, boolean value );
void SelectSurfCPoints ( int id );
void UnselectSurfCPoints ( int id );
void SaveSurfCPoints ( void );
void TransformSurfCPoints ( void );

boolean FindNearestDomainCPoint ( short x, short y, short mindist );
void SetDomainCPoint ( short x, short y );
void SelectDomainCPoints ( void );
void UnselectDomainCPoints ( void );
void SaveDomainCPoints ( void );
void TransformDomainCPoints ( void );

void InvalFinalSurface1 ( void );
void InvalFinalSurface2 ( void );
void InvalFinalSurfaces ( void );
void OutputFinalBezPatch1 ( int n, int m, const double *cp, void *usrptr );
void OutputFinalBezPatch2 ( int n, int m, const double *cp, void *usrptr );
void OutputFinalBSPatch1 ( int n, int lknu, const double *knu,
                           int m, int lknv, const double *knv,
                           const double *cp, void *usrptr );
void OutputFinalBSPatch2 ( int n, int lknu, const double *knu,
                           int m, int lknv, const double *knv,
                           const double *cp, void *usrptr );
boolean UpdateFinalSurfaceNC ( int surfno, GHoleDomaind *domain, GHoptions *options,
                  double **acoeff,
                  void (*OutputBezPatch) ( int n, int m, const double *cp,
                                           void *usrptr ),
                  void (*OutputBSPatch) ( int n, int lknu, const double *knu,
                                          int m, int lknv, const double *knv,
                                          const double *cp, void *usrptr ),
                  int *final_np );
boolean UpdateFinalSurfaceC1 ( int surfno, GHoleDomaind *domain, GHoptions *options,
                  double **acoeff,
                  void (*OutputBezPatch) ( int n, int m, const double *cp,
                                           void *usrptr ),
                  void (*OutputBSPatch) ( int n, int lknu, const double *knu,
                                          int m, int lknv, const double *knv,
                                          const double *cp, void *usrptr ),
                  int *final_np );
boolean UpdateFinalSurfaceC2 ( int surfno, GHoleDomaind *domain, GHoptions *options,
                  double **acoeff,
                  void (*OutputBezPatch) ( int n, int m, const double *cp,
                                           void *usrptr ),
                  void (*OutputBSPatch) ( int n, int lknu, const double *knu,
                                          int m, int lknv, const double *knv,
                                          const double *cp, void *usrptr ),
                  int *final_np );
boolean UpdateFinalSurface1 ( void );
boolean UpdateFinalSurface2 ( void );
boolean SurfFast ( GHoptions *opt );

void SendPatchesToRenderer ( void );

boolean FindNearestLightPoint ( int id, short x, short y, short mindist );
void SetLightPoint ( int id, short x, short y );
void DrawLightPoints ( int id );

void EditSetSelectView ( void );
void EditSetSurface ( void );
void EditSetConstraints ( void );
void EditSetLight ( void );

boolean FindNearestSWinPoint ( int id, short x, short y, short mindist );
void SetSWinPoint ( int id, short x, short y );
void SelectSWinPoints ( int id );
void UnselectSWinPoints ( int id );
void SaveSWinPoints ( void );
void TransformSWinPoints ( void );

void InitConstraintFrame ( char surfno );
boolean GetCurrentConstraintFrame ( char surfno );
void DrawConstraintFrame ( int id, char surfno );
boolean FindNearestConstraintCPoint ( int id, short x, short y, short mindist );
void SetConstraintCPoint ( int id, short x, short y );
void SelectConstraintCPoints ( int id );
void UnselectConstraintCPoints ( int id );
void SaveConstraintCPoints ( void );
void TransformConstraintCPoints ( void );

void CountTheConstraints ( void );
boolean ComputeConstraintMatrices ( int surfno );
boolean ComputeConstraintRightSide ( int surfno, vector3d *rhs );
boolean ComputeConstraintMatricesAlt ( int surfno );
boolean ComputeConstraintRightSideAlt ( int surfno, double *altrhs );

