
#define SCRATCH_MEM_SIZE 2097152  /* 2MB */

extern CameraRecf PPos[3], DomPPos;
extern CameraRecf CPos;

extern int hole_k, hole_np, nconstr;

extern float knots[GH_MAX_K][11];
extern point2f domcp[12*GH_MAX_K+1];
extern point3f surfcp[12*GH_MAX_K+1];
extern point3f constrcp[4*GH_MAX_K+1];

extern vector2f firstpartition[GH_MAX_K];
extern float surfparams[4];
extern float domparams[3];

extern boolean swCoonsPatches, swBezierPatches,
               swConstraintsOn, swPointMode, swZeroDer, swNormalConstr;
extern boolean swConstraint[4*GH_MAX_K+1];
extern boolean swZConstraint[4*GH_MAX_K+1];

extern GHoleDomainf *domain;
extern Box3f RefBBox;

extern boolean FormMatrixValid, ExtFormMatrixValid,
               ConstrMatrixValid, ExtConstrMatrixValid,
               NConstrMatrixValid, ExtNConstrMatrixValid,
               FinalSurfValid, ConstraintsValid,
               NLFinalSurfValid;

extern char constrql1, constrql2;

extern point3f FinalCP[2*GH_MAX_K][(G1H_FINALDEG+1)*(G1H_FINALDEG+1)];       
extern int finaldegu, finaldegv;       
extern int nfinal;       


void InitSurface ( int k );
void RecreateDomain ( void );
boolean UpdateFormMatrix ( void );
boolean UpdateFinalSurf ( void );
boolean UpdateExtFormMatrix ( void );
boolean UpdateNLFinalSurf ( void );

void SetupRefBox ( float x0, float x1, float y0, float y1, float z0, float z1 );
void InitProjections ( void );
void SetupPerspProj ( ed_rect *edr, boolean reset );
void ResetCPos ( void );
void SetupParProj ( ed_rect *edr );
void SetupDomPPos ( int w, int h, int x, int y, float d );

void ProjectDomainCP ( void );
void ProjectSurfCP ( void );
void ProjectDomainCentralPoint ( void );

void MyDrawLine ( point2f *p1, point2f *p2 );
void MyMarkRect ( point2f *p );
void RedrawDomainCP ( boolean bright );
void RedrawDomainSurrPatches ( void );
void RedrawDomainCurves ( void );
void RedrawDomainAuxPatches ( void );
void RedrawDomainPatches ( void );
void RedrawDomainPartition ( void );
void RedrawFirstPartition ( boolean bright );
void RedrawDomainNumbers ( void );
void RedrawDomainCentralPoint ( void );

void RedrawHoleKnots ( int w, int h, int x, int y );

boolean FindNearestKnot ( int w, int h, int xc, int yc, 
                          int x, int y, int *ink, int *jnk );
boolean SetKnot ( int w, int h, int xc, int yc, int i, int j, int x );

int FindNearestDomCP ( int x, int y );
void SetDomCP ( int np, int x, int y );

int FindNearestFirstPartitionVector ( int x, int y );
boolean SetFirstPartitionVector ( int nfp, int x, int y );

int FindNearestCPointPart ( int x, int y );
boolean SetCPointPart ( int nfp, int x, int y );

void RedrawSurfaceCP ( int id, boolean bright );
void RedrawSurface ( int id );
void RedrawFinalPatches ( int id );
void RedrawNLFinalPatches ( int id );
void RedrawSurfNumbers ( int id );
void SetCPoint ( int id, int np, int x, int y );
void ProjectScene ( int id );
int FindNearestPoint ( int id, int x, int y );
void FlattenSurface ( void );
void FitSurfaceToDomain ( void );
void FitDomainToSurface ( void );

void UpdateDerConstraintsHandle ( int np, point3f *hp );
void GetConstraintsHandle ( void );
void ProjectConstraintsHandle ( int id );
void RedrawConstraintsHandle ( int id, boolean bright );
int GetConstraintsNumber ( void );
boolean SetupConstraintsMatrix ( void );
boolean SetupExtConstraintsMatrix ( void );
void SetupConstraintsRHS ( point3f *rhs );
boolean SetupNConstraintsMatrix ( void );
boolean SetupExtNConstraintsMatrix ( void );
void SetupNConstraintsRHS ( float *rhs );
void SetDefaultConstraints ( void );
void TurnConstraints ( void );
void SwitchAConstraint ( int cno );

void WriteInfo ( void );
