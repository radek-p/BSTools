
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_DEGREE   9
#define MAX_KNOTS  100

#define MIN_KNOT_SCALE  0.01
#define MAX_KNOT_SCALE 100.0

#define SCRATCH_MEM_SIZE 1048576  /* I hope that 1MB is enough */

extern xge_3Dwind   cwind;
extern xge_KnotWind kwind;

extern boolean curve;
extern boolean lamana;
extern boolean bezpoly;
extern boolean ticks;
extern boolean convh;
extern boolean nurbs;
extern boolean polarf;
extern boolean curvgr;
extern boolean torsgr;
extern boolean uniform;
extern int curv_graph_dens;
extern double curvgraphscale;
extern double torsgraphscale;

extern point4d cpoints[];
extern double  knots[];
extern int     npoints;
extern byte    cpmark[];

/* ////////////////////////////////////////////////////////////////////////// */
void ProjectCurve ( int id );
void ResetObject ( void );
void ResizeObject ( void );
boolean FindNearestPoint ( int id, int x, int y, int mindist );
void SetCPoint ( int id, int x, int y );
void FindBoundingBox ( Box3d *box );

void SelectCPoints ( int id, short x0, short x1, short y0, short y1 );
void UnSelectCPoints ( int id, short x0, short x1, short y0, short y1 );
void TransformCPoints ( trans3d *tr );
void SaveCPoints ( void );
void ClearPointMarking ( void );

boolean DegreeElevation ( void );
boolean DegreeReduction ( void );
void InsertKnot ( void );
void RemoveKnot ( void );
void SetUniformKnots ( void );

void UstawNURBS ( void );
void UstawPolarf ( void );
void UstawCurvGraph ( void );
void UstawTorsGraph ( void );
boolean UstawZamknieta ( void );
boolean RefineUniform ( void );

void DisplayCurve ( int id );
void DisplayCPolygon ( int id );
void DisplayCPoints ( int id );
void DisplayBezPoly ( int id );
void DisplayAuxPoints ( int id );
void DisplayConvh ( int id );
void DisplayKnots ( int y, int dy );
void DisplayKnotCursorPos ( int x, int y, int h );
void DisplayPolarf ( int id );
void DisplayCurvGraph ( int id );
