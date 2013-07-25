
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_DEGREE   9
#define MAX_KNOTS  100

#define SCRATCH_MEM_SIZE 1048576  /* I hope that 1MB is enough */

#define MIN_KNOT_SCALE  0.01
#define MAX_KNOT_SCALE 100.0

extern boolean curve;
extern boolean lamana;
extern boolean bezpoly;
extern boolean funkcja;
extern boolean uniform;
extern boolean ticks;
extern boolean convh;
extern boolean nurbs;
extern boolean closed;
extern boolean baza;
extern boolean polarf;
extern boolean curvgr;
extern int curv_graph_dens;  /* curvature graph density */

extern xge_2Dwind   cwind;
extern xge_KnotWind kwind;

extern point3d cpoints[];
extern double  knots[];
extern int     npoints;
extern double  curvgraphscale;
extern byte    cpmark[];

void ProjectCurve ( void );

void FindBoundingBox ( Box2d *box );
void ResetObject ( void );
void ResizeObject ( void );
void ClearPointMarking ( void );
boolean FindNearestPoint ( int x, int y, int mindist );
void SetCPoint ( int x, int y );

void SelectCPoints ( short x0, short x1, short y0, short y1 );
void UnSelectCPoints ( short x0, short x1, short y0, short y1 );
void TransformCPoints ( trans2d *tr );
void SaveCPoints ( void );

boolean DegreeElevation ( void );
boolean DegreeReduction ( void );
void InsertKnot ( void );
void RemoveKnot ( void );
void SetUniformKnots ( void );
boolean RefineUniform ( void );

void UstawFunkcje ( void );
void UstawNURBS ( void );
void UstawPolarf ( void );
void UstawCurvGraph ( void );
boolean UstawZamknieta ( void );

void DisplayCurve ( void );
void DisplayCPolygon ( void );
void DisplayCPoints ( void );
void DisplayBezPoly ( void );
void DisplayAuxPoints ( void );
void DisplayKnots ( int y, int dy );
void DisplayKnotCursorPos ( int x, int y, int h );
void DisplayConvh ( void );
void DisplayBasis ( int y, int yu );
void DisplayPolarf ( void );
void DisplayCurvGraph ( void );

void DumpData ( void );
void ExportPovRay ( void );

