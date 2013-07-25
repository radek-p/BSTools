
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */


extern xge_3Dwinf swind;

/* displaying switches */
extern boolean bsnet, bezp, auxcurv, starcurv, auxpatch, targpatch;

extern float beta1, beta2;
extern boolean autocpoint;

void FindBoundingBox ( Box3f *box );

void ResetObject ( void );
boolean FindNearestPoint ( int id, int x, int y, int mindist );
void SetCPoint ( int id, int x, int y );

void SelectCPoints ( int id, short x0, short x1, short y0, short y1 );
void UnSelectCPoints ( int id, short x0, short x1, short y0, short y1 );
void SaveCPoints ( void );
void TransformCPoints ( trans3f *tr );

void UpdateSurface ( void );
void ResizeObject ( void );
void RzutujXPunkt ( int id, point3f *p, XPoint* q );
void RzutujPK ( int id, const point3f *p, point3f *q );
void ProjectScene ( int id );

void DrawBSNets ( int id );
void DrawBezCurves ( int id );
void DrawBezPatches ( int id );
void DrawAuxCurves ( int id );
void DrawStarCurves ( int id );
void DrawAuxPatches ( int id );
void DrawTargetPatches ( int id );
void DrawSurface ( int id );

void UstawAutoCPoint ( void );
void UstawK3 ( void );
void UstawK5 ( void );
void UstawK6 ( void );
void UstawK8 ( void );

