
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

extern int hole_np;

#define MAX_BPTS (12*GH_MAX_K+2)

extern point3d Bpt[MAX_BPTS];
extern point2d Dompt[MAX_BPTS];
extern double  knots[11*GH_MAX_K];

void InitKnots ( int kk );
void ModifyKnots ( int kk );
void InitHole ( int kk, double *ip );
void InitDomain ( int kk, double *ip );
