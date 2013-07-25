
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

extern int hole_np;

#define MAX_BPTS (12*GH_MAX_K+2)

extern point3f Bpt[MAX_BPTS];
extern point2f Dompt[MAX_BPTS];
extern float   knots[11*GH_MAX_K];

void InitKnots ( int kk );
void ModifyKnots ( int kk );
void InitHole ( int kk, float *ip );
void InitDomain ( int kk, float *ip );
