
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

float ComputeLaplacian ( int n, int m, const point2f *xycp, const float *zcp,
                         float u, float v );

void ComputeLaplacianGrad ( int n, int m,
                            const point2f *xycp, const float *zcp,
                            float u, float v, float *jac, vector2f *lgr );

void TestTabLaplacianU ( int n, int m, const point3f *cp, float u, int dd,  
                         point3f *gp );
void TestTabLaplacianV ( int n, int m, const point3f *cp, float v, int dd,  
                         point3f *gp );
void TestTabLapGradU ( int n, int m, const point3f *cp, float u, int dd,
                       char sxy, point3f *gp );
void TestTabLapGradV ( int n, int m, const point3f *cp, float v, int dd,
                       char sxy, point3f *gp );

