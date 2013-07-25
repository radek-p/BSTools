
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

double ComputeLaplacian ( int n, int m, const point2d *xycp, const double *zcp,
                         double u, double v );

void ComputeLaplacianGrad ( int n, int m,
                            const point2d *xycp, const double *zcp,
                            double u, double v, double *jac, vector2d *lgr );

void TestTabLaplacianU ( int n, int m, const point3d *cp, double u, int dd,  
                         point3d *gp );
void TestTabLaplacianV ( int n, int m, const point3d *cp, double v, int dd,  
                         point3d *gp );
void TestTabLapGradU ( int n, int m, const point3d *cp, double u, int dd,
                       char sxy, point3d *gp );
void TestTabLapGradV ( int n, int m, const point3d *cp, double v, int dd,
                       char sxy, point3d *gp );

