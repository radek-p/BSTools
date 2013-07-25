
extern CameraRecd   CPos, PPos;

void DrawBezCurve2d ( int n, const point2d *cp, double t0, double t1 );
void DrawBezCurve3d ( int n, const point3d *cp, double t0, double t1 );

void DrawBSCurve2d ( int n, int lkn, const double *knots, const point2d *cp );
void DrawBSCurve3d ( int n, int lkn, const double *knots, const point3d *cp );

void SetLW ( int i, int i0, int i1 );

void DrawBezPatch2d ( int n, int m, const point2d *cp, int d0, int d1 );
void DrawBezPatch3d ( int n, int m, const point3d *cp, int d0, int d1 );
void DrawBezPatch3Rd ( int n, int m, const point3d *cp,  
                       double u0, double u1, double v0, double v1,  
                       int d0, int d1 );

void DrawBezPatchNet2d ( int n, int m, const point2d *cp );
void DrawBezPatchNet3d ( int n, int m, const point3d *cp );

void DrawBSPatch2d ( int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point2d *cp,
                     int d0, int d1 );
void DrawBSPatch3d ( int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point3d *cp,
                     int d0, int d1 );
void DrawBSPatch3Rd ( int n, int lknu, const double *knu,
                      int m, int lknv, const double *knv, const point3d *cp,
                      double u0, double u1, double v0, double v1,
                      int d0, int d1 );

double FindBSPatchMaxLaplaciand ( int n, int lknu, const double *knu,
                                 int m, int lknv, const double *knv,
                                 const point3d *cp, int d0, int d1 );
void DrawBSPatchLaplaciand ( int n, int lknu, const double *knu,
                             int m, int lknv, const double *knv,
                             const point3d *cp, int d0, int d1 );

