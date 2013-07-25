
extern CameraRecf   CPos, PPos;

void DrawBezCurve2f ( int n, const point2f *cp, float t0, float t1 );
void DrawBezCurve3f ( int n, const point3f *cp, float t0, float t1 );

void DrawBSCurve2f ( int n, int lkn, const float *knots, const point2f *cp );
void DrawBSCurve3f ( int n, int lkn, const float *knots, const point3f *cp );

void SetLW ( int i, int i0, int i1 );

void DrawBezPatch2f ( int n, int m, const point2f *cp, int d0, int d1 );
void DrawBezPatch3f ( int n, int m, const point3f *cp, int d0, int d1 );
void DrawBezPatch3Rf ( int n, int m, const point3f *cp,  
                       float u0, float u1, float v0, float v1,  
                       int d0, int d1 );

void DrawBezPatchNet2f ( int n, int m, const point2f *cp );
void DrawBezPatchNet3f ( int n, int m, const point3f *cp );

void DrawBSPatch2f ( int n, int lknu, const float *knu,
                     int m, int lknv, const float *knv, const point2f *cp,
                     int d0, int d1 );
void DrawBSPatch3f ( int n, int lknu, const float *knu,
                     int m, int lknv, const float *knv, const point3f *cp,
                     int d0, int d1 );
void DrawBSPatch3Rf ( int n, int lknu, const float *knu,
                      int m, int lknv, const float *knv, const point3f *cp,
                      float u0, float u1, float v0, float v1,
                      int d0, int d1 );

float FindBSPatchMaxLaplacianf ( int n, int lknu, const float *knu,
                                 int m, int lknv, const float *knv,
                                 const point3f *cp, int d0, int d1 );
void DrawBSPatchLaplacianf ( int n, int lknu, const float *knu,
                             int m, int lknv, const float *knv,
                             const point3f *cp, int d0, int d1 );

float FindBezPatchMaxLaplacianf ( int n, int m,
                                  const point3f *cp, int d0, int d1 );
void DrawBezPatchLaplacianf ( int n, int m, const point3f *cp, int d0, int d1 );

