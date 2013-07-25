
extern boolean upsd;

void DrawBSNet ( int k, point3f *Bpt, CameraRecf *CPos, int i, boolean last );
void DrawBSPatch ( int k, float *knots, point3f *Bpt, CameraRecf *CPos, int i,
                   boolean first, boolean last );

void DrawBSDomNet ( int k, point2f *Dompt, CameraRecf *PPos, int i, boolean last );
void DrawBSDomNetNum ( int k, point2f *Dompt, CameraRecf *PPos );
void DrawBSDomPatch ( int k, float *knots, point2f *Dompt, CameraRecf *PPos, int i,
                      boolean first, boolean last );

void DrawBezCurve2a ( CameraRecf *PPos, int n, point2f *cp,
                      float t0, float t1 );
void DrawBezPatch2a ( CameraRecf *PPos, int n, int m, const point2f *cp,
                      float u0, float u1, float v0, float v1 );

void DrawBezCurve3a ( CameraRecf *CPos, int n, point3f *cp,
                      float t0, float t1 );
void DrawBezPatch3a ( CameraRecf *CPos, int n, int m, const point3f *cp,
                      float u0, float u1, float v0, float v1 );

