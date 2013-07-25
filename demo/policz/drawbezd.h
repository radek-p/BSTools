
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

void Get2dClipLines ( CameraRecd *CPos, vector3d *cliplines );

void DrawBezPatch2d ( CameraRecd *CPos,
                      int n, int m, const point2d *cp,
                      int d0, int d1, xgecolour_int c0, xgecolour_int c1 );
void DrawBezPatch3d ( CameraRecd *CPos,
                      int n, int m, const point3d *cp,
                      int d0, int d1, xgecolour_int c0, xgecolour_int c1 );
void DrawBSPatch3d ( CameraRecd *CPos,
                     int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point3d *cp,
                     int d0, int d1, xgecolour_int c0, xgecolour_int c1 );

