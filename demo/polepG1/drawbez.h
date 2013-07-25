
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void DrawBezPatch3a ( CameraRecf *CPos,
                      int n, int m, const point3f *p,
                      int dd1, int dd2, int colour1, int colour2 );
void DrawBezCurve3a ( CameraRecf *CPos,
                      int n, const point3f *p, int colour );
void DrawBezPatch3b ( CameraRecf *PPos,
                      int n, int m, const point3f *p,
                      int dd1, int dd2, int colour1, int colour2 );
void DrawBezCurve3b ( CameraRecf *PPos,
                      int n, const point3f *p, int colour );

void DrawBezPatch2 ( CameraRecf *PPos,
                     int n, int m, const point2f *p,
                     int dd1, int dd2, int colour1, int colour2 );
void DrawBezPatch2a ( CameraRecf *PPos,
                      int n, int m, const point2f *p,
                      float u0, float u1, float v0, float v1,
                      int dd1, int dd2, int colour1, int colour2 );

