
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005,2007                             */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void WriteKnots ( int hole_k, const float knots[][11] );
void WriteDomCP ( int hole_np, const point2f *domcp );
void WriteSurfCP ( int hole_np, const point3f *surfcp );
void WriteCamera ( const CameraRecf *CPos );
void WriteLightDir ( const vector3f *ld );
void WriteBezPatch ( int n, int m, const float *p );
void WriteBSPatch ( int n, int lknu, const float *knotsu,
                    int m, int lknv, const float *knotsv, const point3f *cp );
void WriteSurfPatches ( int hole_k, float knots[][11], const point3f *surfcp );
void OpenDataFile ( const char *fn );
void CloseDataFile ( void );

