
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005,2007                             */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void WriteKnots ( int hole_k, const double knots[][11] );
void WriteDomCP ( int hole_np, const point2d *domcp );
void WriteSurfCP ( int hole_np, const point3d *surfcp );
void WriteCamera ( const CameraRecd *CPos );
void WriteLightDir ( const vector3d *ld );
void WriteBezPatch ( int n, int m, const double *p );
void WriteBSPatch ( int n, int lknu, const double *knotsu,
                    int m, int lknv, const double *knotsv, const point3d *cp );
void WriteSurfPatches ( int hole_k, double knots[][11], const point3d *surfcp );
void OpenDataFile ( const char *fn );
void CloseDataFile ( void );

