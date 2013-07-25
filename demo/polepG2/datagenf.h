
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

extern int hole_np;

#define MAX_BPTS (12*GH_MAX_K+2)

void InitKnots ( int k, float *knots );
int InitHole ( int k, float *ip, point3f *surfcp );
int InitDomain ( int k, vector2f *part, float *ip, point2f *domcp );
void StoreDomain ( int k, point2f *domcp );

