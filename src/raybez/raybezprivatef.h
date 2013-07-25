
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean _rbez_ConvexHullTest2f ( int ncp, point2f *mcp );
boolean _rbez_UniquenessTest2f ( int n, int m, int ncp, point2f *mcp,
                                 point2f *p, vector2f *du, vector2f *dv,
                                 float *K1, float *K2 );
boolean _rbez_NewtonMethod2f ( int n, int m, point2f *mcp,
                               point2f *p, vector2f *pu, vector2f *pv,
                               point2f *z );
boolean _rbez_SecondTest2f ( point2f *z, int n, int m, float K1, float K2 );

