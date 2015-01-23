
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#define RBEZ_NEWTON_YES   0
#define RBEZ_NEWTON_NO    1
#define RBEZ_NEWTON_ERROR 2

boolean _rbez_ConvexHullTest1d ( int n, double *ff );
boolean _rbez_UniquenessTest1d ( int n, double *ff );
boolean _rbez_NewtonMethod1d ( int n, double *ff, double *z );

boolean _rbez_SubdividePatch2d ( int n, int m, point2d *p, point2d *q );
boolean _rbez_ConvexHullTest2d ( int ncp, point2d *mcp );
boolean _rbez_UniquenessTest2d ( int n, int m, int ncp, point2d *mcp,
                                 point2d *p, vector2d *du, vector2d *dv,
                                 double *K1, double *K2,
                                 double *workspace );
char _rbez_NewtonMethod2d ( int n, int m, point2d *mcp,
                            point2d *p, vector2d *pu, vector2d *pv,
                            point2d *z, double *workspace );
boolean _rbez_SecondTest2d ( point2d *z, int n, int m, double K1, double K2 );

