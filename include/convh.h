
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libpkgeom library of C procedures - */
/* finding convex hull                                     */

#ifndef CONVH_H
#define CONVH_H

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

void FindConvexHull2f ( int *n, point2f *p );
void FindConvexHull2d ( int *n, point2d *p );

#ifdef __cplusplus
}
#endif

#endif /*CONVH_H*/

