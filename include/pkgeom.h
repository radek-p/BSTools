
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libpkgeom library of C procedures - */
/* basic 2D, 3D and 4D affine geometry                     */ 

#ifndef PKGEOM_H
#define PKGEOM_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#include "pkgeomf.h"
#include "pkgeomd.h"

#ifdef __cplusplus   
extern "C" {
#endif

void TransPoint3df ( const trans3d *tr, const point3f *p, point3f *q );
void TransVector3df ( const trans3d *tr, const vector3f *v, vector3f *w );
void TransContra3df ( const trans3d *tri, const vector3f *v, vector3f *w );
void Trans3Point2df ( const trans3d *tr, const point2f *p, point2f *q );

void TransPoint3fd ( const trans3f *tr, const point3d *p, point3d *q );
void TransVector3fd ( const trans3f *tr, const vector3d *v, vector3d *w );
void TransContra3fd ( const trans3f *tri, const vector3d *v, vector3d *w );
void Trans3Point2fd ( const trans3f *tr, const point2d *p, point2d *q );

void Trans3fTod ( const trans3f *trf, trans3d *trd );
void Trans3dTof ( const trans3d *trd, trans3f *trf );

#ifdef __cplusplus
}
#endif

#endif /*PKGEOM_H*/

