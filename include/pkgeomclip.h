
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef PKGEOMCLIP_H
#define PKGEOMCLIP_H

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define PKGEOM_CLIP_NONE   0
#define PKGEOM_CLIP_PART   1
#define PKGEOM_CLIP_ENTIRE 2

int LiangBarskyClip2f ( point2f *p0, point2f *p1, float t0, float t1,
                        Box2f *box,
                        point2f *q0, point2f *q1 );

int LiangBarskyClip2d ( point2d *p0, point2d *p1, double t0, double t1,
                        Box2d *box,
                        point2d *q0, point2d *q1 );

#ifdef __cplusplus
}
#endif

#endif

