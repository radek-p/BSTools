
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libcamera library of C procedures - 3D camera */
/* for perspective projections                                       */

#ifndef CAMERA_H
#define CAMERA_H

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#define CPLANE_TOP    0
#define CPLANE_BOTTOM 1
#define CPLANE_LEFT   2
#define CPLANE_RIGHT  3
#define CPLANE_NEAR   4
#define CPLANE_FAR    5

#include "cameraf.h"
#include "camerad.h"

#endif /*CAMERA_H*/

