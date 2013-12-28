
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libcamera library of C procedures - 3D camera */
/* for perspective and parallel projections                          */

/* ///////////////////////////////////////////////////////////////////////// */
/* This library implements perspective projections of 3D space to 2D plane.  */
/* and it is intended to be used in graphics. There are 4 sets of functions: */
/* 1. Camera setting in the space and finding the point image,               */
/* 2. Camera manipulations - changes of position and direction,              */
/* 3. Point, line and polygon clipping,                                      */
/* 4. Finding the point in the space on the given plane, that maps to the    */
/*    given pixel.                                                           */
/* ///////////////////////////////////////////////////////////////////////// */
/* There are 4 coordinate systems in use:                                    */
/* - global, in which geometrical objects are defined,                       */
/* - camera, in which the perspective projection is performed,               */
/* - scaled camera, shift by (xi0, eta0, 0) transforms to the image,         */
/* - image,  in which pixel coordinates are given.                           */
/* x and y axes of the image coordinate system, are scaled x and y axes of   */
/* the camera system; the z coordinates in both systems are identical.       */
/* -points and vectors defining objects are specified in global coordinates. */
/* -camera movements may be specified in global or camera coordinates.       */
/* -pixel coordinates are given in the image system.                         */
/* Any exceptions (extensions) to the above rules are commented in the code. */
/* ///////////////////////////////////////////////////////////////////////// */
/* Way of using: Call CameraInitFramed, and then CameraInitPosd in oder to   */
/* initialise; then you may move/rotate the camera and project points etc.   */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef CAMERAD_H
#define CAMERAD_H

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef CAMERA_H
#define CPLANE_TOP    0
#define CPLANE_BOTTOM 1
#define CPLANE_LEFT   2
#define CPLANE_RIGHT  3
#define CPLANE_NEAR   4
#define CPLANE_FAR    5
#endif


typedef struct CameraRecd {
      /* Data initialised before calling the CameraSetMappingd procedure */
  boolean  parallel;        /* if false, a perspective projection.  */
  boolean  upside;          /* if true, upside down.                */
  boolean  c_fixed;         /* if true, then centre is fixed in the */
                            /* camera coordinate system.            */
  byte     magnification;   /* the picture has to be magnified      */
                            /* if there is supersampling.           */
  short    xmin, ymin;      /* pixel margins of the frame.          */
  short    width, height;   /* pixel dimensions of the frame.       */
  double   aspect;          /* aspect factor of the screen.         */

  point3d  position;        /* objective position.                  */
  double   psi, theta, phi; /* Euler angles of the camera.          */
  point3d  g_centre;        /* centre of camera rotations in global */
  point3d  c_centre;        /* and scaled camera coordinates.       */
  double   zmin, zmax;      /* depth range.                  */
  union {
    struct {
      double f;             /* focal length.                 */
      double xi0, eta0;     /* shift after projection.       */
      double dxi0, deta0;   /* additional shifts in x and y. */
    } persp;
    struct {
      double wdt, hgh, diag; /* dimensions of the frame in space */
      char   dim_case;
    } para;
  } vd;                     /* variant data for perspective  */
                            /* and parallel projections      */
      /* Data computed from the above and the union below */
  double   xscale, yscale;  /* scaling of axes between the   */
                            /* camera and image coordinates. */
  trans3d  CTr;             /* affine transformation to the  */
                            /* scaled camera coordinates.    */
  trans3d  CTrInv;          /* inversion of CTr.             */
  char     ncplanes;        /* number of clipping planes.    */
  vector4d cplane[6];       /* clipping planes               */
                            /* currently only 4 of them are used */
} CameraRecd;


void CameraInitFramed ( CameraRecd *CPos,
                        boolean parallel, boolean upside,
                        short width, short height, short xmin, short ymin,
                        double aspect, int ncplanes );
void CameraSetMagd ( CameraRecd *CPos, byte mag );
void CameraSetDepthRanged ( CameraRecd *CPos, double zmin, double zmax );
boolean CameraSetMappingd ( CameraRecd *CPos );

void CameraProjectPoint3d ( CameraRecd *CPos, const point3d *p, point3d *q );
void CameraUnProjectPoint3d ( CameraRecd *CPos, const point3d *p, point3d *q );
void CameraProjectPoint2d ( CameraRecd *CPos, const point2d *p, point2d *q );
void CameraUnProjectPoint2d ( CameraRecd *CPos, const point2d *p, point2d *q );
void CameraProjectPoint3Rd ( CameraRecd *CPos, const point4d *p, point3d *q );
void CameraUnProjectPoint3Rd ( CameraRecd *CPos, const point3d *p, double w,
                               point4d *q );
void CameraProjectPoint2Rd ( CameraRecd *CPos, const point3d *p, point2d *q );
void CameraUnProjectPoint2Rd ( CameraRecd *CPos, const point2d *p, double w,
                               point3d *q );
void CameraProjectPoint3d2s ( CameraRecd *CPos, const point3d *p, point2s *q );
void CameraProjectPoint3Rd2s ( CameraRecd *CPos, const point4d *p, point2s *q );

void CameraRayOfPixeld ( CameraRecd *CPos, double xi, double eta, ray3d *ray );

void CameraInitPosd ( CameraRecd *CPos );
void CameraSetRotCentred ( CameraRecd *CPos, point3d *centre,
                           boolean global_coord, boolean global_fixed );
void CameraMoveToGd ( CameraRecd *CPos, point3d *pos );
void CameraTurnGd ( CameraRecd *CPos, double _psi, double _theta, double _phi );
void CameraMoveGd ( CameraRecd *CPos, vector3d *v );
void CameraMoveCd ( CameraRecd *CPos, vector3d *v );
void CameraRotGd ( CameraRecd *CPos, double _psi, double _theta, double _phi );
#define CameraRotXGd(CPos,angle) \
  CameraRotGd(CPos, 0.0, angle, 0.0)
#define CameraRotYGd(CPos,angle) \
  CameraRotGd(CPos, 0.5 * PI, angle, -0.5 * PI)
#define CameraRotZGd(CPos,angle) \
  CameraRotGd(CPos, angle, 0.0, 0.0)
void CameraRotVGd ( CameraRecd *CPos, vector3d *v, double angle );
void CameraRotCd ( CameraRecd *CPos, double, double, double );
#define CameraRotXCd(CPos,angle) \
  CameraRotCd ( CPos, 0.0, angle, 0.0 )
#define CameraRotYCd(CPos,angle) \
  CameraRotCd ( CPos, 0.5 * PI, angle, -0.5 * PI )
#define CameraRotZCd(CPos,angle) \
  CameraRotCd ( CPos, angle, 0.0, 0.0 )
void CameraRotVCd ( CameraRecd *CPos, vector3d *v, double angle );
void CameraSetFd ( CameraRecd *CPos, double f );
void CameraZoomd ( CameraRecd *CPos, double fchange );


boolean CameraClipPoint3d ( CameraRecd *CPos, point3d *p, point3d *q );

boolean CameraClipLine3d ( CameraRecd *CPos,
                           point3d *p0, double t0, point3d *p1, double t1,
                           point3d *q0, point3d *q1 );

boolean CameraClipPolygon3d ( CameraRecd *CPos, int n, const point3d *p,
                              void (*output)(int n, point3d *p) );

/* the following procedures are intended to be changed later; */
/* their use should be avoided  */

boolean Pixel2Spaced ( CameraRecd *CPos, double x, double y, vector3d *n,
                       point3d *p, point3d *q,
                       boolean global_in, boolean global_out );

typedef struct PixHLine_detd {
  double D, Dt, Dz, Ds, Dts;
} PixHLine_detd;

boolean PixHLine2Spaced ( CameraRecd *CPos, int xi1, int xi2,
                          int eta, vector3d *n, point3d *p, point3d *p1,
                          point3d *p2, PixHLine_detd *det );

#ifdef __cplusplus
}
#endif

#endif /*CAMERAD_H*/

