
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
/* Way of using: Call CameraInitFramef, and then CameraInitPosf in oder to   */
/* initialise; then you may move/rotate the camera and project points etc.   */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef CAMERAF_H
#define CAMERAF_H

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


typedef struct CameraRecf {
      /* Data initialised before calling the CameraSetMappingf procedure */
  boolean  parallel;        /* if false, a perspective projection.  */
  boolean  upside;          /* if true, upside down.                */
  boolean  c_fixed;         /* if true, then centre is fixed in the */
                            /* camera coordinate system.            */
  byte     magnification;   /* the picture has to be magnified      */
                            /* if there is supersampling.           */
  short    xmin, ymin;      /* pixel margins of the frame.          */
  short    width, height;   /* pixel dimensions of the frame.       */
  float    aspect;          /* aspect factor of the screen.         */

  point3f  position;        /* objective position.                  */
  float    psi, theta, phi; /* Euler angles of the camera.          */
  point3f  g_centre;        /* centre of camera rotations in global */
  point3f  c_centre;        /* and scaled camera coordinates.       */
  float    zmin, zmax;      /* depth range.                  */
  union {
    struct {
      float f;              /* focal length.                 */
      float xi0, eta0;      /* shift after projection.       */
      float dxi0, deta0;    /* additional shifts in x and y. */
    } persp;
    struct {
      float wdt, hgh, diag; /* dimensions of the frame in space */
      char  dim_case;
    } para;
  } vd;                     /* variant data for perspective  */
                            /* and parallel projections      */
      /* Data computed from the above */
  float    xscale, yscale;  /* scaling of axes between the   */
                            /* camera and image coordinates. */
  trans3f  CTr;             /* affine transformation to the  */
                            /* scaled camera coordinates.    */
  trans3f  CTrInv;          /* inversion of CTr.             */
  char     ncplanes;        /* number of clipping planes.    */
  vector4f cplane[6];       /* clipping planes               */
                            /* currently only 4 of them are used */
} CameraRecf;


void CameraInitFramef ( CameraRecf *CPos,
                        boolean parallel, boolean upside,
                        short width, short height, short xmin, short ymin,
                        float aspect, int ncplanes );
void CameraSetMagf ( CameraRecf *CPos, byte mag );
void CameraSetDepthRangef ( CameraRecf *CPos, float zmin, float zmax );
boolean CameraSetMappingf ( CameraRecf *CPos );

void CameraProjectPoint3f ( CameraRecf *CPos, const point3f *p, point3f *q );
void CameraUnProjectPoint3f ( CameraRecf *CPos, const point3f *p, point3f *q );
void CameraProjectPoint2f ( CameraRecf *CPos, const point2f *p, point2f *q );
void CameraUnProjectPoint2f ( CameraRecf *CPos, const point2f *p, point2f *q );
void CameraProjectPoint3Rf ( CameraRecf *CPos, const point4f *p, point3f *q );
void CameraUnProjectPoint3Rf ( CameraRecf *CPos, const point3f *p, float w,
                               point4f *q );
void CameraProjectPoint2Rf ( CameraRecf *CPos, const point3f *p, point2f *q );
void CameraUnProjectPoint2Rf ( CameraRecf *CPos, const point2f *p, float w,
                               point3f *q );
void CameraProjectPoint3f2s ( CameraRecf *CPos, const point3f *p, point2s *q );
void CameraProjectPoint3Rf2s ( CameraRecf *CPos, const point4f *p, point2s *q );

void CameraRayOfPixelf ( CameraRecf *CPos, float xi, float eta, ray3f *ray );

void CameraInitPosf ( CameraRecf *CPos );
void CameraSetRotCentref ( CameraRecf *CPos, point3f *centre,
                           boolean global_coord, boolean global_fixed );
void CameraMoveToGf ( CameraRecf *CPos, point3f *pos );
void CameraTurnGf ( CameraRecf *CPos, float _psi, float _theta, float _phi );
void CameraMoveGf ( CameraRecf *CPos, vector3f *v );
void CameraMoveCf ( CameraRecf *CPos, vector3f *v );
void CameraRotGf ( CameraRecf *CPos, float _psi, float _theta, float _phi );
#define CameraRotXGf(CPos,angle) \
  CameraRotGf(CPos, 0.0, (float)(angle), 0.0)
#define CameraRotYGf(CPos,angle) \
  CameraRotGf(CPos, (float)(0.5 * PI), (float)(angle), (float)(-0.5 * PI))
#define CameraRotZGf(CPos,angle) \
  CameraRotGf(CPos, (float)(angle), 0.0, 0.0)
void CameraRotVGf ( CameraRecf *CPos, vector3f *v, float angle );
void CameraRotCf ( CameraRecf *CPos, float, float, float );
#define CameraRotXCf(CPos,angle) \
  CameraRotCf ( CPos, 0.0, (float)(angle), 0.0 )
#define CameraRotYCf(CPos,angle) \
  CameraRotCf ( CPos, (float)(0.5 * PI), (float)(angle), (float)(-0.5 * PI) )
#define CameraRotZCf(CPos,angle) \
  CameraRotCf ( CPos, (float)(angle), 0.0, 0.0 )
void CameraRotVCf ( CameraRecf *CPos, vector3f *v, float angle );
void CameraSetFf ( CameraRecf *CPos, float f );
void CameraZoomf ( CameraRecf *CPos, float fchange );


boolean CameraClipPoint3f ( CameraRecf *CPos, point3f *p, point3f *q );

boolean CameraClipLine3f ( CameraRecf *CPos,
                           point3f *p0, float t0, point3f *p1, float t1,
                           point3f *q0, point3f *q1 );

boolean CameraClipPolygon3f ( CameraRecf *CPos, int n, const point3f *p,
                              void (*output)(int n, point3f *p) );

/* the following procedures are intended to be changed later; */
/* their use should be avoided  */

boolean Pixel2Spacef ( CameraRecf *CPos, float x, float y, vector3f *n,
                       point3f *p, point3f *q,
                       boolean global_in, boolean global_out );

typedef struct PixHLine_detf {
  float D, Dt, Dz, Ds, Dts;
} PixHLine_detf;

boolean PixHLine2Spacef ( CameraRecf *CPos, int xi1, int xi2,
                          int eta, vector3f *n, point3f *p, point3f *p1,
                          point3f *p2, PixHLine_detf *det );

#ifdef __cplusplus
}
#endif

#endif /*CAMERAF_H*/

