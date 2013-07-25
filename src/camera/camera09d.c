
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "camera.h"

#include "msgpool.h"


boolean Pixel2Spaced ( CameraRecd *CPos, double x, double y, vector3d *n,
                       point3d *p, point3d *q,
                       boolean global_in, boolean global_out )
{
  /* this function finds the point q on the plane specified by its normal n and  */
  /* point p, that is transformed by the camera to the point (x,y) on the screen.*/
  /* if no such point or more than one exists, the function returns false.       */
  /* global_in specifies the coordinate system used for n and p (global if true, */
  /* scaled camera if false); global_out specifies the coord. system used for q. */
#define tol             1.0e-6
  double xi, eta, d;
  vector3d nn, pp;

  if (global_in) {
    TransContra3d ( &CPos->CTrInv, n, &nn ); /* n is a contravariant vector */
    TransPoint3d ( &CPos->CTr, p, &pp );
  } else {
    nn = *n;
    pp = *p;
  }
  xi = x - CPos->vd.persp.xi0;
  eta = y - CPos->vd.persp.eta0;
  d = nn.z + xi * nn.x + eta * nn.y;
  if (fabs(d) < tol)
    return false;
  q->z = DotProduct3d ( &nn, &pp ) / d;
  q->y = eta * q->z;
  q->x = xi * q->z;
  if (global_out)
    TransPoint3d ( &CPos->CTrInv, q, q );
  return true;
#undef tol
} /*Pixel2Spaced*/

boolean PixHLine2Spaced ( CameraRecd *CPos, int xi1, int xi2, int eta,
                          vector3d *n, point3d *p, point3d *p1, point3d *p2,
                          PixHLine_detd *det )
{
  /* horizontal line */
  /* plane */
  /* This procedure helps to find a series of points on a plane given by */
  /* n and p, that project on subsequent pixels on a horizontal line     */
  /* between (xi1,y) and (xi2,y). It returns points p1 and p2 that       */
  /* project on the line end points, in global coordinate system, and    */
  /* determinants, which allow to calculate the parameter t (for linear  */
  /* interpolation between p1 and p2) and z - depth of the point in the  */
  /* image/camera coordinate system.                                     */
  /*--------------------
    An example of use: procedure invoked by a scanline algorithm in order
    to fill a horizontal line on a facet image, with z-buffer visibility
    algorithm and texture.

    procedure FillHorizontalLine ( f : facet index;
                                   x1, x2, y : integer );
      var
        p1, p2, p : point3d;
        v         : vector3d;
        det       : PixHLine_detd;
        x         : integer;
        t, z      : double;
    begin {FillHorizontalLine}
      if not PixHLine2Space ( CPos, x1, x2-1, y, f.n, f.p, p1, p2, det )
        then exit;

      { here you can transform points p1, p2 from global to texture coordinates }

      with det do begin
        SubtractPoints3 ( p2, p1, v );
        for x := x1 to x2-1 do begin
          z := Dz/D;
          if ZBufPointVisible ( x, y, z )   { test visibility in the z-buffer }
          then begin
            t := Dt/D;
            AddVector3M ( p1, v, t, p );
            SetPixel ( x, y, Texture ( p ) )
          end;
                            { determinants for the next pixel; Dz is constant }
          D := D + Ds;  Dt := Dt + Ds
        end
      end
    end {FillHorizontalLine};
--------------------*/
  if (!Pixel2Spaced ( CPos, (double)xi1, (double)eta, n, p, p1, true, false ))
    return false;
  if (xi2 == xi1) {
    det->D = 1.0;
    det->Dt = 0.0;
    det->Dz = p1->z;
    det->Ds = 0.0;
    det->Dts = 0.0;
    *p2 = *p1;
  } else {
    if (!Pixel2Spaced ( CPos, (double)xi2, (double)eta, n, p, p2, true, false ))
      return false;
    det->Ds = p2->z - p1->z;
    det->Dts = -p1->z;
    det->D = p1->x - p2->x + (xi1 - CPos->vd.persp.xi0) * det->Ds;
    det->Dt = 0.0;
    det->Dz = p1->x * p2->z - p1->z * p2->x;
  }
  TransPoint3d ( &CPos->CTrInv, p1, p1 );
  TransPoint3d ( &CPos->CTrInv, p2, p2 );
  return true;
} /*PixHLine2Spaced*/

