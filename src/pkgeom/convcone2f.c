
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"

short ExtendConvexCone2f ( ConvexCone2f *cone, float a )
{
  switch ( cone->code ) {
case 0:
    cone->amin = cone->amax = a;
    cone->code = 1;
    return 1;
case 1:
    if ( cone->amin <= cone->amax ) {
      if ( a > cone->amax ) {
        if ( a < cone->amin+2.0 ) cone->amax = a;
        else if ( a > cone->amax+2.0 ) cone->amin = a;
        else cone->code = 2;
      }
      else if ( a < cone->amin ) {
        if ( a >= cone->amin-2.0 ) {
          cone->amin = a;
          if ( cone->amax-cone->amin >= 2.0 )
            cone->code = 2;
        }
        else if ( a <= cone->amax+4.0 ) {
          cone->amax = a;
          if ( cone->amin-cone->amax <= 2.0 )
            cone->code = 2;
        }
        else cone->code = 2;
      }
    }
    else {
      if ( a > cone->amax && a < cone->amin-2.0 ) cone->amax = a;
      else if ( a >= cone->amin-2.0 && a <= cone->amax+2.0 ) cone->code = 2;
      else if ( a > cone->amax+2.0 && a < cone->amin ) cone->amin = a;
    }
default:
    return cone->code;
  }
} /*ExtendConvexCone2f*/

short ExtendConvexCone2fv ( ConvexCone2f *cone, vector2f *v )
{
  if ( v->x == 0.0 && v->y == 0.0 ) {
    cone->code = 2;
    return 2;
  }
  else
    return ExtendConvexCone2f ( cone, (float)pkv_SqAngle ( v->x, v->y ) );
} /*ExtendConvexCone2fv*/

boolean InsideConvexCone2f ( ConvexCone2f *cone, float a )
{
  if ( cone->amax >= cone->amin )
    return a >= cone->amin && a <= cone->amax;
  else
    return a <= cone->amax || a >= cone->amin;
} /*InsideConvexCone2f*/

