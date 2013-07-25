
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

void mbs_SphericalProductd (
                int degree_eq, int lastknot_eq, const point2d *cpoints_eq,
                int degree_mer, int lastknot_mer, const point2d *cpoints_mer,
                int pitch, point3d *spr_cp )
{
  int i, j, k;

        /* the unit of the pitch is one double */
  pitch /= 3;
  for ( i = k = 0;  i < lastknot_eq-degree_eq;  i++, k += pitch )
    for ( j = 0; j < lastknot_mer-degree_mer; j++ )
      SetPoint3d ( &spr_cp[k+j],
                   cpoints_eq[i].x*cpoints_mer[j].x,
                   cpoints_eq[i].y*cpoints_mer[j].x,
                   cpoints_mer[j].y );
} /*mbs_SphericalProductd*/

