
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "raybez.h"
#include "pkrender.h"

#include "pkrenderprivate.h"

void GetTexColour ( pkRenderer *rend,
                    ray3d *ray, renderobj *obj, RayObjectIntersd *roint,
                    vector3d *texcolour )
{
  double t;
  int    id;

  if ( rend->c_shape_func == shapefunc_NONE ) {
    id = roint->object_id;
    memcpy ( texcolour, rend->obj_tab[id].triang.colour, sizeof(vector3d) );
  }
  else {
    t = cShapeFunc ( rend, ray, obj, roint );
    t = (t-rend->minsf)/(rend->maxsf-rend->minsf);
    t = max ( 0.0, t );
    t = min ( 1.0, t );
  /* "rainbow" texture */
    if      ( t < 0.2 ) SetPoint3d ( texcolour, 1.0, 0.0, 1.0-5.0*t );  
    else if ( t < 0.4 ) SetPoint3d ( texcolour, 1.0, 5.0*(t-0.2), 0.0 );
    else if ( t < 0.6 ) SetPoint3d ( texcolour, 1.0-5.0*(t-0.4), 1.0, 0.0 );
    else if ( t < 0.8 ) SetPoint3d ( texcolour, 0.0, 1.0, 5.0*(t-0.6) );
    else                SetPoint3d ( texcolour, 0.0, 1.0-5.0*(t-0.8), 1.0 );
  }
  if ( !dShapeFunc ( rend, ray, obj, roint ) )
    MultVector3d ( 0.85, texcolour, texcolour );
} /*GetTexColour*/

void AmbientLight ( vector3d *tex, double amb_intens, vector3d *rgb )
{
  MultVector3d ( amb_intens, tex, rgb );
} /*AmbientLight*/

void PhongLight ( pkRenderer *rend,
                  vector3d *tex, point3d *p, vector3d *nv,
                  vector3d *vv, vector3d *light, double light_intens,
                  vector3d *rgb )
{
  double    c, d;
  vector3d h;

  NormalizeVector3d ( nv );
  c = DotProduct3d ( nv, light );
  d = DotProduct3d ( nv, vv );
  if ( c*d <= 0.0 ) {  
    if ( rend->swShadows ) {
      if ( IsInShadow ( rend, p, light ) )
        return;
    }
    SubtractPoints3d ( vv, light, &h );
    NormalizeVector3d ( &h );
    d = fabs ( DotProduct3d ( nv, &h ) );
    d = d*d;  d = d*d;  d = d*d;  d = d*d;  d = d*d;  d = d*d;
    c = 0.8*light_intens*fabs(c);
    d *= 0.15*light_intens;
    rgb->x += (c*tex->x + d);
    rgb->y += (c*tex->y + d);
    rgb->z += (c*tex->z + d);
  }
} /*PhongLight*/

void GetPixelColour ( pkRenderer *rend, ray3d *ray,
                      unsigned int *r, unsigned int *g, unsigned int *b )
{
  int              k;
  RayObjectIntersd iint;
  vector3d         texcolour, rgb;
  unsigned int     a;

  *r = *g = *b = 220;
  if ( FindRayTrInters ( rend, ray, &iint, &k ) ) {
    GetTexColour ( rend, ray, &rend->obj_tab[k], &iint, &texcolour );
    AmbientLight ( &texcolour, rend->lightint[R_NLIGHTS], &rgb );
    for ( k = 0; k < R_NLIGHTS; k++ )
      if ( rend->lightint[k] > 0.001 )
        PhongLight ( rend, &texcolour, &iint.p, &iint.nv, &ray->v,
                     &rend->lightdir[k], rend->lightint[k], &rgb );
    a = (unsigned int)(rgb.x*255.0);  *r = min ( a, 255 );
    a = (unsigned int)(rgb.y*255.0);  *g = min ( a, 255 );
    a = (unsigned int)(rgb.z*255.0);  *b = min ( a, 255 );
  }
} /*GetPixelColour*/

