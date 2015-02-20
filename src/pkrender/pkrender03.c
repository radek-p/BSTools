
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

boolean RendEnterCamerad ( pkRenderer *rend, CameraRecd *CPos )
{
  if ( CPos->parallel )
    return false;
      /* convert from double to single precision */
  rend->CPos = *CPos;
  if ( rend->CPos.width > rend->maxwidth ) rend->CPos.width = rend->maxwidth;
  if ( rend->CPos.height > rend->maxheight ) rend->CPos.height = rend->maxheight;
  CameraSetMappingd ( &rend->CPos );
  return true;
} /*RendEnterCamerad*/

void RendEnterLightsd ( pkRenderer *rend,
                        int nlights, const vector3d *lightdir,
                        const double *lightint )
{
  int   i;
  double lint;

  memset ( &rend->lightint[0], 0, (R_NLIGHTS+1)*sizeof(double) );
  for ( i = 0; i < nlights && i < R_NLIGHTS; i++ ) {
    rend->lightdir[i] = lightdir[i];
    NormalizeVector3d ( &rend->lightdir[i] );
    lint = min ( 1.0, lightint[i] );
    rend->lightint[i] = max ( 0.0, lint );
  }
  rend->lightint[R_NLIGHTS] = lightint[nlights];  /* ambient light intensity */
} /*RendEnterLightsd*/

void RendEnterReflectionLinesFramed ( pkRenderer *rend, point3d rf[3] )
{
  rend->rp0 = rf[0];
  SubtractPoints3d ( &rf[1], &rf[0], &rend->rv1 );
  SubtractPoints3d ( &rf[2], &rf[0], &rend->rv2 );
} /*RendEnterReflectionLinesFramed*/

void RendEnterHighlightLinesFramed ( pkRenderer *rend, point3d hf[3] )
{
  rend->hp0 = hf[0];
  SubtractPoints3d ( &hf[1], &hf[0], &rend->hv1 );
  SubtractPoints3d ( &hf[2], &hf[0], &rend->hv2 );
} /*RendEnterHighlightLinesFramed*/

void RendEnterSectionPlanesNormald ( pkRenderer *rend, vector3d *spn )
{
  rend->sectiondir = *spn;
} /*RendEnterSectionPlanesNormald*/

