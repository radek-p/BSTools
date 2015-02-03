
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_PrincipalDirectionsBP3d ( int degreeu, int degreev,
                                      const point3d *ctlpoints, 
                                      double u, double v,
                                      double *k1, vector2d *v1,
                                      double *k2, vector2d *v2 )
{
  void    *sp;
  double  *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemd ( (36+3*(degreev+1))*3 );
  if ( !workspace )
    return false;
  result = _mbs_PrincipalDirectionsBP3d ( degreeu, degreev, ctlpoints, 
                                          u, v, k1, v1, k2, v2, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_PrincipalDirectionsBP3d*/

