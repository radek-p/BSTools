
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g2blendingd.h"

#include "g2blprivated.h"

#define _DEBUG

/* This procedure sets up the constraint equations for a bicubic blending  */
/* patch, given by a set of curves of constant "u" parameter. The system   */
/* is in the form consistent with the optimization method implemented by   */
/* the procedure g2bl_FindBlSurfaceConstrLMTd, i.e. the unknown variables  */
/* are the coordinates of the internal control points and the right hand   */
/* side matrix has one column.                                             */

boolean g2bl_SetupUNLConstraintsd ( int lastknotu, int lastknotv,
                                    int ppitch, point3d *cp,
                                    int nucurv, double *ucknots,
                                    int cpitch, point3d *uccp,
                                    int *nconstr, double *cmat, double *crhs )
{
  void   *sp;
  int    nvars, ncon, bls, i, j;
  double *acmat;

  sp = pkv_GetScratchMemTop ();
  if ( lastknotu <= 9 || lastknotv <= 9 || nucurv < 1 )
    goto failure;
  bls = lastknotv-9;
  nvars = (lastknotu-9)*bls;
  ncon = nucurv*bls;
  acmat = pkv_GetScratchMemd ( nvars*ncon );
  if ( !acmat )
    goto failure;
  if ( !g2bl_SetupULConstraintsd ( lastknotu, lastknotv, 3,
                                   ppitch, (double*)cp,
                                   nucurv, ucknots, cpitch, (double*)uccp,
                                   &ncon, acmat, crhs ) )
    goto failure;

  *nconstr = 3*ncon;
  memset ( cmat, 0, 9*nvars*ncon*sizeof(double) );
  for ( i = 0; i < ncon; i++ )
    for ( j = 0; j < 3; j++ )
      pkv_Selectd ( nvars, 1, 1, 3, &acmat[i*nvars], &cmat[(3*i+j)*3*nvars+j] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_SetupUNLConstraintsd*/

