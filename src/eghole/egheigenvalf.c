
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "egholef.h"

/* below are eigenvalues of the mesh refinement operator, */
/* whose eigenvectors are the meshes given in the         */
/* egheigenmeshf.c source file. */

float egh_eigenvalf[GH_MAX_K-3] =
  {0.4100970508,  /* k == 3 */
   0.5499883545,  /* k == 5 */
   0.5796823261,  /* k == 6 */
   0.5985102835,  /* k == 7 */
   0.6111165267,  /* k == 8 */
   0.6199392206,  /* k == 9 */
   0.6263412675,  /* k == 10 */
   0.6311275867,  /* k == 11 */
   0.6347964124,  /* k == 12 */
   0.6376687290,  /* k == 13 */
   0.6399585303,  /* k == 14 */
   0.6418127555,  /* k == 15 */
   0.6433349225}; /* k == 16 */

