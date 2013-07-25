
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/times.h>
#include <unistd.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"
#include "g2mblmlprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void g2mbl_MLGetTimes ( void *data,
                        float *time_prep, float *time_h, float *time_cg )
{
  mesh_ml_optdata *d;
  float           t;

  d = (mesh_ml_optdata*)data;
  t = (float)sysconf(_SC_CLK_TCK);
  *time_prep = (float)d->time_prep/t;
  *time_h    = (float)d->time_h/t;
  *time_cg   = (float)d->time_cg/t;
} /*g2mbl_MLGetTimes*/

