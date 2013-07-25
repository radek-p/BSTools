
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#undef CONST_
#define CONST_

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

#define TRACE

#include "eg1hcore1d.c"
#include "eg1hcore2d.c"
#include "eg1hcore3d.c"
#include "eg1hcore4d.c"
#include "eg1h9core4d.c"
#include "eg1h9core5d.c"
#include "eg1hcore6d.c"
#include "eg1hcore7d.c"
#include "eg1hcore8d.c"
#include "eg1hcore9d.c"

