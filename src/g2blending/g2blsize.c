
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2010                            */
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

int g2bl_NiSize ( int nkn )
{
  return 16*9*nkn*nkn;
} /*g2bl_NiSize*/

int g2bl_NijSize ( int nkn )
{
  return 136*9*nkn*nkn;
} /*g2bl_NijSize*/

int g2bl_MijSize ( int nkn )
{
  return 120*18*nkn*nkn;
} /*g2bl_MijSize*/

