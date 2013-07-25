
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

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

int g2mbl_NiSize ( int nkn, int k )
{
  if ( k == 4 )
    return 16*9*nkn*nkn;
  else
    return (6*k+1)*9*nkn*nkn;
} /*g2mbl_NiSize*/

int g2mbl_NijSize ( int nkn, int k )
{
  if ( k == 4 )
    return 136*9*nkn*nkn;
  else {
    k = 6*k+1;
    return (((k+1)*k)/2)*9*nkn*nkn;
  }
} /*g2mbl_NijSize*/

int g2mbl_MijSize ( int nkn, int k )
{
  if ( k == 4 )
    return 120*18*nkn*nkn;
  else {
    k = 6*k+1;
    return (((k-1)*k)/2)*18*nkn*nkn;
  }
} /*g2mbl_MijSize*/

