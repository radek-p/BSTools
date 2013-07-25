
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

void pkn_SPMindex2to3 ( unsigned int nnz, index2 *ai, index3 *sai )
{
  int k;

  pkv_Selectc ( nnz, sizeof(index2), sizeof(index2), sizeof(index3),
                (char*)ai, (char*)sai );
  for ( k = 0; k < nnz; k++ )
    sai[k].k = k;
} /*pkn_SPMindex2to3*/

void pkn_SPMindex3to2 ( unsigned int nnz, index3 *sai, index2 *ai )
{
  pkv_Selectc ( nnz, sizeof(index2), sizeof(index3), sizeof(index2),
                (char*)sai, (char*)ai );
} /*pkn_SPMindex3to2*/

