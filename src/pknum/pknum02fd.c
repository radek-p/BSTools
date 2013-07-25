
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

/* A part of the libpknum library of C procedures      */
/* Copyright (C) by Przemyslaw Kiciak, 2005            */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"


void pkn_PrintProfile ( int ncols, const bandm_profile *prof )
{
  int i;

  for ( i = 0; i < ncols; i++ ) {
    printf ( "%3d: %3d, %3d, %3d\n", i, prof[i].firstnz,
             prof[i].firstnz+prof[i+1].ind-prof[i].ind,
             prof[i].ind );
  }
  printf ( "%3d: %3d, %3d, %3d\n\n",
           i, prof[ncols].firstnz, 0, prof[ncols].ind );
} /*pkn_PrintProfile*/

