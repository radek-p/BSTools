
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"


/* /////////////////////////////////////////// */
/* matrix addition and subtraction             */
/* for explicit representation                 */

void pkn_SubtractMatrixd ( int nrows, int rowlen,
                           int inpitch1, const double *indata1,
                           int inpitch2, const double *indata2,
                           int outpitch, double *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = indata1[j]-indata2[j];
} /*pkn_SubtractMatrixd*/

void pkn_AddMatrixd ( int nrows, int rowlen,
                      int inpitch1, const double *indata1,
                      int inpitch2, const double *indata2,
                      int outpitch, double *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = indata1[j]+indata2[j];
} /*pkn_AddMatrixd*/

void pkn_AddMatrixMd ( int nrows, int rowlen,
                       int inpitch1, const double *indata1,
                       int inpitch2, const double *indata2,
                       double a,
                       int outpitch, double *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = indata1[j]+a*indata2[j];
} /*pkn_AddMatrixMd*/

void pkn_MatrixMDifferenced ( int nrows, int rowlen,
                              int inpitch1, const double *indata1,
                              int inpitch2, const double *indata2,
                              double a,
                              int outpitch, double *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = a*(indata1[j]-indata2[j]);
} /*pkn_MatrixMDifferenced*/

void pkn_MatrixLinCombd ( int nrows, int rowlen,
                          int inpitch1, const double *indata1,
                          double a,
                          int inpitch2, const double *indata2,
                          double b,
                          int outpitch, double *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = a*indata1[j]+b*indata2[j];
} /*pkn_MatrixLinCombd*/

