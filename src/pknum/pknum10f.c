
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

void pkn_SubtractMatrixf ( int nrows, int rowlen,
                           int inpitch1, const float *indata1,
                           int inpitch2, const float *indata2,
                           int outpitch, float *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = indata1[j]-indata2[j];
} /*pkn_SubtractMatrixf*/

void pkn_AddMatrixf ( int nrows, int rowlen,
                      int inpitch1, const float *indata1,
                      int inpitch2, const float *indata2,
                      int outpitch, float *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = indata1[j]+indata2[j];
} /*pkn_AddMatrixf*/

void pkn_AddMatrixMf ( int nrows, int rowlen,
                       int inpitch1, const float *indata1,
                       int inpitch2, const float *indata2,
                       double a,
                       int outpitch, float *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = (float)(indata1[j]+a*indata2[j]);
} /*pkn_AddMatrixMf*/

void pkn_MatrixMDifferencef ( int nrows, int rowlen,
                              int inpitch1, const float *indata1,
                              int inpitch2, const float *indata2,
                              double a,
                              int outpitch, float *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = (float)(a*(indata1[j]-indata2[j]));
} /*pkn_MatrixMDifferencef*/

void pkn_MatrixLinCombf ( int nrows, int rowlen,
                          int inpitch1, const float *indata1,
                          double a,
                          int inpitch2, const float *indata2,
                          double b,
                          int outpitch, float *outdata )
{
  int i, j;

  for ( i = 0;
        i < nrows;
        i++, indata1 += inpitch1, indata2 += inpitch2, outdata += outpitch )
    for ( j = 0; j < rowlen; j++ )
      outdata[j] = (float)(a*indata1[j]+b*indata2[j]);
} /*pkn_MatrixLinCombf*/

