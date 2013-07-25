
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "msgpool.h"

void pkv_SortPermute2 ( unsigned int num_data,
                        size_t item_length1, void *data1,
                        size_t item_length2, void *data2,
                        unsigned int *permut )
{
  unsigned int  i, k, m;
  char *tmp1, *tmp2, *data1a, *data1b, *data2a, *data2b;

  tmp1 = (char*)pkv_GetScratchMem ( item_length1+item_length2 );
  if ( !tmp1 ) {
    pkv_SignalError ( LIB_PKVARIA, 2, ERRMSG_2 );
    exit ( 1 );
  }
  tmp2 = &tmp1[item_length1];

  for ( i = 0; i < num_data-1; i++ ) {
    if ( i != permut[i] ) {
      k = i;
      data1a = ((char*)data1) + item_length1*k;
      memcpy ( tmp1, data1a, item_length1 );
      data2a = ((char*)data2) + item_length2*k;
      memcpy ( tmp2, data2a, item_length2 );
      m = permut[k];
      while ( m != i ) {
        data1b = ((char*)data1) + item_length1*m;
        memcpy ( data1a, data1b, item_length1 );
        data2b = ((char*)data2) + item_length2*m;
        memcpy ( data2a, data2b, item_length2 );
        permut[k] = k;
        k = m;
        m = permut[k];
        data1a = data1b;
        data2a = data2b;
      }
      memcpy ( data1a, tmp1, item_length1 );
      memcpy ( data2a, tmp2, item_length2 );
      permut[k] = k;
    }
  }
  pkv_FreeScratchMem ( item_length1+item_length2 );
} /*pkv_SortPermute2*/

