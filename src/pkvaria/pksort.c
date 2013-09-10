
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <endian.h>
#if __BYTE_ORDER != __LITTLE_ENDIAN
#error The big endian version has not been tested yet!
/* if you have this error during the compilation, please comment out the   */
/* #error line above, compile and test the procedure, and then let me know */
/* the test result - I do not have a computer with a big endian processor  */
#endif

#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"

/* The code is written with assumption that sizeof(int) is 4 and sizeof(short) is 2 */

#define BOUND_VALUE  32
#define MAXCOUNT    256

static void _ConvertInts ( size_t key_size, size_t item_length,
                           int num_data, unsigned char *mdata )
{
  int i;

  mdata += key_size-1;
  for ( i = 0;  i < num_data;  i++, mdata += item_length )
    *mdata ^= 0x80;
} /*_ConvertInts*/

static void _ConvertFloats1 ( size_t item_length,
                              int num_data, unsigned char *mdata )
{
  int i;
  int *ap;
#define CVF1 { ap = (int*)mdata; \
    if ( *ap >= 0 ) *ap ^= 0x80000000; \
               else *ap = ~(*ap); \
  }

  for ( i = 0;  i < num_data;  i++, mdata += item_length )
    CVF1
#undef CVF1
} /*_ConvertFloats1*/

static void _ConvertFloats2 ( size_t item_length,
                              int num_data, unsigned char *mdata )
{
  int i;
  int *ap;
#define CVF2 { ap = (int*)mdata; \
    if ( *ap < 0 ) *ap ^= 0x80000000; \
               else *ap = ~(*ap); \
  }

  for ( i = 0;  i < num_data;  i++, mdata += item_length )
    CVF2
#undef CVF2
} /*_ConvertFloats2*/

static void _ConvertDoubles1 ( size_t item_length,
                              int num_data, unsigned char *mdata )
{
  int i;
  int *ap;
#define CVD1 { ap = (int*)mdata; \
    if ( ap[1] >= 0 ) ap[1] ^= 0x80000000; \
               else { ap[0] = ~ap[0];  ap[1] = ~ap[1]; } \
  }

  for ( i = 0;  i < num_data;  i++, mdata += item_length )
    CVD1
#undef CVD1
} /*_ConvertDoubles1*/

static void _ConvertDoubles2 ( size_t item_length,
                              int num_data, unsigned char *mdata )
{
  int i;
  int *ap;
#define CVD2 { ap = (int*)mdata; \
    if ( ap[1] < 0 ) ap[1] ^= 0x80000000; \
              else { ap[0] = ~ap[0];  ap[1] = ~ap[1]; } \
  }

  for ( i = 0;  i < num_data;  i++, mdata += item_length )
    CVD2
#undef CVD2
} /*_ConvertDoubles2*/

/* The trick with IEEE754 works, provided that the sequence */
/* does not contain unnormals (denormals are o.k.). */
/* The application is responsible for that. */
static void _ConvertLDoubles1 ( size_t item_length,
                               int num_data, unsigned char *mdata )
{
  int   i;
  short *ap;
  int   *bp;
#define CVLD1 { ap = (short*)mdata; \
    if ( ap[4] >= 0 ) ap[4] ^= 0x8000; \
               else { bp = (int*)mdata; \
                      bp[0] = ~bp[0];  bp[1] = ~bp[1];  ap[4] = ~ap[4]; } \
  }

  for ( i = 0;  i < num_data;  i++, mdata += item_length )
    CVLD1
#undef CVLD1
} /*_ConvertLDoubles1*/

static void _ConvertLDoubles2 ( size_t item_length,
                               int num_data, unsigned char *mdata )
{
  int   i;
  short *ap;
  int   *bp;
#define CVLD2 { ap = (short*)mdata; \
    if ( ap[4] < 0 ) ap[4] ^= 0x8000; \
               else { bp = (int*)mdata; \
                      bp[0] = ~bp[0];  bp[1] = ~bp[1];  ap[4] = ~ap[4]; } \
  }

  for ( i = 0;  i < num_data;  i++, mdata += item_length )
    CVLD2
#undef CVLD2
} /*_ConvertLDoubles2*/

static size_t _ConvertSigned1 ( size_t key_size, char key_type,
                                size_t item_length,
                                unsigned char *mdata, unsigned int num_data )
{
  switch ( key_type ) {
case ID_UNSIGNED:
    return key_size;
case ID_SIGNED_INT:
    _ConvertInts ( key_size, item_length, num_data, mdata );
    return key_size;
case ID_IEEE754_FLOAT:
    if ( key_size != sizeof(float) )
      return 0;
    _ConvertFloats1 ( item_length, num_data, mdata );
    return key_size;
case ID_IEEE754_DOUBLE:
    if ( key_size != sizeof(double) )
      return 0;
    _ConvertDoubles1 ( item_length, num_data, mdata );
    return key_size;
case ID_IEEE754_LONG_DOUBLE:
    if ( key_size != sizeof(long double) )
      return 0;
    _ConvertLDoubles1 ( item_length, num_data, mdata );
    return 10;  /* most significant bytes are unused */
default:
    return 0;
  }
} /*_ConvertSigned1*/

char pkv_SortKernel ( size_t key_size, char key_type,
                      size_t item_length, size_t key_offset,
                      unsigned int num_data, void *data,
                      unsigned int *permut )
{
  unsigned char *mdata;
  unsigned int  *count, *iptr0, *auxpermut;
  unsigned char *InData, *InDataBase;
  int           i, j, k, m;
  int           itlnumbyt;
  void          *top;
  int           countsize;
  boolean       odd_swap;

  if ( num_data <= 1 )
    return SORT_OK;

  mdata = (unsigned char*)data;
  mdata += key_offset;

  if ( num_data < BOUND_VALUE ) {   /* insertion sort of a short sequence */
        /* convert the sequence with signed numbers */
        /* the key size may change in case of long doubles */
    key_size = _ConvertSigned1 ( key_size, key_type, item_length, mdata, num_data );

    switch ( key_size ) {
case 0:
      return SORT_BAD_DATA;
case sizeof(short):
      {
        unsigned short *ptr0, *ptr1;

        for ( i = 1; i < num_data; i++ ) {
          for ( j = i, k = j-1;  j > 0;  j = k, k-- ) {
            ptr0 = (unsigned short*)(mdata + permut[j]*item_length);
            ptr1 = (unsigned short*)(mdata + permut[k]*item_length);
            if ( *ptr0 < *ptr1 )
              { m = permut[j];  permut[j] = permut[k];  permut[k] = m; }
            else
              break;
          }
        }
      }
      break;
case sizeof(int):
      {
        unsigned int *ptr0, *ptr1;

        for ( i = 1; i < num_data; i++ ) {
          for ( j = i, k = j-1;  j > 0;  j = k, k-- ) {
            ptr0 = (unsigned int*)(mdata + permut[j]*item_length);
            ptr1 = (unsigned int*)(mdata + permut[k]*item_length);
            if ( *ptr0 < *ptr1 )
              { m = permut[j];  permut[j] = permut[k];  permut[k] = m; }
            else
              break;
          }
        }
      }
      break;
case 2*sizeof(int):
      {
#if __WORDSIZE == 64
        unsigned long int *ptr0, *ptr1;

        for ( i = 1; i < num_data; i++ ) {
          for ( j = i, k = j-1;  j > 0;  j = k, k-- ) {
            ptr0 = (unsigned long int*)(mdata + permut[j]*item_length);
            ptr1 = (unsigned long int*)(mdata + permut[k]*item_length);
            if ( *ptr0 < *ptr1 )
              { m = permut[j];  permut[j] = permut[k];  permut[k] = m; }
            else
              break;
          }
        }
#else /* __WORDSIZE == 32 */
        unsigned int *ptr0, *ptr1;

        for ( i = 1; i < num_data; i++ ) {
          for ( j = i, k = j-1;  j > 0;  j = k, k-- ) {
            ptr0 = (unsigned int*)(mdata + permut[j]*item_length);
            ptr1 = (unsigned int*)(mdata + permut[k]*item_length);
#if __BYTE_ORDER == __LITTLE_ENDIAN
            if ( ptr0[1] < ptr1[1] || (ptr0[1] == ptr1[1] && ptr0[0] < ptr1[0]) )
#else
#if __BYTE_ORDER == __BIG_ENDIAN
            if ( ptr0[0] < ptr1[0] || (ptr0[0] == ptr1[0] && ptr0[1] < ptr1[1]) )
#else
#error Either little endian or big endian byte ordering is assumed
#endif
#endif
              { m = permut[j];  permut[j] = permut[k];  permut[k] = m; }
            else
              break;
          }
        }
#endif
      }
      break;
case 10:  /* most probably long double */
      {
#if __WORDSIZE == 64
        unsigned long int  *ptr0, *ptr1;
        unsigned short int *ptr0s, *ptr1s;

        for ( i = 1; i < num_data; i++ ) {
          for ( j = i, k = j-1;  j > 0;  j = k, k-- ) {
            ptr0 = (unsigned long int*)(mdata + permut[j]*item_length);
            ptr0s = (unsigned short int*)ptr0;
            ptr1 = (unsigned long int*)(mdata + permut[k]*item_length);
            ptr1s = (unsigned short int*)ptr1;
#if __BYTE_ORDER == __LITTLE_ENDIAN
            if ( ptr0s[4] < ptr1s[4] || (ptr0s[4] == ptr1s[4] && *ptr0 < *ptr1) )
#else
#if __BYTE_ORDER == __BIG_ENDIAN
            if ( *ptr0 < *ptr1 || (*ptr0 == *ptr1 && ptr0s[4] < ptr1s[4]) )
#else
#error Either little endian or big endian byte ordering is assumed
#endif
#endif
              { m = permut[j];  permut[j] = permut[k];  permut[k] = m; }
            else
              break;

          }
        }
#else /* __WORDSIZE == 32 */
        unsigned int *ptr0, *ptr1;
        unsigned short *ptr0s, *ptr1s;

        for ( i = 1; i < num_data; i++ ) {
          for ( j = i, k = j-1;  j > 0;  j = k, k-- ) {
            ptr0 = (unsigned int*)(mdata + permut[j]*item_length);
            ptr0s = (unsigned short*)ptr0;
            ptr1 = (unsigned int*)(mdata + permut[k]*item_length);
            ptr1s = (unsigned short*)ptr1;
#if __BYTE_ORDER == __LITTLE_ENDIAN
            if ( ptr0s[4] < ptr1s[4] )      goto swap1;
            else if ( ptr0s[4] > ptr1s[4] ) break;
            if ( ptr0[1] < ptr1[1] )        goto swap1;
            else if ( ptr0[1] > ptr1[1] )   break;
            if ( ptr0[0] >= ptr1[0] )       break;
#else
#if __BYTE_ORDER == __BIG_ENDIAN
            if ( ptr0[0] < ptr1[0] )      goto swap1;
            else if ( ptr0[0] > ptr1[0] ) break;
            if ( ptr0[1] < ptr1[1] )      goto swap1;
            else if ( ptr0[1] > ptr1[1] ) break;
            if ( ptr0s[4] >= ptr1s[4] )   break;
#else
#error Either little endian or big endian byte ordering is assumed
#endif
#endif
            else
              break;
swap1:
            { m = permut[j];  permut[j] = permut[k];  permut[k] = m; }
          }
        }
#endif
      }
      break;
default:
      {
/* if the key size is even or big, it is still possible to accelerate */
/* the comparison of the keys; by comparing them machine word-by word; */
/* now it is done byte-by-byte */
        unsigned char *ptr0, *ptr1;

        for ( i = 1; i < num_data; i++ ) {
          for ( j = i, k = j-1;  j > 0;  j = k, k-- ) {
            ptr0 = mdata + permut[j]*item_length;
            ptr1 = mdata + permut[k]*item_length;
#if __BYTE_ORDER == __LITTLE_ENDIAN
            for ( m = key_size-1; m >= 0; m-- )
#else
#if __BYTE_ORDER == __BIG_ENDIAN
            for ( m = 0; < key_size; m++ )
#else
#error Either little endian or big endian byte ordering is assumed
#endif
#endif
              if ( ptr0[m] < ptr1[m] )
                goto swap2;
              else if ( ptr0[m] > ptr1[m] )
                break;
            break;
swap2:
            m = permut[j];  permut[j] = permut[k];  permut[k] = m;
          }
        }
      }
      break;
    }

  }
  else {                            /* count sort */
        /* allocate the working arrays */
    top = pkv_GetScratchMemTop ();
    auxpermut = pkv_GetScratchMem ( num_data*sizeof(int) );
    count = pkv_GetScratchMem ( countsize = MAXCOUNT*key_size*sizeof(int) );
    if ( !auxpermut || !count )
      return SORT_NO_MEMORY;
    
        /* convert the sequence with signed numbers */
        /* the key size may change in case of long doubles */
    key_size = _ConvertSigned1 ( key_size, key_type, item_length, mdata, num_data );
    if ( !key_size ) {
      pkv_SetScratchMemTop ( top );
      return SORT_BAD_DATA;
    }

    memset ( count, 0, countsize );
    itlnumbyt = item_length - key_size;
    InData = mdata;
    for ( j = 0; j < num_data; j++ ) {
      for ( i = 0; i < key_size; i++ ) {
        count[MAXCOUNT*i+(*InData)]++;
        InData++;
      }
      InData += itlnumbyt;
    }
    odd_swap = false;
#if __BYTE_ORDER == __LITTLE_ENDIAN
    for ( i = 0, InDataBase = mdata;
          i < key_size;
          i++, count += MAXCOUNT, InDataBase++ )
#else
#if __BYTE_ORDER == __BIG_ENDIAN
    for ( i = key_size-1, count += i*MAXCOUNT, InDataBase += i;
          i >= 0;
          i--, count -= MAXCOUNT, InDataBase-- )
#else
#error Either little endian or big endian byte ordering is assumed
#endif
#endif
    {
      if ( count[*InDataBase] < num_data ) {
        for ( j = 1, k = 0;  j < MAXCOUNT; k = j++ )
          count[j] += count[k];
        for ( j = num_data-1;  j >= 0;  j-- ) {
          InData = InDataBase + item_length*permut[j];
          k = *InData;
          count[k]--;
          auxpermut[count[k]] = permut[j];
        }
        iptr0 = auxpermut;  auxpermut = permut;  permut = iptr0;
        odd_swap = !odd_swap;
      }
    }
    if ( odd_swap )
      memcpy ( auxpermut, permut, num_data*sizeof(int) );
    pkv_SetScratchMemTop ( top );
  }
        /* convert the sequence with signed numbers back */
  switch ( key_type ) {
case ID_UNSIGNED:
    break;
case ID_SIGNED_INT:
    _ConvertInts ( key_size, item_length, num_data, mdata );
    break;
case ID_IEEE754_FLOAT:
    _ConvertFloats2 ( item_length, num_data, mdata );
    break;
case ID_IEEE754_DOUBLE:
    _ConvertDoubles2 ( item_length, num_data, mdata );
    break;
case ID_IEEE754_LONG_DOUBLE:
    _ConvertLDoubles2 ( item_length, num_data, mdata );
  }
  return SORT_OK;
} /*pkv_SortKernel*/

char pkv_SortPermute ( size_t item_length, unsigned int num_data, void *data,
                       unsigned int *permut )
{
  void *sp;
  int  i, k, m;
  char *tmp, *data1, *data2;

  sp = pkv_GetScratchMemTop ();
  tmp = (char*)pkv_GetScratchMem ( item_length );
  if ( !tmp )
    return SORT_NO_MEMORY;

  for ( i = 0; i < num_data-1; i++ ) {
    if ( i != permut[i] ) {
      k = i;
      data1 = ((char*)data) + item_length*k;
      memcpy ( tmp, data1, item_length );
      m = permut[k];
      while ( m != i ) {
        data2 = ((char*)data) + item_length*m;
        memcpy ( data1, data2, item_length );
        permut[k] = k;
        k = m;
        m = permut[k];
        data1 = data2;
      }
      memcpy ( data1, tmp, item_length );
      permut[k] = k;
    }
  }
  pkv_SetScratchMemTop ( sp );
  return SORT_OK;
} /*pkv_SortPermute*/

char pkv_SortFast ( size_t key_size, char key_type,
                    size_t item_length, size_t key_offset,
                    unsigned int num_data, void *data )
{
  void         *sp;
  unsigned int i;
  unsigned int *permut;
  char         result;

  sp = pkv_GetScratchMemTop ();
  permut = pkv_GetScratchMem ( num_data*sizeof(int) );
  if ( !permut )
    return SORT_NO_MEMORY;
  for ( i = 0; i < num_data; i++ )
    permut[i] = i;
  if ( (result = pkv_SortKernel ( key_size, key_type, item_length, key_offset,
                                  num_data, data, permut )) == SORT_OK )
    result = pkv_SortPermute ( item_length, num_data, data, permut );
  else
    result = SORT_NO_MEMORY;
  pkv_SetScratchMemTop ( sp );
  return result;
} /*pkv_SortFast*/

