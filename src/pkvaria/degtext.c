
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
#include <ctype.h>
#include <math.h>

#include "pkvaria.h"

void pkv_RadToDegreeStr ( double angle, char *txt )
{
  int deg, min, sec;

  if ( angle < 0 )
    *txt++ = '-';
  angle = fabs(angle)*180.0/PI;
  deg = (int)angle;
  angle = 60.0*(angle-deg);
  min = (int)angle;
  sec = 60.0*(angle-min);
  sprintf ( txt, "%d.%02d'%02d\"", deg, min, sec );
} /*pkv_RadToDegreeStr*/

boolean pkv_DegreeStrToRad ( char *txt, double *angle )
{
  int     deg, min, sec;
  boolean neg;

  while ( *txt == ' ' ) txt++;
  if ( !(*txt) )
    return false;
  neg = false;
  if ( *txt == '+' ) txt++;
  else if ( *txt == '-' ) { txt++;  neg = true; }
  if ( !isdigit( *txt ) )
    return false;
  deg = min = sec = 0;
  deg = (int)strtol ( txt, &txt, 10 );
  if ( *txt == '.' ) {
    txt++;
    if ( isdigit ( *txt ) ) {
      min = (int) strtol ( txt, &txt, 10 );
      if ( *txt == '\'' ) {
        txt++;
        if ( isdigit ( *txt ) ) {
          sec = (int) strtol ( txt, &txt, 10 );
          if ( *txt != '\"' )
            return false;
        }
      }
      else if ( *txt == '\"' ) {
        sec = min;
        min = 0;
      }
    }
  }
  if ( min > 59 || sec > 59 )
    return false;
  *angle = ((double)deg + (double)min/60.0 + (double)sec/3600.0)*PI/180.0;
  if ( neg )
    *angle = -(*angle);
  return true;
} /*pkv_DegreeStrToRad*/

