
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h> 
#include <math.h> 
#include <malloc.h>
#include <string.h> 
#include <ctype.h>  

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "bsfile.h"

#include "bsfprivate.h"

boolean bsf_WriteTrimmedDomaind ( int nelem, const mbs_polycurved *bound )
{
  int     i;
  int     sci;
  boolean one_contour, new_contour;
  boolean result;

  if ( nelem < 1 )
    return false;

  result = true;
        /* find out, whether there is just one contour, or more */
        /* and make some simple data correctness tests */
  one_contour = true;
  for ( i = 0; i < nelem; i++ ) {
    if ( i < nelem-1 && bound[i].closing )
      one_contour = false;
          /* each element must be planar, may be rational */
    if ( bound[i].cpdimen != 2 && bound[i].cpdimen != 3 )
      return false;
  }
  if ( !bound[nelem-1].closing )
    return false;

        /* now the actual writing */
  sci = bsf_current_indentation;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s %s {",
            bsf_keyword[BSF_SYMB_TRIMMED-BSF_FIRST_KEYWORD],
            bsf_keyword[BSF_SYMB_DOMAIN-BSF_FIRST_KEYWORD] );
  BSFeol
  bsf_current_indentation += 2;
  new_contour = !one_contour;
  for ( i = 0; i < nelem; i++ ) {
        /* if more than one contour than write additional braces around */
        /* each contour */
    if ( new_contour ) {
      BSFwci
      bsf_current_length += fprintf ( bsf_output, "{" );
      BSFeol
      new_contour = false;
      bsf_current_indentation += 2;
    }
        /* write the curves and polylines; they have no attributes */
        /* in the current version of the .bs files language */
    if ( bound[i].knots ) {              /* write a B-spline curve */
      result &= bsf_WriteBSplineCurved ( 2, bound[i].cpdimen,
                                         bound[i].cpdimen > 2, bound[i].degree,
                                         bound[i].lastknot, bound[i].knots,
                                         false, bound[i].points,
                                         NULL, bound[i].ident, NULL, NULL );
    }
    else if ( bound[i].lastknot <= 0 ) { /* write a Bezier curve */
      result &= bsf_WriteBezierCurved ( 2, bound[i].cpdimen,
                                        bound[i].cpdimen > 2, bound[i].degree,
                                        bound[i].points, NULL, bound[i].ident,
                                        NULL, NULL );
    }
    else {                               /* write a polyline */
      result &= bsf_WritePolylined ( 2, bound[i].cpdimen,
                                     bound[i].cpdimen > 2, false,
                                     bound[i].lastknot+1, bound[i].points,
                                     NULL, bound[i].ident, NULL, NULL );
    }
        /* if more than one contour then close the additional brace */
    if ( !one_contour && bound[i].closing ) {
      bsf_current_indentation -= 2;
      BSFwci
      bsf_current_length += fprintf ( bsf_output, "%s", "} %end contour" );
      BSFeol
      new_contour = true;
    }
  }
  bsf_current_indentation = sci;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "}" );
  BSFeol
  return result;
} /*bsf_WriteTrimmedDomaind*/

