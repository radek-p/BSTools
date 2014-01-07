
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
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

boolean bsf_WriteCamera ( CameraRecd *camera )
{
  fprintf ( bsf_output, "%s {\n",
            bsf_keyword[BSF_SYMB_CAMERA-BSF_FIRST_KEYWORD] );
  fprintf ( bsf_output, "  %s { ",
            bsf_keyword[BSF_SYMB_FRAME-BSF_FIRST_KEYWORD] );
  fprintf ( bsf_output, "%d, %d, %d, %d }",
            camera->width, camera->height, camera->xmin, camera->ymin );
  fprintf ( bsf_output, "\n  %s ",
            bsf_keyword[BSF_SYMB_POSITION-BSF_FIRST_KEYWORD] );
  bsf_WritePointd ( 3, &camera->position.x );
  fprintf ( bsf_output, "\n  %s { ",
            bsf_keyword[BSF_SYMB_EULERANGLES-BSF_FIRST_KEYWORD] );
  bsf_WriteDoubleNumber ( camera->psi );
  fprintf ( bsf_output, ", " );
  bsf_WriteDoubleNumber ( camera->theta );
  fprintf ( bsf_output, ", " );
  bsf_WriteDoubleNumber ( camera->phi );
  fprintf ( bsf_output, " }\n" );
  if ( camera->zmin > 0.0 && camera->zmax > camera->zmin ) {
    fprintf ( bsf_output, "  %s { ",
              bsf_keyword[BSF_SYMB_DEPTH-BSF_FIRST_KEYWORD] );
    bsf_WriteDoubleNumber ( camera->zmin );
    fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->zmax );
    fprintf ( bsf_output, " }\n" );
  }
  if ( camera->parallel ) {  /* parallel projection camera */
    fprintf ( bsf_output, "  %s { ",
              bsf_keyword[BSF_SYMB_PARALLEL-BSF_FIRST_KEYWORD] );
    bsf_WriteDoubleNumber ( camera->vd.para.wdt );
    fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->vd.para.hgh );
    fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->vd.para.diag );
    fprintf ( bsf_output, ", %d", (int)camera->vd.para.dim_case );
  }
  else {  /* perspective projection camera */
    fprintf ( bsf_output, "  %s { ",
              bsf_keyword[BSF_SYMB_PERSPECTIVE-BSF_FIRST_KEYWORD] );
    bsf_WriteDoubleNumber ( camera->vd.persp.f );
    fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->vd.persp.xi0 );
    fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->vd.persp.eta0 );
  }
  fprintf ( bsf_output, " }\n}\n" );
  return true;
} /*bsf_WriteCamera*/

