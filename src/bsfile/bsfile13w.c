
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013, 2014                            */
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

boolean bsf_WriteCamera ( CameraRecd *camera, int ident )
{
  int sci;

  sci = bsf_current_indentation;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s {",
            bsf_keyword[BSF_SYMB_CAMERA-BSF_FIRST_KEYWORD] );
  BSFeol
  bsf_current_indentation += 2;
  BSFwci
  bsf_WriteIdent ( ident );
  bsf_current_length += fprintf ( bsf_output, "%s {",
            bsf_keyword[BSF_SYMB_FRAME-BSF_FIRST_KEYWORD] );
  bsf_current_length += fprintf ( bsf_output, "%d, %d, %d, %d}",
            camera->width, camera->height, camera->xmin, camera->ymin );
  BSFeol
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s ",
            bsf_keyword[BSF_SYMB_POSITION-BSF_FIRST_KEYWORD] );
  bsf_WritePointd ( 3, &camera->position.x );
  BSFeol
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "%s {",
            bsf_keyword[BSF_SYMB_EULERANGLES-BSF_FIRST_KEYWORD] );
  bsf_WriteDoubleNumber ( camera->psi );
  bsf_current_length += fprintf ( bsf_output, ", " );
  bsf_WriteDoubleNumber ( camera->theta );
  bsf_current_length += fprintf ( bsf_output, ", " );
  bsf_WriteDoubleNumber ( camera->phi );
  bsf_current_length += fprintf ( bsf_output, "}" );
  BSFeol
  if ( camera->zmin > 0.0 && camera->zmax > camera->zmin ) {
    BSFwci
    bsf_current_length += fprintf ( bsf_output, "%s { ",
              bsf_keyword[BSF_SYMB_DEPTH-BSF_FIRST_KEYWORD] );
    bsf_WriteDoubleNumber ( camera->zmin );
    bsf_current_length += fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->zmax );
    bsf_current_length += fprintf ( bsf_output, " }" );
    BSFeol
  }
  BSFwci
  if ( camera->parallel ) {  /* parallel projection camera */
    bsf_current_length += fprintf ( bsf_output, "%s { ",
              bsf_keyword[BSF_SYMB_PARALLEL-BSF_FIRST_KEYWORD] );
    bsf_WriteDoubleNumber ( camera->vd.para.wdt );
    bsf_current_length += fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->vd.para.hgh );
    bsf_current_length += fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->vd.para.diag );
    bsf_current_length += fprintf ( bsf_output, ", %d}",
                                    (int)camera->vd.para.dim_case );
  }
  else {  /* perspective projection camera */
    bsf_current_length += fprintf ( bsf_output, "%s { ",
              bsf_keyword[BSF_SYMB_PERSPECTIVE-BSF_FIRST_KEYWORD] );
    bsf_WriteDoubleNumber ( camera->vd.persp.f );
    bsf_current_length += fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->vd.persp.xi0 );
    bsf_current_length += fprintf ( bsf_output, ", " );
    bsf_WriteDoubleNumber ( camera->vd.persp.eta0 );
    bsf_current_length += fprintf ( bsf_output, "}" );
  }
  BSFeol
  bsf_current_indentation = sci;
  BSFwci
  bsf_current_length += fprintf ( bsf_output, "}" );
  BSFeol
  return true;
} /*bsf_WriteCamera*/

