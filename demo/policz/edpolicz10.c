
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "splhole.h"

#define INFO_N_LINES 17
#define INFO_TXT_LENGTH 1025

static char *basis_info[INFO_N_LINES];
static char basis_info_text[INFO_TXT_LENGTH];

static char **FormatBasisInfo ( GHoleDomaind *domain, GHoptions *options,
                                char *title, char **info_ptr )
{
  int a, b, c, d;

  b = 6*hole_k+1;
  switch ( options->order ) {
case 1:
    if ( options->coons ) {
      a = g1h_V0SpaceDimd ( domain );
      sprintf ( info_ptr[0], "%s A = %d, B = %d; dim V0 = %d ",
                title, a, b, a );
    }
    else if ( options->bezier ) {
      a = g1h_V0SpaceDimd ( domain );
      c = g1h_ExtV0SpaceDimd ( domain ) - a;
      sprintf ( info_ptr[0], "%s A = %d, B = %d, C = %d; dim V0 = %d",
                title, a, b, c, a+c );
    }
    else if ( options->spline ) {
      if ( !options->spl_basis_valid )
        options->spl_basis_valid = g1h_ComputeSplBasisd ( domain, options->nk,
                              options->m1, options->m2 );
      if ( options->spl_basis_valid ) {
        g1h_DrawSplBasFuncNumd ( domain, &a, &b, &c, &d );
        sprintf ( info_ptr[0], "%s A = %d, B = %d, C = %d, D = %d; dim V0 = %d",
                  title, a, b, c, d, a+c+d );
      }
      else
        sprintf ( info_ptr[0], "%s %s", title, "Cannot create valid basis!" );
    }
    break;
case 2:
    if ( options->coons ) {
      a = g2h_V0SpaceDimd ( domain );
      sprintf ( info_ptr[0], "%s A = %d, B = %d; dim V0 = %d",
                title, a, b, a );
    }
    else if ( options->bezier ) {
      a = g2h_V0SpaceDimd ( domain );
      c = g2h_ExtV0SpaceDimd ( domain ) - a;
      sprintf ( info_ptr[0], "%s A = %d, B = %d, C = %d; dim V0 = %d",
                title, a, b, c, a+c );
    }
    else if ( options->spline ) {
      if ( !options->spl_basis_valid )
        options->spl_basis_valid = g2h_ComputeSplBasisd ( domain, options->nk,
                              options->m1, options->m2 );
      if ( options->spl_basis_valid ) {
        g2h_DrawSplBasFuncNumd ( domain, &a, &b, &c, &d );
        sprintf ( info_ptr[0], "%s A = %d, B = %d, C = %d, D = %d; dim V0 = %d",
                  title, a, b, c, d, a+c+d );
      }
      else
        sprintf ( info_ptr[0], "%s %s", title, "Cannot create valid basis!" );
    }
    break;
  }
  info_ptr[1] = info_ptr[0] + strlen(info_ptr[0])+1;
  return &info_ptr[1];
} /*FormatBasisInfo*/

void DisplayBasisInfo ( void )
{
  char **ptr;

  memset ( basis_info, 0, INFO_N_LINES*sizeof(char*) );
  memset ( basis_info_text, 0, INFO_TXT_LENGTH );
  basis_info[0] = &basis_info_text[0];
  ptr = &basis_info[0];
  ptr = FormatBasisInfo ( domain1, &options1, "1st basis:", ptr );
  *(ptr+1) = (*ptr)+1;
  ptr = FormatBasisInfo ( domain2, &options2, "2nd basis:", ptr+1 );
  *ptr = NULL;
  xge_DisplayInfoMessage ( basis_info, -1 );
} /*DisplayBasisInfo*/

