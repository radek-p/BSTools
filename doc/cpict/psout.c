
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include "psout.h"

int main ()
{
  ps_WriteBBox ( 12, 13, 264, 80 );
  ps_OpenFile ( "psout.ps", 600 );
  psl_SetLine ( 200, 600, 2200, 100, 0.0, 4.0 );
  psl_Draw ( 0.5, 1.5, 2.0 );
  psl_Draw ( 2.0, 3.0, 6.0 );
  psl_ADraw ( 3.5, 4.0, 0.0, -arrowl, 1.0 );
  psl_DrawEye ( 0.0, 1, 1.2, 0.15 );
  psl_HTick ( 1.0, false );
  psl_HTick ( 1.5, true );
  psl_Tick ( 2.0 );
  psl_BlackHighLTrMark ( 2.5 );
  psl_LTrMark ( 3.0 );
  psl_Dot ( 3.5 );
  psl_Arrow ( 4.0, true );
  ps_CloseFile ();
  exit ( 0 );
} /*main*/

