
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pkvaria.h"
#include "psout.h"

char fn1[12] = "memory.ps";

int main ()
{
  int i;

  ps_WriteBBox ( 12, 0, 311, 100 );
  ps_OpenFile ( fn1, 600 );

  ps_Set_Gray ( 0.68 );
  for ( i = 0; i < 4; i++ )
    ps_Fill_Rect ( 280.0, 80.0, (float)(200.0+500.0*i), 500.0 );
  ps_Set_Line_Width ( 3 );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < 4; i++ )
    ps_Draw_Rect ( 280.0, 80.0, (float)(200.0+500.0*i), 500.0 );
  ps_Draw_Line ( 200.0, 500.0, 2600.0, 500.0 );
  ps_Draw_Line ( 200.0, 580.0, 2600.0, 580.0 );

  ps_Set_Gray ( 0.68 );
  for ( i = 0; i < 4; i++ )
    ps_Fill_Rect ( 280.0, 80.0, (float)(200.0+600.0*i), 200.0 );
  ps_Set_Line_Width ( 3 );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < 4; i++ )
    ps_Draw_Rect ( 280.0, 80.0, (float)(200.0+600.0*i), 200.0 );
  ps_Draw_Line ( 200.0, 200.0, 2600.0, 200.0 );
  ps_Draw_Line ( 200.0, 280.0, 2600.0, 280.0 );

  for ( i = 0; i < 4; i++ ) {
    psl_SetLine ( (float)(200.0+500.0*i), 500.0, (float)(200.0+600.0*i), 280.0,
                  0.0, 1.0 );
    psl_Draw ( 0.1, 0.85, 1.5 );
    psl_Arrow ( 0.9, true );
    psl_SetLine ( (float)(480.0+500.0*i), 500.0, (float)(480.0+600.0*i), 280.0,
                  0.0, 1.0 );
    psl_Draw ( 0.1, 0.85, 1.5 );
    psl_Arrow ( 0.9, true );
  }
  ps_Draw_Line ( 200.0, 590.0, 200.0, 790.0 );
  ps_Draw_Line ( 480.0, 590.0, 480.0, 790.0 );
  psl_SetLine ( 200.0, 750.0, 480.0, 750.0, 0.0, 1.0 );
  psl_Draw ( 0.05, 0.95, 1.0 );
  psl_Arrow ( 0.0, false );
  psl_Arrow ( 1.0, true );

  ps_Draw_Line ( 700.0, 590.0, 700.0, 790.0 );
  ps_Draw_Line ( 1200.0, 590.0, 1200.0, 790.0 );
  psl_SetLine ( 700.0, 750.0, 1200.0, 750.0, 0.0, 1.0 );
  psl_Draw ( 0.05, 0.95, 1.0 );
  psl_Arrow ( 0.0, false );
  psl_Arrow ( 1.0, true );

  ps_Draw_Line ( 800.0, 190.0, 800.0, -10 );
  ps_Draw_Line ( 1400.0, 190.0, 1400.0, -10 );
  psl_SetLine ( 800.0, 30.0, 1400.0, 30.0, 0.0, 1.0 );
  psl_Draw ( 0.05, 0.95, 1.0 );
  psl_Arrow ( 0.0, false );
  psl_Arrow ( 1.0, true );

  psl_SetLine ( 200.0, 50.0, 200.0, 200.0, 0.0, 1.0 );
  psl_Draw ( 0.2, 0.95, 1.0 );
  psl_Arrow ( 1.0, true );

  ps_Write_Command ( "/Courier findfont 83 scalefont setfont" );
  ps_Write_Command ( "/Center { /s exch def s stringwidth pop -2 div 0 rmoveto s show } def" );
  ps_Write_Command ( "200 20 moveto (data) Center" );
  ps_Write_Command ( "1100 50 moveto (outpitch) Center" );
  ps_Write_Command ( "950 770 moveto (inpitch) Center" );
  ps_Write_Command ( "340 770 moveto (rowlen) Center" );

  printf ( "%s\n", fn1 );
  ps_CloseFile ();

  exit ( 0 );
} /*main*/

