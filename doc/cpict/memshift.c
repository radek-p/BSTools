
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h> 
#include <string.h>

#include "pkvaria.h"
#include "psout.h"

char fn2[12] = "memshift.ps";

int main ()
{
  int i;

  ps_WriteBBox ( 12, 0, 311, 100 );
  ps_OpenFile ( fn2, 600 );

  ps_Set_Gray ( 0.68 );
  for ( i = 0; i < 4; i++ )
    ps_Fill_Rect ( 280, 80, (float)(200+500*i), 500 );
  ps_Set_Line_Width ( 3 );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < 4; i++ )
    ps_Draw_Rect ( 280, 80, (float)(200+500*i), 500 );
  ps_Draw_Line ( 200, 500, 2600, 500 );
  ps_Draw_Line ( 200, 580, 2600, 580 );

  ps_Set_Gray ( 0.68 );
  for ( i = 0; i < 4; i++ )
    ps_Fill_Rect ( 280, 80, (float)(300+500*i), 200 );
  ps_Set_Line_Width ( 3 );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < 4; i++ )
    ps_Draw_Rect ( 280, 80, (float)(300+500*i), 200 );
  ps_Draw_Line ( 200, 198.5, 200, 281.5 );
  ps_Draw_Line ( 200, 200, 2600, 200 );
  ps_Draw_Line ( 200, 280, 2600, 280 );

  for ( i = 0; i < 4; i++ ) {
    psl_SetLine ( (float)(200+500*i), 500, (float)(300+500*i), 280, 0.0, 1.0 );
    psl_Draw ( 0.1, 0.85, 1.5 );
    psl_Arrow ( 0.9, true );
    psl_SetLine ( (float)(480+500*i), 500, (float)(580+500*i), 280, 0.0, 1.0 );
    psl_Draw ( 0.1, 0.85, 1.5 );
    psl_Arrow ( 0.9, true );
  }
  ps_Draw_Line ( 200, 590, 200, 790 );
  ps_Draw_Line ( 480, 590, 480, 790 );
  psl_SetLine ( 200, 750, 480, 750, 0.0, 1.0 );
  psl_Draw ( 0.05, 0.95, 1.0 );
  psl_Arrow ( 0.0, false );
  psl_Arrow ( 1.0, true );

  ps_Draw_Line ( 700, 590, 700, 790 );
  ps_Draw_Line ( 1200, 590, 1200, 790 );
  psl_SetLine ( 700, 750, 1200, 750, 0.0, 1.0 );
  psl_Draw ( 0.05, 0.95, 1.0 );
  psl_Arrow ( 0.0, false );
  psl_Arrow ( 1.0, true );

  ps_Draw_Line ( 800, 190, 800, -10 );
  ps_Draw_Line ( 1300, 190, 1300, -10 );
  psl_SetLine ( 800, 30, 1300, 30, 0.0, 1.0 );
  psl_Draw ( 0.05, 0.95, 1.0 );
  psl_Arrow ( 0.0, false );
  psl_Arrow ( 1.0, true );

  psl_SetLine ( 200, 50, 200, 200, 0.0, 1.0 );
  psl_Draw ( 0.2, 0.95, 1.0 );
  psl_Arrow ( 1.0, true );

  ps_Write_Command ( "/Courier findfont 83 scalefont setfont" );
  ps_Write_Command ( "/Center { /s exch def s stringwidth pop -2 div 0 rmoveto s show } def" );
  ps_Write_Command ( "200 20 moveto (data) Center" );
  ps_Write_Command ( "1050 50 moveto (pitch) Center" );
  ps_Write_Command ( "950 770 moveto (pitch) Center" );
  ps_Write_Command ( "340 770 moveto (rowlen) Center" );

  ps_CloseFile ();
  printf ( "%s\n", fn2 );

  exit ( 0 );
} /*main*/

