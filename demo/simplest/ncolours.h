
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

void DrawColourSample ( int nc, short x, short y );
void RysujOkno ( xge_widget *er, boolean onscreen );
void WriteColourInfo ( short x, short y );
boolean OknoMsg ( xge_widget *er, int msg, int key, short x, short y );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );
void init_edwin ( void );
void destroy_edwin ( void );

