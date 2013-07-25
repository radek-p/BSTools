
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

void RysujOkno ( xge_widget *er, boolean onscreen );
boolean OknoMsg ( xge_widget *er, int msg, int key, short x, short y );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );
void init_edwin ( void );
void destroy_edwin ( void );

void ProcessChildMessage ( int msg, int size );
