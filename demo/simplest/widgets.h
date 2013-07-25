
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */
boolean ColourWMsg ( xge_widget *er, int msg, int key, short x, short y );
void ColourWRedraw ( xge_widget *er, boolean onscreen );
void RysujOkno ( xge_widget *er, boolean onscreen );
boolean OknoMsg ( xge_widget *er, int msg, int key, short x, short y );
boolean MyMenuMsg ( xge_widget *er, int msg, int key, short x, short y );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );
void PrepareCursors ( void );
void init_edwin ( void );
void destroy_edwin ( void );

