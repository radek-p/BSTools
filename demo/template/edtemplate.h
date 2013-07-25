
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#define MAX_PATH_LGT       1024
#define MAX_FILENAME_LGT   64

#define STATUS_LINE_LENGTH ((xge_MAX_WIDTH-110)/6)

extern int win0, win1;
extern xge_widget *menu0, *menu1, *menu2;
extern xge_widget *menu01alist;
extern xge_widget *popup00, *popup01, *popup02;
extern xge_widget *status0, *status0sw;
extern char statustext0[];
extern boolean status_0_sw;

extern xge_listbox   filelist, dirlist;
extern xge_string_ed filename_editor;
extern char filename[];
extern char current_directory[];
extern char file_filter[];

extern xge_widget *menu3, *menu4, *menu5;
extern xge_widget *menu14alist;
extern xge_widget *knwind;
extern xge_widget *status1, *status1sw;
extern char statustext1[];
extern boolean status_1_sw;


void RysujOkno ( xge_widget *er, boolean onscreen );
void RysujPOkno ( xge_widget *er, boolean onscreen );

void RysujDomOkno ( xge_widget *er, boolean onscreen );
void RysujKnotOkno ( xge_widget *er, boolean onscreen );
boolean KnotOknoMsg ( xge_widget *er, int msg, int key, short x, short y );

void init_edwin ( void );
void destroy_edwin ( void );

int ProcessKey ( int key );
int Win0CallBack ( xge_widget *er, int msg, int key, short x, short y );
int Win1CallBack ( xge_widget *er, int msg, int key, short x, short y );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );

boolean FilenameCorrect ( const char *fn );
void OpenFile ( const char *fn );
void SaveFile ( const char *fn );

void ResizeWinStatus ( int win );
void SetStatusLineText ( int win, const char *text );
void StatusLineOnOff ( int win );
void Notify2DTrans ( int win, xge_2Dwind *_2Dwin );
void Notify2DTransChange ( int win, xge_2Dwind *_2Dwin );
void Notify3DTrans ( int win, xge_3Dwind *_3Dwin );
void Notify3DTransChange ( int win, xge_3Dwind *_3Dwin );

