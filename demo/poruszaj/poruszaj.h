
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_PATH_LGT      1024
#define MAX_PATH_SHRT       63
#define MAX_FILENAME_LGT    64

#define TOP_MENU_HEIGHT     20
#define BOTTOM_MENU_HEIGHT  42
#define SIDE_MENU_WIDTH    123


#define cmdANIMATE 1

/* ///////////////////////////////////////////////////////////////////////// */
extern xge_3Dwind    wind;
extern xge_KnotWind  knwind;

extern xge_widget    *topmenu0, *sidemenu0, *geommenu0, *bottommenu0,
                     *popup0, *popup1, *popup2, *popup3;
extern xge_widget    *side0a, *side0b, *side0c;
extern xge_widget    *bottom0a, *bottom0b;

extern xge_listbox   dirlist1, filelist1, dirlist2, filelist2;
extern xge_string_ed filename_editor;
extern const char    file_filter[];
extern const char    file_ext[];
extern char          initial_directory[MAX_PATH_LGT+1];
extern char          current_directory[MAX_PATH_LGT+1], current_dir[MAX_PATH_SHRT+1];
extern char          filename[MAX_FILENAME_LGT+1];

extern xge_scroll_widget scrollw0;

extern boolean sw_palmlod1, sw_palmlod2, sw_palmlod3;

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *TopMenu0Init ( xge_widget *prev );
int TopMenu0CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *SideMenu0Init ( xge_widget *prev );
int SideMenu0aCallBack ( xge_widget *er, int msg, int key, short x, short y );
int SideMenu0bCallBack ( xge_widget *er, int msg, int key, short x, short y );
int SideMenu0cCallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *GeomMenu0Init ( xge_widget *prev );
void RysujOkno ( xge_widget *er, boolean onscreen );
int GeomMenu0CallBack ( xge_widget *er, int msg, int key, short x, short y );

void MyDrawKnotWind ( xge_widget *er, boolean onscreen );
void MyDrawButton ( xge_widget *er, boolean onscreen );
xge_widget *BottomMenu0Init ( xge_widget *prev );
void AnimaOn ( void );
double AnimaCurrentTime ( clock_t *ct );
boolean Animate ( boolean nudge );
int BottomMenu0CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *Popup0Init ( void );
void SetCurrentDir ( void );
boolean FilenameCorrect ( char *fn );
void OpenSaveAsPopup ( void );
int Popup0CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *Popup1Init ( void );
void Popup1ChangeDir ( void );
boolean Popup1OpenFile ( void );
int Popup1CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *Popup2Init ( void );
void Popup2ChangeDir ( void );
boolean Popup2SaveFile ( void );
int Popup2CallBack ( xge_widget *er, int msg, int key, short x, short y );

xge_widget *Popup3Init ( void );
int Popup3CallBack ( xge_widget *er, int msg, int key, short x, short y );

int NonWidgetCallBack ( int msg, int key, short x, short y );

int CallBack ( xge_widget *er, int msg, int key, short x, short y );

void init_edwin ( void );
void destroy_edwin ( void );

int main ( int argc, char *argv[] );

