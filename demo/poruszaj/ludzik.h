
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

/* parametry artykulacji */
#define N_PARAMS 21

/* tu musza byc kolejne liczby calkowite od 0, bo to indeksy do tablicy */
#define PAR_GLOWA_1  0
#define PAR_GLOWA_2  1
#define PAR_GLOWA_3  2

#define PAR_LBARK_1  3
#define PAR_LBARK_2  4
#define PAR_LBARK_3  5
#define PAR_LLOK_1   6
#define PAR_LLOK_2   7

#define PAR_PBARK_1  8
#define PAR_PBARK_2  9
#define PAR_PBARK_3 10
#define PAR_PLOK_1  11
#define PAR_PLOK_2  12

#define PAR_LBIO_1  13
#define PAR_LBIO_2  14
#define PAR_LBIO_3  15
#define PAR_LKOL_1  16

#define PAR_PBIO_1  17
#define PAR_PBIO_2  18
#define PAR_PBIO_3  19
#define PAR_PKOL_1  20

extern double art_param[N_PARAMS];
extern boolean sw_spline_hands;
extern boolean sw_draw_bsnets;

void ResetParameters ( void );
void InitCharacter ( void );
void DisplayCharacter ( void );

