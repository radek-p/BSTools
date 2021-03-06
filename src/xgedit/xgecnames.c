
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "xgedit.h"

const char *xge_colour_name[] =
 {"AliceBlue",            /*   0 */
  "AntiqueWhite",         /*   1 */
  "AntiqueWhite1",        /*   2 */
  "AntiqueWhite2",        /*   3 */
  "AntiqueWhite3",        /*   4 */
  "AntiqueWhite4",        /*   5 */
  "Aquamarine",           /*   6 */
  "Aquamarine1",          /*   7 */
  "Aquamarine2",          /*   8 */
  "Aquamarine3",          /*   9 */
  "Aquamarine4",          /*  10 */
  "Azure",                /*  11 */
  "Azure1",               /*  12 */
  "Azure2",               /*  13 */
  "Azure3",               /*  14 */
  "Azure4",               /*  15 */
  "Beige",                /*  16 */
  "Bisque",               /*  17 */
  "Bisque1",              /*  18 */
  "Bisque2",              /*  19 */
  "Bisque3",              /*  20 */
  "Bisque4",              /*  21 */
  "Black",                /*  22 */
  "BlanchedAlmond",       /*  23 */
  "Blue",                 /*  24 */
  "Blue1",                /*  25 */
  "Blue2",                /*  26 */
  "Blue3",                /*  27 */
  "Blue4",                /*  28 */
  "Blue5",                /*  29 */
  "Blue6",                /*  30 */
  "Blue7",                /*  31 */
  "BlueViolet",           /*  32 */
  "Brown",                /*  33 */
  "Brown1",               /*  34 */
  "Brown2",               /*  35 */
  "Brown3",               /*  36 */
  "Brown4",               /*  37 */
  "Burlywood",            /*  38 */
  "Burlywood1",           /*  39 */
  "Burlywood2",           /*  40 */
  "Burlywood3",           /*  41 */
  "Burlywood4",           /*  42 */
  "CadetBlue",            /*  43 */
  "CadetBlue1",           /*  44 */
  "CadetBlue2",           /*  45 */
  "CadetBlue3",           /*  46 */
  "CadetBlue4",           /*  47 */
  "Chartreuse",           /*  48 */
  "Chartreuse1",          /*  49 */
  "Chartreuse2",          /*  50 */
  "Chartreuse3",          /*  51 */
  "Chartreuse4",          /*  52 */
  "Chocolate",            /*  53 */
  "Chocolate1",           /*  54 */
  "Chocolate2",           /*  55 */
  "Chocolate3",           /*  56 */
  "Chocolate4",           /*  57 */
  "Coral",                /*  58 */
  "Coral1",               /*  59 */
  "Coral2",               /*  60 */
  "Coral3",               /*  61 */
  "Coral4",               /*  62 */
  "CornflowerBlue",       /*  63 */
  "Cornsilk",             /*  64 */
  "Cornsilk1",            /*  65 */
  "Cornsilk2",            /*  66 */
  "Cornsilk3",            /*  67 */
  "Cornsilk4",            /*  68 */
  "Cyan",                 /*  69 */
  "Cyan1",                /*  70 */
  "Cyan2",                /*  71 */
  "Cyan3",                /*  72 */
  "Cyan4",                /*  73 */
  "DarkBlue",             /*  74 */
  "DarkCyan",             /*  75 */
  "DarkGoldenrod",        /*  76 */
  "DarkGoldenrod1",       /*  77 */
  "DarkGoldenrod2",       /*  78 */
  "DarkGoldenrod3",       /*  79 */
  "DarkGoldenrod4",       /*  80 */
  "DarkGreen",            /*  81 */
  "DarkGrey",             /*  82 */
  "DarkKhaki",            /*  83 */
  "DarkMagenta",          /*  84 */
  "DarkOliveGreen",       /*  85 */
  "DarkOliveGreen1",      /*  86 */
  "DarkOliveGreen2",      /*  87 */
  "DarkOliveGreen3",      /*  88 */
  "DarkOliveGreen4",      /*  89 */
  "DarkOrange",           /*  90 */
  "DarkOrange1",          /*  91 */
  "DarkOrange2",          /*  92 */
  "DarkOrange3",          /*  93 */
  "DarkOrange4",          /*  94 */
  "DarkOrchid",           /*  95 */
  "DarkOrchid1",          /*  96 */
  "DarkOrchid2",          /*  97 */
  "DarkOrchid3",          /*  98 */
  "DarkOrchid4",          /*  99 */
  "DarkRed",              /* 100 */
  "DarkSalmon",           /* 101 */
  "DarkSeaGreen",         /* 102 */
  "DarkSeaGreen1",        /* 103 */
  "DarkSeaGreen2",        /* 104 */
  "DarkSeaGreen3",        /* 105 */
  "DarkSeaGreen4",        /* 106 */
  "DarkSlateBlue",        /* 107 */
  "DarkSlateGrey",        /* 108 */
  "DarkSlateGrey1",       /* 109 */
  "DarkSlateGrey2",       /* 110 */
  "DarkSlateGrey3",       /* 111 */
  "DarkSlateGrey4",       /* 112 */
  "DarkTurquoise",        /* 113 */
  "DarkViolet",           /* 114 */
  "DeepPink",             /* 115 */
  "DeepPink1",            /* 116 */
  "DeepPink2",            /* 117 */
  "DeepPink3",            /* 118 */
  "DeepPink4",            /* 119 */
  "DeepSkyBlue",          /* 120 */
  "DeepSkyBlue1",         /* 121 */
  "DeepSkyBlue2",         /* 122 */
  "DeepSkyBlue3",         /* 123 */
  "DeepSkyBlue4",         /* 124 */
  "DimGrey",              /* 125 */
  "DodgerBlue",           /* 126 */
  "DodgerBlue1",          /* 127 */
  "DodgerBlue2",          /* 128 */
  "DodgerBlue3",          /* 129 */
  "DodgerBlue4",          /* 130 */
  "Firebrick",            /* 131 */
  "Firebrick1",           /* 132 */
  "Firebrick2",           /* 133 */
  "Firebrick3",           /* 134 */
  "Firebrick4",           /* 135 */
  "FloralWhite",          /* 136 */
  "ForestGreen",          /* 137 */
  "Gainsboro",            /* 138 */
  "GhostWhite",           /* 139 */
  "Gold",                 /* 140 */
  "Gold1",                /* 141 */
  "Gold2",                /* 142 */
  "Gold3",                /* 143 */
  "Gold4",                /* 144 */
  "Goldenrod",            /* 145 */
  "Goldenrod1",           /* 146 */
  "Goldenrod2",           /* 147 */
  "Goldenrod3",           /* 148 */
  "Goldenrod4",           /* 149 */
  "Green",                /* 150 */
  "Green1",               /* 151 */
  "Green2",               /* 152 */
  "Green3",               /* 153 */
  "Green4",               /* 154 */
  "Green5",               /* 155 */
  "Green6",               /* 156 */
  "Green7",               /* 157 */
  "GreenYellow",          /* 158 */
  "Grey",                 /* 159 */
  "Grey1",                /* 160 */
  "Grey2",                /* 161 */
  "Grey3",                /* 162 */
  "Grey4",                /* 163 */
  "Grey5",                /* 164 */
  "Grey6",                /* 165 */
  "Grey7",                /* 166 */
  "Honeydew",             /* 167 */
  "Honeydew1",            /* 168 */
  "Honeydew2",            /* 169 */
  "Honeydew3",            /* 170 */
  "Honeydew4",            /* 171 */
  "HotPink",              /* 172 */
  "HotPink1",             /* 173 */
  "HotPink2",             /* 174 */
  "HotPink3",             /* 175 */
  "HotPink4",             /* 176 */
  "IndianRed",            /* 177 */
  "IndianRed1",           /* 178 */
  "IndianRed2",           /* 179 */
  "IndianRed3",           /* 180 */
  "IndianRed4",           /* 181 */
  "Ivory",                /* 182 */
  "Ivory1",               /* 183 */
  "Ivory2",               /* 184 */
  "Ivory3",               /* 185 */
  "Ivory4",               /* 186 */
  "Khaki",                /* 187 */
  "Khaki1",               /* 188 */
  "Khaki2",               /* 189 */
  "Khaki3",               /* 190 */
  "Khaki4",               /* 191 */
  "Lavender",             /* 192 */
  "LavenderBlush",        /* 193 */
  "LavenderBlush1",       /* 194 */
  "LavenderBlush2",       /* 195 */
  "LavenderBlush3",       /* 196 */
  "LavenderBlush4",       /* 197 */
  "LawnGreen",            /* 198 */
  "LemonChiffon",         /* 199 */
  "LemonChiffon1",        /* 200 */
  "LemonChiffon2",        /* 201 */
  "LemonChiffon3",        /* 202 */
  "LemonChiffon4",        /* 203 */
  "LightBlue",            /* 204 */
  "LightBlue1",           /* 205 */
  "LightBlue2",           /* 206 */
  "LightBlue3",           /* 207 */
  "LightBlue4",           /* 208 */
  "LightCoral",           /* 209 */
  "LightCyan",            /* 210 */
  "LightCyan1",           /* 211 */
  "LightCyan2",           /* 212 */
  "LightCyan3",           /* 213 */
  "LightCyan4",           /* 214 */
  "LightGoldenrod",       /* 215 */
  "LightGoldenrod1",      /* 216 */
  "LightGoldenrod2",      /* 217 */
  "LightGoldenrod3",      /* 218 */
  "LightGoldenrod4",      /* 219 */
  "LightGoldenrodYellow", /* 220 */
  "LightGreen",           /* 221 */
  "LightGrey",            /* 222 */
  "LightPink",            /* 223 */
  "LightPink1",           /* 224 */
  "LightPink2",           /* 225 */
  "LightPink3",           /* 226 */
  "LightPink4",           /* 227 */
  "LightSalmon",          /* 228 */
  "LightSalmon1",         /* 229 */
  "LightSalmon2",         /* 230 */
  "LightSalmon3",         /* 231 */
  "LightSalmon4",         /* 232 */
  "LightSeaGreen",        /* 233 */
  "LightSkyBlue",         /* 234 */
  "LightSkyBlue1",        /* 235 */
  "LightSkyBlue2",        /* 236 */
  "LightSkyBlue3",        /* 237 */
  "LightSkyBlue4",        /* 238 */
  "LightSlateBlue",       /* 239 */
  "LightSlateGrey",       /* 240 */
  "LightSteelBlue",       /* 241 */
  "LightSteelBlue1",      /* 242 */
  "LightSteelBlue2",      /* 243 */
  "LightSteelBlue3",      /* 244 */
  "LightSteelBlue4",      /* 245 */
  "LightYellow",          /* 246 */
  "LightYellow1",         /* 247 */
  "LightYellow2",         /* 248 */
  "LightYellow3",         /* 249 */
  "LightYellow4",         /* 250 */
  "LimeGreen",            /* 251 */
  "Linen",                /* 252 */
  "Magenta",              /* 253 */
  "Magenta1",             /* 254 */
  "Magenta2",             /* 255 */
  "Magenta3",             /* 256 */
  "Magenta4",             /* 257 */
  "Maroon",               /* 258 */
  "Maroon1",              /* 259 */
  "Maroon2",              /* 260 */
  "Maroon3",              /* 261 */
  "Maroon4",              /* 262 */
  "MediumAquamarine",     /* 263 */
  "MediumBlue",           /* 264 */
  "MediumOrchid",         /* 265 */
  "MediumOrchid1",        /* 266 */
  "MediumOrchid2",        /* 267 */
  "MediumOrchid3",        /* 268 */
  "MediumOrchid4",        /* 269 */
  "MediumPurple",         /* 270 */
  "MediumPurple1",        /* 271 */
  "MediumPurple2",        /* 272 */
  "MediumPurple3",        /* 273 */
  "MediumPurple4",        /* 274 */
  "MediumSeaGreen",       /* 275 */
  "MediumSlateBlue",      /* 276 */
  "MediumSpringGreen",    /* 277 */
  "MediumTurquoise",      /* 278 */
  "MediumVioletRed",      /* 279 */
  "MidnightBlue",         /* 280 */
  "MintCream",            /* 281 */
  "MistyRose",            /* 282 */
  "MistyRose1",           /* 283 */
  "MistyRose2",           /* 284 */
  "MistyRose3",           /* 285 */
  "MistyRose4",           /* 286 */
  "Moccasin",             /* 287 */
  "NavajoWhite",          /* 288 */
  "NavajoWhite1",         /* 289 */
  "NavajoWhite2",         /* 290 */
  "NavajoWhite3",         /* 291 */
  "NavajoWhite4",         /* 292 */
  "NavyBlue",             /* 293 */
  "OldLace",              /* 294 */
  "OliveDrab",            /* 295 */
  "OliveDrab1",           /* 296 */
  "OliveDrab2",           /* 297 */
  "OliveDrab3",           /* 298 */
  "OliveDrab4",           /* 299 */
  "Orange",               /* 300 */
  "Orange1",              /* 301 */
  "Orange2",              /* 302 */
  "Orange3",              /* 303 */
  "Orange4",              /* 304 */
  "OrangeRed",            /* 305 */
  "OrangeRed1",           /* 306 */
  "OrangeRed2",           /* 307 */
  "OrangeRed3",           /* 308 */
  "OrangeRed4",           /* 309 */
  "Orchid",               /* 310 */
  "Orchid1",              /* 311 */
  "Orchid2",              /* 312 */
  "Orchid3",              /* 313 */
  "Orchid4",              /* 314 */
  "PaleGoldenrod",        /* 315 */
  "PaleGreen",            /* 316 */
  "PaleGreen1",           /* 317 */
  "PaleGreen2",           /* 318 */
  "PaleGreen3",           /* 319 */
  "PaleGreen4",           /* 320 */
  "PaleTurquoise",        /* 321 */
  "PaleTurquoise1",       /* 322 */
  "PaleTurquoise2",       /* 323 */
  "PaleTurquoise3",       /* 324 */
  "PaleTurquoise4",       /* 325 */
  "PaleVioletRed",        /* 326 */
  "PaleVioletRed1",       /* 327 */
  "PaleVioletRed2",       /* 328 */
  "PaleVioletRed3",       /* 329 */
  "PaleVioletRed4",       /* 330 */
  "PapayaWhip",           /* 331 */
  "PeachPuff",            /* 332 */
  "PeachPuff1",           /* 333 */
  "PeachPuff2",           /* 334 */
  "PeachPuff3",           /* 335 */
  "PeachPuff4",           /* 336 */
  "Peru",                 /* 337 */
  "Pink",                 /* 338 */
  "Pink1",                /* 339 */
  "Pink2",                /* 340 */
  "Pink3",                /* 341 */
  "Pink4",                /* 342 */
  "Plum",                 /* 343 */
  "Plum1",                /* 344 */
  "Plum2",                /* 345 */
  "Plum3",                /* 346 */
  "Plum4",                /* 347 */
  "PowderBlue",           /* 348 */
  "Purple",               /* 349 */
  "Purple1",              /* 350 */
  "Purple2",              /* 351 */
  "Purple3",              /* 352 */
  "Purple4",              /* 353 */
  "Red",                  /* 354 */
  "Red1",                 /* 355 */
  "Red2",                 /* 356 */
  "Red3",                 /* 357 */
  "Red4",                 /* 358 */
  "Red5",                 /* 359 */
  "Red6",                 /* 360 */
  "Red7",                 /* 361 */
  "RosyBrown",            /* 362 */
  "RosyBrown1",           /* 363 */
  "RosyBrown2",           /* 364 */
  "RosyBrown3",           /* 365 */
  "RosyBrown4",           /* 366 */
  "RoyalBlue",            /* 367 */
  "RoyalBlue1",           /* 368 */
  "RoyalBlue2",           /* 369 */
  "RoyalBlue3",           /* 370 */
  "RoyalBlue4",           /* 371 */
  "SaddleBrown",          /* 372 */
  "Salmon",               /* 373 */
  "Salmon1",              /* 374 */
  "Salmon2",              /* 375 */
  "Salmon3",              /* 376 */
  "Salmon4",              /* 377 */
  "SandyBrown",           /* 378 */
  "SeaGreen",             /* 379 */
  "SeaGreen1",            /* 380 */
  "SeaGreen2",            /* 381 */
  "SeaGreen3",            /* 382 */
  "SeaGreen4",            /* 383 */
  "Seashell",             /* 384 */
  "Seashell1",            /* 385 */
  "Seashell2",            /* 386 */
  "Seashell3",            /* 387 */
  "Seashell4",            /* 388 */
  "Sienna",               /* 389 */
  "Sienna1",              /* 390 */
  "Sienna2",              /* 391 */
  "Sienna3",              /* 392 */
  "Sienna4",              /* 393 */
  "SkyBlue",              /* 394 */
  "SkyBlue1",             /* 395 */
  "SkyBlue2",             /* 396 */
  "SkyBlue3",             /* 397 */
  "SkyBlue4",             /* 398 */
  "SlateBlue",            /* 399 */
  "SlateBlue1",           /* 400 */
  "SlateBlue2",           /* 401 */
  "SlateBlue3",           /* 402 */
  "SlateBlue4",           /* 403 */
  "SlateGrey",            /* 404 */
  "SlateGrey1",           /* 405 */
  "SlateGrey2",           /* 406 */
  "SlateGrey3",           /* 407 */
  "SlateGrey4",           /* 408 */
  "Snow",                 /* 409 */
  "Snow1",                /* 410 */
  "Snow2",                /* 411 */
  "Snow3",                /* 412 */
  "Snow4",                /* 413 */
  "SpringGreen",          /* 414 */
  "SpringGreen1",         /* 415 */
  "SpringGreen2",         /* 416 */
  "SpringGreen3",         /* 417 */
  "SpringGreen4",         /* 418 */
  "SteelBlue",            /* 419 */
  "SteelBlue1",           /* 420 */
  "SteelBlue2",           /* 421 */
  "SteelBlue3",           /* 422 */
  "SteelBlue4",           /* 423 */
  "Tan",                  /* 424 */
  "Tan1",                 /* 425 */
  "Tan2",                 /* 426 */
  "Tan3",                 /* 427 */
  "Tan4",                 /* 428 */
  "Thistle",              /* 429 */
  "Thistle1",             /* 430 */
  "Thistle2",             /* 431 */
  "Thistle3",             /* 432 */
  "Thistle4",             /* 433 */
  "Tomato",               /* 434 */
  "Tomato1",              /* 435 */
  "Tomato2",              /* 436 */
  "Tomato3",              /* 437 */
  "Tomato4",              /* 438 */
  "Turquoise",            /* 439 */
  "Turquoise1",           /* 440 */
  "Turquoise2",           /* 441 */
  "Turquoise3",           /* 442 */
  "Turquoise4",           /* 443 */
  "Violet",               /* 444 */
  "VioletRed",            /* 445 */
  "VioletRed1",           /* 446 */
  "VioletRed2",           /* 447 */
  "VioletRed3",           /* 448 */
  "VioletRed4",           /* 449 */
  "Wheat",                /* 450 */
  "Wheat1",               /* 451 */
  "Wheat2",               /* 452 */
  "Wheat3",               /* 453 */
  "Wheat4",               /* 454 */
  "White",                /* 455 */
  "WhiteSmoke",           /* 456 */
  "Yellow",               /* 457 */
  "Yellow1",              /* 458 */
  "Yellow2",              /* 459 */
  "Yellow3",              /* 460 */
  "Yellow4",              /* 461 */
  "YellowGreen"};         /* 462 */

