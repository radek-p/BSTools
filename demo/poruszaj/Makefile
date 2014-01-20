
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2011, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

LDFLAGS =

LIBS = -lxgedit -lbsfile -lbsmesh -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -lGLU -lGL -lX11 -lXext -lpthread -lm

OBJECTS = linkage.o ludzik.o glutki.o palm.o lighting.o widgets.o anima.o texture.o

default: $(BINDIR)/poruszaj

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

ludzik.o: $(SRCDIR)/ludzik.c $(SRCDIR)/ludzik.h $(SRCDIR)/palm.h \
  $(SRCDIR)/linkage.h $(SRCDIR)/glutki.h $(SRCDIR)/texture.h

linkage.o: $(SRCDIR)/linkage.c $(SRCDIR)/linkage.h

glutki.o: $(SRCDIR)/glutki.c $(SRCDIR)/glutki.h

palm.o: $(SRCDIR)/palm.c $(SRCDIR)/palm.h

lighting.o: $(SRCDIR)/lighting.c $(SRCDIR)/lighting.h

widgets.o: $(SRCDIR)/widgets.c $(SRCDIR)/widgets.h

anima.o: $(SRCDIR)/anima.c $(SRCDIR)/anima.h

texture.o: $(SRCDIR)/texture.c $(SRCDIR)/texture.h

$(BINDIR)/poruszaj: $(SRCDIR)/poruszaj.c $(SRCDIR)/poruszaj.h $(SRCDIR)/widgets.h \
  $(SRCDIR)/palm.h $(SRCDIR)/lighting.h $(SRCDIR)/anima.h $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/poruszaj.c \
  $(OBJECTS) $(LIBDIR) $(LIBS) -o $(BINDIR)/poruszaj
	rm -f $(SRCDIR)/poruszaj
	ln -s $(BINDIR)/poruszaj $(SRCDIR)/poruszaj

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/poruszaj $(SRCDIR)/poruszaj

