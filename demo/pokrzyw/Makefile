
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2009
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

LIBS = -lxgedit -lbsfile -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -lX11 -lpthread -lm

OBJECTS = $(SRCDIR)/spline2d.o $(SRCDIR)/ed2dwidgets.o

default: $(BINDIR)/pokrzyw

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(BINDIR)/pokrzyw: $(SRCDIR)/ed2dspl.c $(SRCDIR)/ed2dspl.h \
  $(SRCDIR)/ed2dwidgets.h $(SRCDIR)/spline2d.h $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/ed2dspl.c $(OBJECTS) \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/pokrzyw
	rm -f $(SRCDIR)/pokrzyw
	ln -s $(BINDIR)/pokrzyw $(SRCDIR)/pokrzyw

$(SRCDIR)/spline2d.o: $(SRCDIR)/spline2d.c $(SRCDIR)/spline2d.h

$(SRCDIR)/ed2dwidgets.o: $(SRCDIR)/ed2dwidgets.c $(SRCDIR)/ed2dwidgets.h

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/pokrzyw $(SRCDIR)/pokrzyw

