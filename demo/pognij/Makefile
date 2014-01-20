
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

#LDFLAGS = -static
LDFLAGS =

LIBS = -lxgedit -lbsfile -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -lX11 -lXext -lpthread -lm

OBJECTS = $(SRCDIR)/spline3d.o $(SRCDIR)/ed3dwidgets.o

default: $(BINDIR)/pognij

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/spline3d.o: $(SRCDIR)/spline3d.c $(SRCDIR)/spline3d.h

$(SRCDIR)/ed3dwidgets.o: $(SRCDIR)/ed3dwidgets.c $(SRCDIR)/ed3dwidgets.h

$(BINDIR)/pognij: $(SRCDIR)/ed3dspl.c $(SRCDIR)/ed3dspl.h \
  $(SRCDIR)/spline3d.h $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/ed3dspl.c $(OBJECTS) \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/pognij
	rm -f $(SRCDIR)/pognij
	ln -s $(BINDIR)/pognij $(SRCDIR)/pognij

clean:
	rm -f $(SRCDIR)/*.o

mrproper:
	rm -f $(SRCDIR)/*.o
	rm -f $(BINDIR)/pognij $(SRCDIR)/pognij

