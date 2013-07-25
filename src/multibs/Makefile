
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
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = mbs01fd.o \
  mbs01f.o mbs02f.o mbs03f.o mbs04f.o mbs05f.o mbs06f.o mbs07f.o mbs08f.o \
  mbs09f.o mbs10f.o mbs11f.o mbs12f.o mbs13f.o mbs14f.o mbs15f.o mbs16f.o \
  mbs17f.o mbs18f.o mbs19f.o mbs20f.o mbs21f.o mbs22f.o mbs23f.o mbs24f.o \
  mbs25f.o mbs26f.o mbs27f.o mbs28f.o mbs29f.o mbs30f.o mbs31f.o mbs32f.o \
  mbs33f.o mbs34f.o mbs35f.o mbs36f.o mbs37f.o mbs38f.o mbs40f.o mbs41f.o \
  mbs42f.o mbs43f.o mbs44f.o mbs45f.o mbs46f.o mbs47f.o mbs48f.o mbs49f.o \
  mbs50f.o mbs51f.o mbs52f.o mbs53f.o mbs54f.o mbs55f.o mbs56f.o mbs57f.o \
  mbs58f.o mbs59f.o mbs60f.o mbs61f.o mbs62f.o \
  mbs01d.o mbs02d.o mbs03d.o mbs04d.o mbs05d.o mbs06d.o mbs07d.o mbs08d.o \
  mbs09d.o mbs10d.o mbs11d.o mbs12d.o mbs13d.o mbs14d.o mbs15d.o mbs16d.o \
  mbs17d.o mbs18d.o mbs19d.o mbs20d.o mbs21d.o mbs22d.o mbs23d.o mbs24d.o \
  mbs25d.o mbs26d.o mbs27d.o mbs28d.o mbs29d.o mbs30d.o mbs31d.o mbs32d.o \
  mbs33d.o mbs34d.o mbs35d.o mbs36d.o mbs37d.o mbs38d.o mbs40d.o mbs41d.o \
  mbs42d.o mbs43d.o mbs44d.o mbs45d.o mbs46d.o mbs47d.o mbs48d.o mbs49d.o \
  mbs50d.o mbs51d.o mbs52d.o mbs53d.o mbs54d.o mbs55d.o mbs56d.o mbs57d.o \
  mbs58d.o mbs59d.o mbs60d.o mbs61d.o mbs62d.o

TARGET = $(LIBDIR)/libmultibs.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o: $< $(SRCDIR)/msgpool.h $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h \
  $(INCDIR)/multibs.h $(INCDIR)/multibsf.h $(INCDIR)/multibsd.h
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)
