
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2015
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

AR = ar
ARFLAGS = rv

OBJECTS = pknum01fd.o pknum02fd.o pknum03fd.o pknum04fd.o pknum05fd.o \
  pknum06fd.o pknum07fd.o pknum08fd.o pknum09fd.o pknum10fd.o pknum11fd.o \
  pknum12fd.o pknum13fd.o pknum14fd.o pknum15fd.o pknum16fd.o pknum17fd.o \
  pknum18fd.o pknum19fd.o pknum20fd.o pknum21fd.o pknum22fd.o pknum23fd.o \
  pknum24fd.o pknum25fd.o \
  pknum01f.o pknum02f.o pknum03f.o pknum04f.o pknum05f.o pknum06f.o \
  pknum07f.o pknum08f.o pknum09f.o pknum10f.o pknum11f.o pknum12f.o \
  pknum13f.o pknum14f.o pknum15f.o pknum16f.o pknum17f.o pknum18f.o \
  pknum19f.o pknum20f.o pknum21f.o pknum22f.o pknum23f.o pknum24f.o \
  pknum25f.o pknum26f.o pknum27f.o pknum28f.o pknum29f.o pknum30f.o \
  pknum31f.o pknum32f.o pknum33f.o pknum34f.o pknum35f.o pknum36f.o \
  pknum37f.o pknum38f.o pknum39f.o pknum40f.o pknum41f.o pknum42f.o \
  pknum44f.o pknum47f.o pknum48f.o pknum49f.o pknum50f.o pknum51f.o \
  pknum52f.o pknum53f.o pknum54f.o pknum55f.o pknum56f.o pknum57f.o \
  pknum58f.o pknum59f.o pknum60f.o pknum61f.o \
  pknum01d.o pknum02d.o pknum03d.o pknum04d.o pknum05d.o pknum06d.o \
  pknum07d.o pknum08d.o pknum09d.o pknum10d.o pknum11d.o pknum12d.o \
  pknum13d.o pknum14d.o pknum15d.o pknum16d.o pknum17d.o pknum18d.o \
  pknum19d.o pknum20d.o pknum21d.o pknum22d.o pknum23d.o pknum24d.o \
  pknum25d.o pknum26d.o pknum27d.o pknum28d.o pknum29d.o pknum30d.o \
  pknum31d.o pknum32d.o pknum33d.o pknum34d.o pknum35d.o pknum36d.o \
  pknum37d.o pknum38d.o pknum39d.o pknum40d.o pknum41d.o pknum42d.o \
  pknum44d.o pknum47d.o pknum48d.o pknum49d.o pknum50d.o pknum51d.o \
  pknum52d.o pknum53d.o pknum54d.o pknum55d.o pknum56d.o pknum57d.o \
  pknum58d.o pknum59d.o pknum60d.o pknum61d.o

TARGET = $(LIBDIR)/libpknum.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o: $< $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pknumf.h $(INCDIR)/pknumd.h
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)

