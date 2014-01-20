
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
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = scratchmem.o sqangle.o pksort.o pksort2.o arraysc.o writearray.o \
  signf.o signd.o exchange.o selectfd.o zeromat.o reverse.o rpower.o \
  sort2f.o sort2d.o transpose.o drawline.o errmsg.o dmalloc.o hexbyte.o \
  quicksort.o heap.o binom.o degtext.o queue.o tictoc.o pkvthreads.o \
  findncpu.o multicounter.o

TARGET = $(LIBDIR)/libpkvaria.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

convh.o: convh.c $(INCDIR)/pkvaria.h $(INCDIR)/msgpool.h $(INCDIR)/convh.h \
  $(INCDIR)/geom.h
scratchmem.o: scratchmem.c pkvprivate.h $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
sqangle.o: sqangle.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
pksort.o: pksort.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
pksort2.o: pksort2.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
arraysc.o: arraysc.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
writearray.o: writearray.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
signf.o: signf.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
signd.o: signd.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
exchange.o: exchange.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
selectfd.o: selectfd.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
zeromat.o: zeromat.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
rpower.o: rpower.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
sort2f.o: sort2f.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
sort2d.o: sort2d.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
transpose.o: transpose.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
drawline.o: drawline.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
dmalloc.o: dmalloc.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
hexbyte.o: hexbyte.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
quicksort.o: quicksort.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
heap.o: heap.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
binom.o: binom.c $(INCDIR)/pkvaria.h
degtext.o: degtext.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
queue.o: queue.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
reverse.o: reverse.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
errmsg.o: errmsg.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
tictoc.o: tictoc.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h
multicounter.o: multicounter.c $(INCDIR)/pkvaria.h
pkvthreads.o: pkvthreads.c pkvprivate.h \
  $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h $(INCDIR)/pkvthreads.h
findncpu.o: findncpu.c $(INCDIR)/msgpool.h $(INCDIR)/pkvaria.h \
  $(INCDIR)/pkvthreads.h

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)
