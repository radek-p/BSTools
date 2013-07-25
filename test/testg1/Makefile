
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2007, 2008
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib

include ../../options.mak

OBJECTSF = datagenf.o testgraphf.o

OBJECTSD = datagend.o testgraphd.o

LIBS = -leghole -lmultibs -lcamera -lpsout -lpkgeom -lpknum -lpkvaria -lm

default: test0f test1f test2f test3f test4f test5f test6f \
  test7f test8f test9f test10f test11f test12f test13f test14f test15f \
  test0d test1d test2d test3d test4d test5d test6d \
  test7d test8d test9d test10d test11d test12d test13d test14d test15d

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

datagenf.o: datagenf.c datagenf.h

datagend.o: datagend.c datagend.h

drawitf.o: drawitf.c drawitf.h bslapf.h

drawitd.o: drawitd.c drawitd.h bslapd.h

bslapf.o: bslapf.c bslapf.h

bslapd.o: bslapd.c bslapd.h

testgraphf.o: testgraphf.c testgraphf.h

testgraphd.o: testgraphd.c testgraphd.h

readdatd.o: readdatd.c readdatd.h

test0f: test0f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test0f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test0f
	./test0f

test0d: test0d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test0d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test0d
	./test0d

test1f: test1f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test1f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test1f
	./test1f

test1d: test1d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test1d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test1d
	./test1d

test2f: test2f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test2f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test2f
	./test2f

test2d: test2d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test2d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test2d
	./test2d

test3f: test3f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test3f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test3f
	./test3f

test3d: test3d.c datagend.o testgraphd.o \
  readdatd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test3d.c $(OBJECTSD) readdatd.o \
  $(LIBDIR) $(LIBS) -o test3d
	./test3d

test4f: test4f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test4f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test4f
	./test4f

test4d: test4d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test4d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test4d
	./test4d

test5f: test5f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test5f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test5f
	./test5f

test5d: test5d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test5d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test5d
	./test5d

test6f: test6f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test6f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test6f
	./test6f

test6d: test6d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test6d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test6d
	./test6d

test7f: test7f.c $(OBJECTSF) bslapf.o drawitf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test7f.c $(OBJECTSF) bslapf.o drawitf.o \
  $(LIBDIR) $(LIBS) -o test7f
	./test7f

test7d: test7d.c $(OBJECTSD) bslapd.o drawitd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test7d.c $(OBJECTSD) bslapd.o drawitd.o \
  $(LIBDIR) $(LIBS) -o test7d
	./test7d

test8f: test8f.c $(OBJECTSF) writdatf.o bslapf.o drawitf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test8f.c $(OBJECTSF) \
  writdatf.o bslapf.o drawitf.o $(LIBDIR) $(LIBS) -o test8f
	./test8f

test8d: test8d.c $(OBJECTSD) writdatd.o bslapd.o drawitd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test8d.c $(OBJECTSD) \
  writdatd.o bslapd.o drawitd.o $(LIBDIR) $(LIBS) -o test8d
	./test8d

test9f: test9f.c $(OBJECTSF) bslapf.o drawitf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test9f.c $(OBJECTSF) bslapf.o drawitf.o \
  $(LIBDIR) $(LIBS) -o test9f
	./test9f

test9d: test9d.c $(OBJECTSD) bslapd.o drawitd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test9d.c $(OBJECTSD) bslapd.o drawitd.o \
  $(LIBDIR) $(LIBS) -o test9d
	./test9d

test10f: test10f.c $(OBJECTSF) bslapf.o drawitf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test10f.c $(OBJECTSF) bslapf.o drawitf.o \
  $(LIBDIR) $(LIBS) -o test10f
	./test10f

test10d: test10d.c $(OBJECTSD) bslapd.o drawitd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test10d.c $(OBJECTSD) bslapd.o drawitd.o \
  $(LIBDIR) $(LIBS) -o test10d
	./test10d

test11f: test11f.c $(OBJECTSF) bslapf.o drawitf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test11f.c $(OBJECTSF) bslapf.o drawitf.o \
  $(LIBDIR) $(LIBS) -o test11f
	./test11f

test11d: test11d.c $(OBJECTSD) bslapd.o drawitd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test11d.c $(OBJECTSD) bslapd.o drawitd.o \
  $(LIBDIR) $(LIBS) -o test11d
	./test11d

test12f: test12f.c $(OBJECTSF) bslapf.o drawitf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test12f.c $(OBJECTSF) bslapf.o drawitf.o \
  $(LIBDIR) $(LIBS) -o test12f
	./test12f

test12d: test12d.c $(OBJECTSD) bslapd.o drawitd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test12d.c $(OBJECTSD) bslapd.o drawitd.o \
  $(LIBDIR) $(LIBS) -o test12d
	./test12d

test13f: test13f.c $(OBJECTSF) bslapf.o drawitf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test13f.c $(OBJECTSF) bslapf.o drawitf.o \
  $(LIBDIR) $(LIBS) -o test13f
	./test13f

test13d: test13d.c $(OBJECTSD) bslapd.o drawitd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test13d.c $(OBJECTSD) bslapd.o drawitd.o \
  $(LIBDIR) $(LIBS) -o test13d
	./test13d

test14f: test14f.c $(OBJECTSF) bslapf.o drawitf.o writdatf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test14f.c $(OBJECTSF) \
  bslapf.o drawitf.o writdatf.o $(LIBDIR) $(LIBS) -o test14f
	./test14f

test14d: test14d.c $(OBJECTSD) bslapd.o drawitd.o writdatd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test14d.c $(OBJECTSD) \
  bslapd.o drawitd.o writdatd.o $(LIBDIR) $(LIBS) -o test14d
	./test14d

test15f: test15f.c $(OBJECTSF) datagenf.new.o bslapf.o drawitf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test15f.c $(OBJECTSF) \
  datagenf.new.o bslapf.o drawitf.o $(LIBDIR) $(LIBS) -o test15f
	./test15f

test15d: test15d.c $(OBJECTSD) datagend.new.o bslapd.o drawitd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test15d.c $(OBJECTSD) \
  datagend.new.o bslapd.o drawitd.o $(LIBDIR) $(LIBS) -o test15d
	./test15d

clean:
	rm -f *.o
	rm -f *.ps
	rm -f test?f test?d test??f test??d

mrproper:
	rm -f *.o
	rm -f *.ps
	rm -f test?f test?d test??f test??d

