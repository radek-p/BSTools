
CC = gcc

SRCDIR = .

INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

LIBS = -lbsfile -lcamera -lmengerc -lraybez -lmultibs \
  -lpkgeom -lpknum -lpkvaria -lm


TARGET = optmenger

default: $(TARGET)

optmenger: optmenger.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) optmenger.c $(LIBS) -o optmenger
	./optmenger -i c23_39.bs -e 4 -o3 -o c23_39_4.bs

gen: gen.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) gen.c $(LIBS) -o gen

clean:
	rm -f *.o *.log

mrproper:
	rm -f $(TARGET) ./gen *.o *.log

