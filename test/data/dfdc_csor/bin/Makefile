SRC1 = ../src-lib
SRC2 = ../src-driver

## Basic compiler setup
AR = ar

##--------------------------
## OSX gfortran
FC = gfortran
FFLAGS = -O -std=legacy
FFLOPT = -O
LIBS =  -L/usr/X11R6/lib -lX11
FTNLIB =
LFLAGS = -static-libgfortran -static-libgcc

##--------------------------

OBJSLIB = hello.o, csor.o, test.o

PROGS = test

#===========================================================

all: $(PROGS)

install:
	$(INSTALLCMD) $(PROGS) $(BINDIR)

clean:
	-/bin/rm *.o
	-/bin/rm $(PROGS)


cleanall:
	-/bin/rm *.o
	-/bin/rm *.a
	-/bin/rm $(PROGS)

#===========================================================

hello.o: $(SRC2)/hello.f
	$(FC) -c $(FFLAGS) -I$(SRC2) $(SRC2)/hello.f

csor.o: $(SRC1)/csor.f
	$(FC) -c $(FFLAGS) -I$(SRC1) $(SRC1)/csor.f

test.o: $(SRC2)/test.f $(FC) -c $(FFLAGS) -I$(SRC1) $(SRC2)/test.f

#===========================================================
