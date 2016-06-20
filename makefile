CC := g++
nagcdir :=/n/home13/gobied/viniciuscode/setup/cll6i25dcl/naglib/cll6i25dcl
nagclink= -L${nagcdir}/lib/ -lnagc_nag -L${nagcdir}/rtl/intel64/ -lifcoremt -lsvml -limf -lirc -lirng -lintlc 
CFLAGS := -pg -g -std=gnu++11 -lgsl -lgslcblas -I${nagcdir}/include/ -DNAGLIB $(nagclink)

all: main.o

main.o: main.cpp bkg.o muk.o approx.o parameters.o libbkg.a libmuk.a libapprox.a libparameters.a
	$(CC) $(CFLAGS) main.cpp -o main.o -L. -lbkg -lmuk -lapprox -lparameters 

libbkg.a: bkg.o
	ar rcs libbkg.a bkg.o

libmuk.a: muk.o
	ar rcs libmuk.a muk.o

libapprox.a: approx.o
	ar rcs libapprox.a approx.o

libparameters.a: parameters.o
	ar rcs libparameters.a parameters.o

bkg.o: bkg.cpp bkg.h parameters.h utils.h constants.h
	$(CC) $(CFLAGS) -c bkg.cpp

muk.o: muk.cpp muk.h bkg.h parameters.h utils.h constants.h
	$(CC) $(CFLAGS) -c muk.cpp

approx.o: approx.cpp approx.h utils.h constants.h
	$(CC) $(CFLAGS) -c approx.cpp

parameters.o: parameters.cpp parameters.h
	$(CC) $(CFLAGS) -c parameters.cpp

clean:
	$(RM) main.o *o *~


