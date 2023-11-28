MADAI_GSLPATH = /usr/local
MADAI_CPP = /usr/bin/clang++
#MADAI_CPP = /usr/local/bin/clang-omp++
#MADAI_CPP = /usr/bin/g++
#compiler

MADAI_CFLAGS = -Oz -std=c++14 -Wall

INCLUDE = -I include -I /usr/local/include -I /usr/local/include/eigen

LIBRARY = -L lib -L /usr/local/lib

OBJFILES = build/addmults.o build/pq.o build/randy.o build/su3_misc.o build/tableaux.o build/trajectory.o

lib/libsu3.a : ${OBJFILES}
	ar -ru lib/libsu3.a ${OBJFILES}

build/addmults.o : include/su3.h src/addmults.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/addmults.cc -o build/addmults.o

build/pq.o : include/su3.h src/pq.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/pq.cc -o build/pq.o

build/randy.o : include/randy.h src/randy.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/randy.cc -o build/randy.o

build/su3_misc.o : include/su3.h src/su3_misc.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/su3_misc.cc -o build/su3_misc.o

build/tableaux.o : include/su3.h src/tableaux.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/tableaux.cc -o build/tableaux.o

build/trajectory.o : include/su3.h src/trajectory.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INCLUDE} src/trajectory.cc -o build/trajectory.o

clean:
	rm lib/*.a build/*.o
