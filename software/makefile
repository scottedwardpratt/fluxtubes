include makefile_defs.mk

OBJFILES = build/randy.o build/trajectory.o build/tableaux.o build/pq.o build/addmults.o build/su3_misc.o
HFILES = include/constants.h include/su3.h include/randy.h include/defs.h 

install: lib/libsu3.a

build/randy.o : src/randy.cc ${HFILES}
	${CC} -c -o build/randy.o src/randy.cc -I include -I${EIGEN_INCLUDE_PATH}

build/trajectory.o : src/trajectory.cc ${HFILES}
	${CC} -c -o build/trajectory.o src/trajectory.cc -I include -I${EIGEN_INCLUDE_PATH}

build/tableaux.o : src/tableaux.cc ${HFILES}
	${CC} -c -o build/tableaux.o src/tableaux.cc -I include -I${EIGEN_INCLUDE_PATH}

build/pq.o : src/pq.cc ${HFILES}
	${CC} -c -o build/pq.o src/pq.cc -I include -I${EIGEN_INCLUDE_PATH}

build/addmults.o : src/addmults.cc ${HFILES}
	${CC} -c -o build/addmults.o src/addmults.cc -I include -I${EIGEN_INCLUDE_PATH}

build/su3_misc.o : src/su3_misc.cc ${HFILES}
	${CC} -c -o build/su3_misc.o src/su3_misc.cc -I include -I${EIGEN_INCLUDE_PATH}

lib/libsu3.a : ${OBJFILES}
	ar -ru lib/libsu3.a ${OBJFILES}

clean:
	rm -f ${OBJFILES}

uninstall:
	rm -f lib/libsu3.a
