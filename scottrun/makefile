MADAI_GSLPATH = /usr/local
MADAI_CPP = /usr/bin/clang++
#MADAI_CPP = /usr/local/bin/clang-omp++
#MADAI_CPP = /usr/bin/g++
#compiler

MADAI_CFLAGS = -Oz -std=c++14 -Wall

INCLUDE = -I ../software/include -I /usr/local/include -I /usr/local/include/eigen -I /usr/local/include/Eigen

LIBRARY = -L ../software/lib -L /usr/local/lib

../software/lib/libsu3.a :
	make -B -C ../software lib/libsu3.a

traj : ../software/include/su3.h traj.cc ../software/lib/libsu3.a
	make -C ../software
	${MADAI_CPP} ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lsu3 traj.cc -o traj
	
slope : ../software/include/su3.h slope.cc ../software/lib/libsu3.a
	make -C ../software
	${MADAI_CPP} ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lsu3 slope.cc -o slope
	
findmults_list : ../software/include/su3.h findmults_list.cc ../software/lib/libsu3.a
	make -C ../software
	${MADAI_CPP} ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lsu3 findmults_list.cc -o findmults_list
	
quark_traj : quark_traj.cc ../software/lib/libsu3.a
	make -C ../software lib/libsu3.a
	${MADAI_CPP} ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lsu3 quark_traj.cc -o quark_traj
	
bcorr : bcorr.cc ../software/lib/libsu3.a
	make -C ../software lib/libsu3.a
	${MADAI_CPP} ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lsu3 bcorr.cc -o bcorr
