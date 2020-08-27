EIGEN_INCLUDE_PATH = /usr/local/include/eigen3
CC = clang++ -std=c++11 -O2

all :
	make -C ../software
	make edens
	make edens_vsA
	make edist

edens : edens.cc ../software/lib/libsu3.a
	${CC} -o edens edens.cc -I  ../software/include -L  ../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/usr/local/lib -lsu3

edens_vsA : edens_vsA.cc ../software/lib/libsu3.a
	${CC} -o edens_vsA edens_vsA.cc -I  ../software/include -L  ../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/usr/local/lib -lsu3

edist : edist.cc ../software/lib/libsu3.a
	${CC} -o edist edist.cc -I  ../software/include -L  ../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/usr/local/lib -lsu3

traj : traj.cc ../software/lib/libsu3.a
	${CC} -o traj traj.cc -I  ../software/include -L  ../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/usr/local/lib -lsu3

corr : corr.cc ../software/lib/libsu3.a
	${CC} -o corr corr.cc -I  ../software/include -L  ../software/lib -I/usr/local/include -L/usr/local/lib -lsu3

slope_smeared : slope_smeared.cc ../software/lib/libsu3.a
	${CC} -o slope_smeared slope_smeared.cc -I  ../software/include -L  ../software/lib -I/usr/local/include -L/usr/local/lib -lsu3

slope : slope.cc ../software/lib/libsu3.a
	${CC} -o slope slope.cc -I  ../software/include -L  ../software/lib -I/usr/local/include -L/usr/local/lib -lsu3

exact : exact.cc ../software/lib/libsu3.a
	${CC} -o exact exact.cc -I  ../software/include -L  ../software/lib -I/usr/local/include -L/usr/local/lib -lsu3

findmults_list : findmults_list.cc ../software/lib/libsu3.a
	${CC} -o findmults_list findmults_list.cc -I  ../software/include -L  ../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/usr/local/lib -lsu3
