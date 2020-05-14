EIGEN_INCLUDE_PATH = /usr/include/eigen3
CC = clang++ -std=c++11 -O2

gluon_walker : gluon_walker.cc su3libs
	${CC} -o gluon_walker gluon_walker.cc -I  ../../software/include -L  ../../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/local/lib -lsu3

su3libs : 
	make -C  ../../software install

