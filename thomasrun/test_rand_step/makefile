EIGEN_INCLUDE_PATH = /usr/include/eigen3
CC = clang++ -std=c++11 -O2

test_rand_traj : test_rand_traj.cc su3libs
	${CC} -o test_rand_traj test_rand_traj.cc -I  ../../software/include -L  ../../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/local/lib -lsu3

calc_test_weight : calc_test_weight.cc su3libs
	${CC} -o calc_test_weight calc_test_weight.cc -I  ../../software/include -L  ../../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/local/lib -lsu3

su3libs : 
	make -C  ../../software install

