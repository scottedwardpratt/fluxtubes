EIGEN_INCLUDE_PATH = /usr/include/eigen3
CC = clang++ -std=c++11 -O2

traj : traj.cc su3libs
	${CC} -o traj traj.cc -I  ../software/include -L  ../software/lib -I${EIGEN_INCLUDE_PATH} -L/usr/local/lib -lsu3

corr : corr.cc su3libs
	${CC} -o corr corr.cc -I  ../software/include -L  ../software/lib -I/usr/local/include -L/usr/local/lib -lsu3

slope_smeared : slope_smeared.cc su3libs
	${CC} -o slope_smeared slope_smeared.cc -I  ../software/include -L  ../software/lib -I/usr/local/include -L/usr/local/lib -lsu3
	
slope : slope.cc su3libs
	${CC} -o slope slope.cc -I  ../software/include -L  ../software/lib -I/usr/local/include -L/usr/local/lib -lsu3

exact : exact.cc su3libs
	${CC} -o exact exact.cc -I  ../software/include -L  ../software/lib -I/usr/local/include -L/usr/local/lib -lsu3

findmults_list : findmults_list.cc su3libs
	${CC} -o findmults_list findmults_list.cc -I  ../software/include -L  ../software/lib -I${EIGEN_INCLUDE_PATH} -I/usr/local/include -L/usr/local/lib -lsu3
	
su3libs :
	make -C  ../software install
