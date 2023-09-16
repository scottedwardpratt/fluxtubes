#include <cmath>
#include <cstdio>
#include <complex>
#include "msu_commonutils/misc.h"

const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

using namespace std;

#include "subs/misc.cc"
#include "subs/pq.cc"
#include "subs/tableaux.cc"
#include "subs/addmults.cc"

int main(){
	int p,q;
	double C1,C2,x,y;
	printf("Enter p,q: ");
	scanf("%d %d",&p,&q);
	Casimir12(p,q,C1,C2);
	x=p+q; y=p-q;
	
	printf("C1=%g, C2=%g\n",C1,C2);
	printf("x=%g, y=%g\n",x,y);
	

	double a=1.0,b=0.0,c=-9.0*C1-9.0,d=18.0*C2;
	printf("c=%g, d=%g\n",c,d);
	
	complex<double> z1,z2,z3;
	double y1,y2,y3,x1=-999,x2=-999,x3=-999,rarg;
	Misc::CubicComplex(d,c,b,a,z1,z2,z3);
	
	y1=real(z1);
	rarg=4.0-y1*y1/3.0+4*C1;
	if(rarg>0)
		x1=-2.0+sqrt(rarg);
	
	y2=real(z2);
	rarg=4.0-y2*y2/3.0+4*C1;
	if(rarg>0)
		x2=-2.0+sqrt(rarg);
	
	y3=real(z3);
	rarg=4.0-y3*y3/3.0+4*C1;
	if(rarg>0)
		x3=-2.0+sqrt(rarg);
	
	if(x1>-0.001 && fabs(imag(z1))<0.0001 && x1-fabs(y1)>-0.0001){
		printf("y1=(%g,%g), x1=%g\n",y1,imag(z1),x1);
		printf("------- TEST1:  p=%g, q=%g\n",0.5*(x1+y1),0.5*(x1-y1));
	}
	if(x2>-0.001 && fabs(imag(z2))<0.0001 && x2-fabs(y2)>-0.0001){
		printf("y2=(%g,%g), x2=%g\n",y2,imag(z2),x2);
		printf("------- TEST2:  p=%g, q=%g\n",0.5*(x2+y2),0.5*(x2-y2));
	}
	if(x3>-0.001 && fabs(imag(z3))<0.0001 && x3-fabs(y3)>-0.0001){
		printf("y3=(%g,%g), x3=%g\n",y3,imag(z3),x3);
		printf("------- TEST3:  p=%g, q=%g\n",0.5*(x3+y3),0.5*(x3-y3));
	}
	
	double t,arg1,arg2;
	for(int k=0;k<3;k++){
		arg1=(1.5*d/c)*sqrt(-3.0/c);
		arg2=acos(arg1)/3.0-k*2.0*PI/3.0;
		t=2.0*sqrt(-c/3.0)*cos(arg2);
		if(fabs(t-y)<.0001)
			printf("t_%d=%g, test=%g\n",k,t,t*t*t+c*t+d);
	}
	
	
}