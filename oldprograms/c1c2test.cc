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
	int p,q,xx,yy;
	double C1,C2,x,y,pp,qq;
	double t,arg1,arg2;
	//complex<double> z1,z2,z3;
	//double y1,y2,y3,x1,x2,x3;
	//double a=1.0,b=0.0;
	double c,d,rarg;
	
	//printf("Enter p,q: ");
	//scanf("%d %d",&p,&q);
	for(xx=0;xx<20;xx++){
		printf("----------- p+q=%d ------------\n",xx);
		for(yy=xx;yy>=-xx;yy--){
			p=(xx+yy)/2;
			q=(xx-yy)/2;
			Casimir12(p,q,C1,C2);
			//x=p+q; y=p-q;
	
			//printf("C1=%g, C2=%g\n",C1,C2);
			//printf("x=%g, y=%g\n",x,y);
	

			c=-9.0*C1-9.0;
			d=18.0*C2;
			//printf("c=%g, d=%g\n",c,d);
			/*
			x1=-999;
			x2=-999;
			x3=-999;
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
			printf("------- TEST1:  p=%g =? %d, q=%g =? %d\n",0.5*(x1+y1),p,0.5*(x1-y1),q);
			}
			if(x2>-0.001 && fabs(imag(z2))<0.0001 && x2-fabs(y2)>-0.0001){
			printf("y2=(%g,%g), x2=%g\n",y2,imag(z2),x2);
			printf("------- TEST1:  p=%g =? %d, q=%g =? %d\n",0.5*(x2+y2),p,0.5*(x2-y2),q);
			}
			if(x3>-0.001 && fabs(imag(z3))<0.0001 && x3-fabs(y3)>-0.0001){
			printf("y3=(%g,%g), x3=%g\n",y3,imag(z3),x3);
			printf("------- TEST1:  p=%g =? %d, q=%g =? %d\n",0.5*(x3+y3),p,0.5*(x3-y3),q);
			}*/
	
			for(int k=0;k<3;k++){
				arg1=(1.5*d/c)*sqrt(-3.0/c);
				arg2=acos(arg1)/3.0-k*2.0*PI/3.0;
				y=2.0*sqrt(-c/3.0)*cos(arg2);
				rarg=4.0-y*y/3.0+4*C1;
				if(rarg>-0.00001){
					x=-2.0+sqrt(rarg);
					if(x>-0.0001 && x-fabs(y)>-0.0001){
						if(k!=1 || fabs(y*y*y+c*y+d)>0.001){
							printf("huh?????? k=%d, test=%g\n",k,y*y*y+c*y+d);
							exit(1);
						}
						pp=0.5*(x+y);
						qq=0.5*(x-y);
						if(fabs(pp)<0.001)
							pp=0;
						if(fabs(qq)<0.001)
							qq=0;
						if(fabs(C1)<0.001)
							C1=0;
						if(fabs(C2)<0.001)
							C2=0;
						printf("p=%g =? %d, q=%g =? %d,   C1=%10.6f, 3*C2=%10.6f\n",pp,p,qq,q,C1,3*C2);
					}
				}
			}
		}
		
	}
	
	
}