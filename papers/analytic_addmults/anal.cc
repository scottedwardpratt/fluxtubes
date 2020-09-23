#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <boost/math/special_functions.hpp>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;
//using namespace boost::math;

long long int degen(int p,int q){
	return (p+1)*(q+1)*(p+q+2)/2;
}
long long int Casimir(int p,int q){
	return (p*p+q*q+3*p+3*q+p*q)/3;
}
double nmults(int p,int q,double sigma){
	double a,b,g,theta,x,y,r,sigma2=sigma*sigma,M=1.5;
	a=sigma*sqrt(3.0);
	b=sigma;
	x=(p+q+1.0)/a;
	y=(p-q)/b;
	r=sqrt(x*x+y*y);
	if(r<40.0*sigma){
		//printf("x=%g, y=%g\n",x,y);
		g=(1.0/sigma)*pow(r/sigma2,M)*exp(-0.5*r*r);
		theta=atan2(y,x);
		return g;
		//return g*cos(M*theta);
	}
	else
		return 0.0;
}

int main(int argc,char *argv[]){
	const int AMAX=4000;
	int p,q;
	double x,y,r,theta;
	double nm,sigma,norm,Ebar;
	for(sigma=0.001*AMAX;sigma<=0.0200001*AMAX;sigma+=0.001*AMAX){
		norm=0.0;
		Ebar=0.0;
		for(p=0;p<=AMAX;p++){
			for(q=0;q<=AMAX;q++){
				nm=nmults(p,q,sigma);
				norm+=nm*degen(p,q);
				Ebar+=nm*degen(p,q)*Casimir(p,q);
			}
		}
		Ebar=Ebar/norm;
		printf("AMAX=%3d, sigma=%g, Ebar=%g, Ebar/sigma^4=%g\n",AMAX,sigma,Ebar,Ebar/pow(sigma,2));
	}
	return 0;
}


