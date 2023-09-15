#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cmath>
#include "subs/randy.h"

using namespace std;

#include "subs/misc.cc"
#include "subs/pq.cc"
#include "subs/tableaux.cc"
#include "subs/addmults.cc"
#include "subs/randy.cc"

int main(){
  int A,Amax,dela,p,q,ntraj,itraj,A1,A2,acmin,acmax;
	double weight,wtot,cw;
	char filename[100];
	printf("Enter Amax and ntraj\n");
	scanf("%d %d",&Amax,&ntraj);
	acmin=Amax/2 -10; // ranges for correlation
	acmax=Amax/2 +10;
	vector<int> ptraj(Amax+1);
	vector<int> qtraj(Amax+1);
	vector<double> denom(Amax+1,0.0);
	vector<double> cbar(Amax+1,0.0);
	vector<long long int> c(Amax+1);
	vector<double> corr(21);
	vector<double> cdenom(21);
	wtot=0.0;
	for(dela=0;dela<=20;dela++){
		cdenom[dela]=corr[dela]=0.0;
	}
	sprintf(filename,"trajectories/Amax%d_ntraj%d.dat",Amax,ntraj);
	ntraj=ntraj/1000;
	FILE *fptr=fopen(filename,"r");
	for(itraj=0;itraj<ntraj;itraj++){
		fscanf(fptr,"%lf",&weight);
		if(weight<0){
			printf("negative weight\n");
			exit(1);
		}
		wtot+=weight;
		for(A=0;A<=Amax;A++){
			fscanf(fptr,"%d %d",&ptraj[A],&qtraj[A]);
			c[A]=Casimir(ptraj[A],qtraj[A]);
			cbar[A]+=c[A]*weight;
			denom[A]+=weight;
		}	
		for(A1=acmin;A1<acmax;A1++){
			for(A2=acmin+1;A2<=acmax;A2++){
				corr[A2-A1]+=weight*c[A1]*c[A2];
			}
		}
	}
	fclose(fptr);
	
	for(A=0;A<=Amax;A++){
		cbar[A]=cbar[A]/denom[A];
		printf("Cbar(%d)=%g\n",A,cbar[A]);
	}
	for(A1=acmin;A1<acmax;A1++){
		for(A2=acmin+1;A2<=acmax;A2++){
			cdenom[A2-A1]+=cbar[A1]*cbar[A2];
		}
	}
	printf("----- Correlations -------\n");
	for(dela=1;dela<=20;dela++){
		corr[dela]=corr[dela]/wtot;
		printf("%2d  %12.4e %12.4e %12.4e\n",
		dela,corr[dela],cdenom[dela],corr[dela]/cdenom[dela]);
	}
	
	return 0;
}
