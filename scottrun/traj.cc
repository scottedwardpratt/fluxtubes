#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	int Amax,p,q,a1,a2,pmax;
	char filename[150];
	double CZbar,Cbar,w;
	CRandy *randy=new CRandy(-time(NULL));
	vector<vector<vector<double>>> pqcount;
  printf("Enter Amax: ");
  scanf("%d",&Amax);
	pmax=GetPmax(Amax);
	CalcPQCount(Amax,pqcount);
	
	Ctraj traj(Amax);
	int Ntraj,itraj,a;
	vector<double> C(Amax+1);
	printf("Enter Ntraj: ");
	scanf("%d",&Ntraj);
	sprintf(filename,"trajectories/A%d.dat",Amax);
	FILE *fptr=fopen(filename,"w");
	for(int itraj=0;itraj<Ntraj;itraj++){
		traj.FindTrajectory(Amax,pqcount,randy);
		traj.CalcCasimirs();
		traj.Write(fptr);
		for(a1=0;a1<=Amax;a1++){
			C[a1]+=traj.casimir[a1];
		}
		traj.Print();
	}
	fclose(fptr);
	for(a1=0;a1<=Amax;a1++){
		printf("Cbar[%d]=%g\n",a1,C[a1]/double(Ntraj));
	}
	
	return 0;
}
