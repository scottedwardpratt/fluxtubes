#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	int Amax,a1,p0=0,q0=0,pf=0,qf=0,Ntraj,itraj;
	char filename[150];
	CRandy *randy=new CRandy(-time(NULL));
	vector<vector<vector<double>>> pqcount;
	
  printf("Enter Amax: ");
  scanf("%d",&Amax);
	printf("Enter p0 and q0: ");
	scanf("%d %d",&p0,&q0);
	printf("Enter pf and qf: ");
	scanf("%d %d",&pf,&qf);
	printf("Enter Ntraj: ");
	scanf("%d",&Ntraj);
	snprintf(filename,100,"trajectories/A%d.dat",Amax);
	FILE *fptr=fopen(filename,"w");
	Ctraj traj(Amax);
	vector<double> Cquad(Amax+1),Ccubic(Amax+1);
	for(a1=0;a1<=Amax;a1++){
		Cquad[a1]=Ccubic[a1]=0.0;
	}
	
	
	CalcPQCount(pf,qf,Amax,pqcount);
	for(itraj=0;itraj<Ntraj;itraj++){
		traj.FindTrajectory(p0,q0,Amax,pqcount,randy);
		traj.CalcCasimirs();
		traj.Write(fptr);
		for(a1=0;a1<=Amax;a1++){
			Cquad[a1]+=traj.Cquad[a1];
			Ccubic[a1]+=traj.Ccubic[a1];
		}
		//traj.Print();
	}
	fclose(fptr);
	for(a1=0;a1<=Amax;a1++){
		printf("Cquad[%d]=%g, Ccubic[%d]=%g\n",a1,Cquad[a1]/double(Ntraj),a1,Ccubic[a1]/double(Ntraj));
	}

	return 0;
}
