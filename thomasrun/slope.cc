#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	int Amax,p,q;
	int Ntraj,itraj,a,ahalf,dela,a1,a2;
	double alpha,alpha2=0.0;
	CRandy *randy=new CRandy(-time(NULL));
	vector<vector<vector<double>>> pqcount;
	printf("Enter Amax: ");
	scanf("%d",&Amax);
	CalcPQCount(Amax,pqcount);
	dela=Amax/5;	
	ahalf=floorl((Amax+1.0)/2.0);
	a1=ahalf-dela;
	a2=ahalf+dela;
	Ctraj traj(Amax);
	vector<double> corr(ahalf,0.0);
	vector<double> Cbar(Amax+1);
	printf("Enter Ntraj: ");
	scanf("%d",&Ntraj);
	for(int itraj=0;itraj<Ntraj;itraj++){
		do{
			traj.FindTrajectory(Amax,pqcount,randy);
			traj.CalcCasimirs();
		}while(traj.casimir[ahalf]==0);
		for(a=0;a<=Amax;a++){
			Cbar[a]+=traj.casimir[a];
		}
		alpha=(traj.casimir[a2]-traj.casimir[a1])/(2.0*dela);
		if(traj.casimir[ahalf]!=0)
			alpha2+=alpha*alpha;
		
	}
	for(a=0;a<=Amax;a++){
		Cbar[a]=Cbar[a]/double(Ntraj);
		printf("Cbar[%d]=%g\n",a,Cbar[a]);
	}
	alpha2=alpha2/(double(Ntraj)*pow(Cbar[ahalf],0.0));
	printf("<alpha^2=%g>\n",alpha2);
	
	return 0;
}
