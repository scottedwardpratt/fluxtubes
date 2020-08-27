#include "su3.h"
#define NBINS 50
using namespace NS_SU3;
using namespace std;

int main(){
	int Amax,p,q,pmax,a,Ntraj,itraj,Nsample=10,n;
	double YB=7.0, WG=0.8;  // beam rapidity and gaussian thermal smearing width
	vector<double> eta;
	CRandy *randy=new CRandy(-time(NULL));
	vector<vector<vector<double>>> pqcount;
	double edens=0,y=0,eta_temp;
	printf("Enter Ntraj: ");
	scanf("%d",&Ntraj);
	printf("Enter n: ");
	scanf("%d",&n);
	printf("Enter Amax: ");
	scanf("%d",&Amax);
	CalcPQCount(Amax,pqcount);
	pmax=GetPmax(Amax+2);
	eta.resize(Amax+2);
	Ctraj traj(Amax+1);
	multimap<double,int> etamap;
	multimap<double,int> amap;
	multimap<double,int>::iterator iter;
	string filename="edist.dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"edens\tfrequency\n");
	int ibin;
	double Emax=100;
	double dE=Emax/double(NBINS);
	int edist[NBINS]={};
	printf("------------------------------------------------------\n");
	for(int isample=0;isample<Nsample;isample++){
		printf("y=%lf, beginning sample %d\n",y,isample);
		for(int i=0; i<Ntraj/n; i++) {
			for(int itraj=0;itraj<n;itraj++){
				etamap.clear();
				etamap.insert(pair<double,int>(-YB,0));
				etamap.insert(pair<double,int>(YB,Amax+1));
				for(a=1;a<=Amax;a++){
					eta_temp=YB*(1.0-2.0*randy->ran());
					etamap.insert(pair<double,int>(eta_temp,a));
				}
				a=0;
				for(iter=etamap.begin();iter!=etamap.end();++iter){
					eta[a]=iter->first;
					a=a+1;
				}
				traj.FindTrajectory(Amax+1,pqcount,randy);
				traj.CalcCasimirs();
				//traj.PrintCasimirs();
				for(a=0;a<=Amax;a++){
					edens+=.5*traj.casimir[a]*(erf((eta[a+1]-y)/(sqrt(2)*WG))-erf((eta[a]-y)/(sqrt(2)*WG)));
				}
			}
			ibin=floorl((edens/double(n))/dE);
            edist[ibin]+=1;
			//printf("edens=%lf\n",edens/double(n));
			edens=0;
		}
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fprintf(fptr,"%lf\t%lf\n",ibin*dE+dE/2.,edist[ibin]/double(Nsample*Ntraj/double(n)));
		printf("%lf\t%lf\n",ibin*dE+dE/2.,edist[ibin]/double(Nsample*Ntraj/double(n)));
	}
	delete randy;
	ClearPQCount(Amax,pqcount);
	fclose(fptr);
	return 0;
}
