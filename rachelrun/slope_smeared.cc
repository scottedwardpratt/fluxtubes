#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	int Amax,p,q,pmax,a,a1,a2,Ntraj,itraj,dela,abar,ahalf;
	char filename[150];
	double w,eta_temp,alpha,alpha2,etabar;
	double YB=7.0, WG=0.8;  // beam rapidity and gaussian thermal smearing width
	vector<double> eta;
	CRandy *randy=new CRandy(-time(NULL));
	vector<vector<vector<double>>> pqcount;
	double deleta=2.0;
	FILE *fptr=fopen(filename,"w");
	double eta1=-deleta,eta2=deleta,rho1,rho2,E,Etot,Etotbar;
	printf("Enter Amax: ");
	scanf("%d",&Amax);
	ahalf=lrint(0.5*Amax);
	pmax=GetPmax(Amax+2);
	eta.resize(Amax+2);
	CalcPQCount(Amax,pqcount);
	multimap<double,int> etamap;
	multimap<double,int> amap;
	multimap<double,int>::iterator iter;
	Ctraj traj(Amax+1);
	vector<double> C(Amax);
	printf("Enter Ntraj: ");
	scanf("%d",&Ntraj);
	
	for(int isample=0;isample<10;isample++){
		alpha2=Etotbar=0.0;
		for(int itraj=0;itraj<Ntraj;itraj++){
			etamap.clear();
			etamap.insert(pair<double,int>(-YB,0));
			etamap.insert(pair<double,int>(YB,Amax+1));
			//make sure on charge is either at +Amax or -Amax
			if(randy->ran()<0.5){
				etamap.insert(pair<double,int>(-YB+0.0001,Amax));
			}
			else{
				etamap.insert(pair<double,int>(YB-0.0001,Amax));
			}
			for(a=1;a<Amax;a++){
				eta_temp=YB*(1.0-2.0*randy->ran());
				etamap.insert(pair<double,int>(eta_temp,a));
			}
			a=0;
			for(iter=etamap.begin();iter!=etamap.end();++iter){
				eta[a]=iter->first;
				//printf("%2d: %g\n",a,eta[a]);
				a=a+1;
			}
			traj.FindTrajectory(Amax+1,pqcount,randy);
			traj.CalcCasimirs();
			//traj.PrintCasimirs();
			rho1=rho2=Etot=0.0;
			for(a=1;a<Amax;a++){
				E=traj.casimir[a]*(eta[a+1]-eta[a]);
				etabar=0.5*(eta[a+1]+eta[1]);
				rho1+=E*exp(-0.5*pow(etabar-eta1,2)/(WG*WG));
				rho2+=E*exp(-0.5*pow(etabar-eta2,2)/(WG*WG));
				Etot+=E;
			}
			Etotbar+=Etot;
			//printf("itraj=%d, eta1=%g, eta2=%g, rho1=%g, rho2=%g\n",itraj,eta1,eta2,rho1,rho2);
			alpha=(rho1-rho2)/(Etot);
			alpha2+=alpha*alpha;		
		}
		alpha2=alpha2/double(Ntraj);
		Etotbar=Etotbar/double(Ntraj);
		printf("<alpha^2>=%g, <Etot>=%g, <alpha^2*Etot>=%g\n",alpha2,Etotbar,alpha2*Etotbar);
	}
	
	return 0;
}
