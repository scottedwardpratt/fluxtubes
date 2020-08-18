#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	int Amax,p,q,pmax,a,Ntraj,itraj,Nsample=10;
	double YB=7.0, WG=0.8;  // beam rapidity and gaussian thermal smearing width
	vector<double> eta;
	CRandy *randy=new CRandy(-time(NULL));
	vector<vector<vector<double>>> pqcount;
	double edens,y;
	double Ebar=0,E2bar=0,E3bar=0,E4bar=0,kappa2,kappa3,kappa4,Ssigma_avg=0,Ksigma2_avg=0,omega_avg=0,Serror,Kerror,werror;
	vector<double> Ssigma,Ksigma2,omega;
	//printf("Enter Amax: ");
	//scanf("%d",&Amax);
	printf("Enter Ntraj: ");
	scanf("%d",&Ntraj);
	printf("Enter y: ");
	scanf("%lf",&y);
	multimap<double,int> etamap;
	multimap<double,int> amap;
	multimap<double,int>::iterator iter;
	string filename="moments.dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"Amax\tSsigma_avg\tSerror\tKsigma2_avg\tKerror\tomega_avg\twerror\n");
	for(Amax=10;Amax<=90;Amax+=10){
		CalcPQCount(Amax,pqcount);
		pmax=GetPmax(Amax+2);
		eta.resize(Amax+2);
		Ctraj traj(Amax+1);
		Ssigma.clear();
		Ksigma2.clear();
		omega.clear();
		Ssigma_avg=Ksigma2_avg=omega_avg=0;
		Serror=Kerror=werror=0;
		//vector<double> C(Amax);
		for(int isample=0;isample<Nsample;isample++){
			Ebar=E2bar=E3bar=E4bar=0;
			printf("A=%d, beginning sample %d\n",Amax,isample);
			for(int itraj=0;itraj<Ntraj;itraj++){
				etamap.clear();
				etamap.insert(pair<double,int>(-YB,0));
				etamap.insert(pair<double,int>(YB,Amax+1));
				//make sure on charge is either at +Amax or -Amax
				for(a=1;a<=Amax;a++){
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
				edens=0;
				for(a=0;a<=Amax;a++){
					edens+=.5*sqrt(traj.casimir[a])*(erf((eta[a+1]-y)/(sqrt(2)*WG))-erf((eta[a]-y)/(sqrt(2)*WG)));
				}
				Ebar+=edens; E2bar+=edens*edens; E3bar+=pow(edens,3); E4bar+=pow(edens,4);
			}
			Ebar/=double(Ntraj);
			E2bar/=double(Ntraj);
			E3bar/=double(Ntraj);
			E4bar/=double(Ntraj);

			kappa2=E2bar-Ebar*Ebar;
			kappa3=E3bar-3*kappa2*Ebar-Ebar*Ebar*Ebar;
			kappa4=E4bar-4*kappa3*Ebar-3*kappa2*kappa2-6*kappa2*Ebar*Ebar-Ebar*Ebar*Ebar*Ebar;

			Ssigma.push_back(kappa3/kappa2);
			Ksigma2.push_back(kappa4/kappa2);
			omega.push_back(kappa2/Ebar);
			Ssigma_avg+=Ssigma.back();
			Ksigma2_avg+=Ksigma2.back();
			omega_avg+=omega.back();
			//printf("omega=%lf Ssigma=%lf Ksigma2=%lf \n",omega.back(),Ssigma.back(),Ksigma2.back());
		}
		Ssigma_avg/=double(Nsample);
		Ksigma2_avg/=double(Nsample);
		omega_avg/=double(Nsample);
		for(int isample=0;isample<Nsample;isample++){
	        Serror+=(Ssigma[isample]-Ssigma_avg)*(Ssigma[isample]-Ssigma_avg);
	        Kerror+=(Ksigma2[isample]-Ksigma2_avg)*(Ksigma2[isample]-Ksigma2_avg);
			werror+=(omega[isample]-omega_avg)*(omega[isample]-omega_avg);
		}
		Serror=sqrt(Serror)/double(Nsample);
		Kerror=sqrt(Kerror)/double(Nsample);
		werror=sqrt(werror)/double(Nsample);
		fprintf(fptr,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Amax,Ssigma_avg,Serror,Ksigma2_avg,Kerror,omega_avg,werror);
		ClearPQCount(Amax,pqcount);
	}
	return 0;
}
