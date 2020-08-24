#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	int Amax,p,q,pmax,a,Ntraj,itraj,Nsample=10;
	double YB=7.0, WG=0.8;  // beam rapidity and gaussian thermal smearing width
	vector<double> eta;
	CRandy *randy=new CRandy(-time(NULL));
	vector<vector<vector<double>>> pqcount;
	double edens,y,eta_temp;
	double Ebar=0,E2bar=0,E3bar=0,E4bar=0,kappa1_avg=0,kappa2_avg=0,kappa3_avg=0,kappa4_avg=0,Ssigma_avg=0,Ksigma2_avg=0,c3c1_avg=0,c4c3_avg=0;
	double k1error,k2error,k3error,k4error,omega_avg=0,Serror,Kerror,werror,c3c1error,c4c3error;
	vector<double> kappa1,kappa2,kappa3,kappa4,Ssigma,Ksigma2,omega,c3c1,c4c3;
	printf("Enter Ntraj: ");
	scanf("%d",&Ntraj);
	printf("Enter y: ");
	scanf("%lf",&y);
	multimap<double,int> etamap;
	multimap<double,int> amap;
	multimap<double,int>::iterator iter;
	string filename="moments_vsA.dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"y\tkappa1\tk1error\tkappa2\tk2error\tkappa3\tk3error\tkappa4\tk4error\tSsigma_avg\tSerror\tKsigma2_avg\tKerror\tomega_avg\twerror\n");
	for(Amax=10;Amax<=100;Amax+=10){
		CalcPQCount(Amax,pqcount);
		pmax=GetPmax(Amax+2);
		eta.resize(Amax+2);
		Ctraj traj(Amax+1);
		Ssigma.clear(); Ksigma2.clear(); omega.clear();
		kappa1.clear(); kappa2.clear(); kappa3.clear(); kappa4.clear();
		c3c1.clear(); c4c3.clear();
		Ssigma_avg=Ksigma2_avg=omega_avg=0;
		kappa1_avg=kappa2_avg=kappa3_avg=kappa4_avg=0;
		c3c1_avg=c4c3_avg=0;
		Serror=Kerror=werror=0;
		k1error=k2error=k3error=k4error=0;
		c3c1error=c4c3error=0;
		printf("------------------------------------------------------\n");
		for(int isample=0;isample<Nsample;isample++){
			Ebar=E2bar=E3bar=E4bar=0;
			printf("Amax=%d, beginning sample %d\n",Amax,isample);
			for(int itraj=0;itraj<Ntraj;itraj++){
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
				edens=0;
				for(a=0;a<=Amax;a++){
					edens+=.5*traj.casimir[a]*(erf((eta[a+1]-y)/(sqrt(2)*WG))-erf((eta[a]-y)/(sqrt(2)*WG)));
				}
				Ebar+=edens; E2bar+=edens*edens; E3bar+=pow(edens,3); E4bar+=pow(edens,4);
			}
			Ebar/=double(Ntraj);
			E2bar/=double(Ntraj);
			E3bar/=double(Ntraj);
			E4bar/=double(Ntraj);

			kappa1.push_back(Ebar);
			kappa2.push_back(E2bar-Ebar*Ebar);
			kappa3.push_back(E3bar-3*kappa2.back()*Ebar-Ebar*Ebar*Ebar);
			kappa4.push_back(E4bar-4*kappa3.back()*Ebar-3*kappa2.back()*kappa2.back()-6*kappa2.back()*Ebar*Ebar-Ebar*Ebar*Ebar*Ebar);

			Ssigma.push_back(kappa3.back()/kappa2.back());
			Ksigma2.push_back(kappa4.back()/kappa2.back());
			omega.push_back(kappa2.back()/kappa1.back());
			c3c1.push_back(kappa3.back()/kappa1.back());
			c4c3.push_back(kappa4.back()/kappa3.back());

			Ssigma_avg+=Ssigma.back(); Ksigma2_avg+=Ksigma2.back(); omega_avg+=omega.back();
			kappa1_avg+=kappa1.back(); kappa2_avg+=kappa2.back(); kappa3_avg+=kappa3.back(); kappa4_avg+=kappa4.back();
			c3c1_avg+=c3c1.back(); c4c3_avg+=c4c3.back();
		}
		//printf("kappa1_avg=%lf kappa2_avg=%lf kappa3_avg=%lf kappa4_avg=%lf\n",kappa1_avg,kappa2_avg,kappa3_avg,kappa4_avg);
		Ssigma_avg/=double(Nsample); Ksigma2_avg/=double(Nsample); omega_avg/=double(Nsample);
		kappa1_avg/=double(Nsample); kappa2_avg/=double(Nsample); kappa3_avg/=double(Nsample); kappa4_avg/=double(Nsample);
		c3c1_avg/=double(Nsample); c4c3_avg/=double(Nsample);
		for(int isample=0;isample<Nsample;isample++){
	        Serror+=(Ssigma[isample]-Ssigma_avg)*(Ssigma[isample]-Ssigma_avg);
	        Kerror+=(Ksigma2[isample]-Ksigma2_avg)*(Ksigma2[isample]-Ksigma2_avg);
			werror+=(omega[isample]-omega_avg)*(omega[isample]-omega_avg);
			k1error+=(kappa1[isample]-kappa1_avg)*(kappa1[isample]-kappa1_avg);
			k2error+=(kappa2[isample]-kappa2_avg)*(kappa2[isample]-kappa2_avg);
			k3error+=(kappa3[isample]-kappa3_avg)*(kappa3[isample]-kappa3_avg);
			k4error+=(kappa4[isample]-kappa4_avg)*(kappa4[isample]-kappa4_avg);
			c3c1error+=(c3c1[isample]-c3c1_avg)*(c3c1[isample]-c3c1_avg);
			c4c3error+=(c4c3[isample]-c4c3_avg)*(c4c3[isample]-c4c3_avg);
		}
		Serror=sqrt(Serror)/double(Nsample); Kerror=sqrt(Kerror)/double(Nsample); werror=sqrt(werror)/double(Nsample);
		k1error=sqrt(k1error)/double(Nsample); k2error=sqrt(k2error)/double(Nsample); k3error=sqrt(k3error)/double(Nsample); k4error=sqrt(k4error)/double(Nsample);
		c3c1error=sqrt(c3c1error)/double(Nsample); c4c3error=sqrt(c4c3error)/double(Nsample);
		fprintf(fptr,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
					Amax,kappa1_avg,k1error,kappa2_avg,k2error,kappa3_avg,k3error,kappa4_avg,k4error,
					omega_avg,werror,Ssigma_avg,Serror,Ksigma2_avg,Kerror,c3c1_avg,c3c1error,c4c3_avg,c4c3error);
		ClearPQCount(Amax,pqcount);
	}
	delete randy;
	return 0;
}
