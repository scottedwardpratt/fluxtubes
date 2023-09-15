#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <cmath>
#include "subs/randy.h"

using namespace std;

#include "subs/misc.cc"
#include "subs/pq.cc"
#include "subs/tableaux.cc"
#include "subs/addmults.cc"
#include "subs/randy.cc"

class Ctraj{
public:
	vector<int> ptraj,qtraj;
	int A;
	Ctraj(int Amax){
		ptraj.resize(Amax+1);
		qtraj.resize(Amax+1);
	}
};

class Cweight{
public:
	double weight,ransum;
	int p,q;
};

int main(){
	CRandy *randy=new CRandy(-time(NULL));
	double ransum=0.0,ranthresh=0.0;
	int Amax,success=0;
	printf("Enter Amax: ");
	scanf("%d",&Amax);
	bool reading;
	int A,pmax=int(Amax*1.5),p,q;
	int itraj,jtraj,ntraj,ipq;
	double w,wtraj,weight,CZ=0.0;
	vector<double> cbar;
	vector<Cweight *> weightlist;
	double wtot=0.0;
	Cweight *cw;
	vector<int> qtraj,ptraj,casimir;
	cbar.resize(2*Amax+1);
	ptraj.resize(2*Amax+1);
	qtraj.resize(2*Amax+1);
	casimir.resize(2*Amax+1);
	FILE *fptr;
	char filename[100];
	vector<Ctraj *> traj;  
	vector<vector<long long int>> Npq;
	Npq.resize(pmax+1);
	for(p=0;p<=pmax;p++){
		Npq[p].resize(p+1);
		for(q=0;q<=p;q++)
			Npq[p][q]=0;
	}
	sprintf(filename,"halftrajdata/Npq_A%d.dat",Amax);
	fptr=fopen(filename,"r");
	int pp,qq;
	for(p=0;p<=pmax;p++){
		for(q=0;q<=p;q++){
			fscanf(fptr,"%d %d %lld\n",&pp,&qq,&Npq[p][q]);
			printf("Npq[%d][%d]=%lld\n",p,q,Npq[p][q]);
			if(Npq[p][q]>1){
				cw=new Cweight();
				cw->p=p; cw->q=q;
				cw->weight=(double(Npq[p][q])*double(Npq[p][q]-1)/2)/(degen(p,q)*degen(p,q));
				if(p!=q) cw->weight=cw->weight*0.5;
				wtot+=cw->weight;
				weightlist.push_back(cw);
			}
		}
	}
	fclose(fptr);
	for(int i=0;i<weightlist.size();i++){
		cw=weightlist[i];
		cw->weight=cw->weight/wtot;
		ransum+=cw->weight;
		cw->ransum=ransum;
	}
	
	// NOW LET'S READ AND ANALYZE TRAJECTORIES
	
	vector<double> corr(Amax+1,0.0);
	int delA;
	
	for(A=0;A<=2*Amax;A++)
		cbar[A]=0.0;
	
	printf("Enter Ntraj: ");
	scanf("%d",&ntraj);

	ransum=0.0;
	ranthresh=-log(randy->ran())/double(ntraj);
	for(ipq=0;ipq<weightlist.size();ipq++){
		cw=weightlist[ipq];
		p=cw->p; q=cw->q;
		ransum+=cw->weight;
		reading=false;
		while(ransum>ranthresh){
			if(reading==false){
				sprintf(filename,"halftrajdata/A%d_p%d_q%d.dat",Amax,p,q);
				printf("reading %s, Npq=%lld\n",filename,Npq[p][q]);
				fptr=fopen(filename,"r");
				if(traj.size()<Npq[p][q]){
					itraj=traj.size();
					traj.resize(Npq[p][q]);
					while(itraj<Npq[p][q]){
						traj[itraj]=new Ctraj(Amax);
						itraj+=1;
					}
				}
				for(itraj=0;itraj<Npq[p][q];itraj++){
					for(A=0;A<=Amax;A++){
						fscanf(fptr,"%d %d",&traj[itraj]->ptraj[A],&traj[itraj]->qtraj[A]);
					}
				}
				fclose(fptr);
				reading=true;
			}
			CZ+=1.0;
			do{
				itraj=floorl(randy->ran()*double(Npq[p][q]));
				jtraj=floorl(randy->ran()*double(Npq[p][q]));
			}while(itraj==jtraj);
			for(A=0;A<=Amax;A++){
				ptraj[A]=traj[itraj]->ptraj[A];
				qtraj[A]=traj[itraj]->qtraj[A];
			}
			for(A=0;A<Amax;A++){
				ptraj[2*Amax-A]=traj[jtraj]->ptraj[A];
				qtraj[2*Amax-A]=traj[jtraj]->qtraj[A];
			}
			for(A=0;A<=2*Amax;A++){
				casimir[A]=Casimir(ptraj[A],qtraj[A]);
				cbar[A]+=casimir[A];
			}
			for(delA=0;delA<=Amax;delA++){
				corr[delA]+=casimir[Amax+delA]*casimir[Amax-delA];
			}
			ranthresh-=log(randy->ran())/double(ntraj);
			success+=1;
		}
	}
		
	sprintf(filename,"figs/Cbar_A%d.dat",Amax);
	fptr=fopen(filename,"w");
	for(A=0;A<=2*Amax;A++){
		cbar[A]=cbar[A]/CZ;
		printf("%3d %g\n",A,cbar[A]);
		fprintf(fptr,"%2d %g\n",A,cbar[A]);
	}
	fclose(fptr);
	sprintf(filename,"figs/corr_A%d.dat",Amax);
	fptr=fopen(filename,"w");
	for(delA=0;delA<Amax;delA++){
		corr[delA]=corr[delA]/(CZ*cbar[Amax+delA]*cbar[Amax-delA]);
		printf("corr[%2d]=%15.12f\n",2*delA,corr[delA]);
		fprintf(fptr,"%2d %g\n",2*delA,corr[delA]);
	}
	fclose(fptr);
	printf("success=%d, success/ntraj=%g\n",success,double(success)/double(ntraj));
}
