#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	int A,Amax=0,p0=0,q0=0,pf=0,qf=0,Ntraj,itraj,Btarget,Bproj,Ngluons,Nquarks,Nanti;
	CRandy *randy=new CRandy(-time(NULL));
	vector<vector<vector<double>>> pqcount;
	vector<Tpqlist *> pqlist0,pqlistf;
	double dtot,r,dsum;
	Tpq *pq;
	Ctraj *traj;	
 
	printf("Enter B of target: ");
	scanf("%d",&Btarget);
	printf("Enter B of projectile: ");
	scanf("%d",&Bproj);
	printf("Enter Ngluons: ");
	scanf("%d",&Ngluons);
	printf("Enter Ntraj: ");
	scanf("%d",&Ntraj);
	
	Nquarks=3*Btarget;
	Nanti=3*Bproj;
	Amax=Nquarks+Nanti+Ngluons;
	pqlist0.resize(Amax+1);
	for(A=0;A<=Amax;A++){
		pqlist0[A]=new Tpqlist(Amax+1);
		pqlist0[A]->clear();
	}
	
	pqlist0[0]->add(0,0,1);
	for(A=1;A<=Nquarks;A++){
		for(pq=pqlist0[A-1]->first; pq!= NULL; pq=pq->next){
			addquark_list(pq->p,pq->q,pq->n,pqlist0[A]);
		}
	}
	for(A=Nquarks+1;A<=Nquarks+Ngluons;A++){
		for(pq=pqlist0[A-1]->first; pq!= NULL; pq=pq->next){
			addgluon_list(pq->p,pq->q,pq->n,pqlist0[A]);
		}
	}
	for(A=Nquarks+Ngluons+1;A<=Amax;A++){
		for(pq=pqlist0[A-1]->first; pq!= NULL; pq=pq->next){
			addantiquark_list(pq->p,pq->q,pq->n,pqlist0[A]);
		}
	}
	pqlist0[Amax]->print();
	printf("-------------------\n");
	
	for(itraj=0;itraj<Ntraj;itraj++){
		
		dtot = 0.0;
		for (pq = pqlist0[Amax]->first; pq != NULL; pq = pq->next){
			dtot+= degen(pq->p,pq->q)*pq->n;
		}
		r=randy->ran();
		pq = pqlist0[Amax]->first;
		dsum=degen(pq->p,pq->q)*pq->n;
		p0=pq->p;
		q0=pq->q;
		while(pq!=NULL && r>dsum/dtot){
			pq=pq->next;
			if(pq!=NULL){
				p0=pq->p;
				q0=pq->q;
				dsum+=degen(pq->p,pq->q)*pq->n;
			}
		}while(pq!=NULL && r>dsum/dtot);
		printf("p0=%d, q0=%d\n",p0,q0);
	
		dtot = 0.0;
		for (pq = pqlistf[Amax]->first; pq != NULL; pq = pq->next){
			dtot+= degen(pq->p,pq->q)*pq->n;
		}
		r=randy->ran();
		pq = pqlistf[Amax]->first;
		dsum=degen(pq->p,pq->q)*pq->n;
		p0=pq->p;
		q0=pq->q;
		while(pq!=NULL && r<dsum/dtot){
			pq=pq->next;
			if(pq!=NULL){
				p0=pq->p;
				q0=pq->q;
				dsum+=degen(pq->p,pq->q)*pq->n;
			}
		}while(pq!=NULL && r<dsum/dtot);
		
		traj=new Ctraj(Ngluons);
		traj->FindOpenTrajectory(p0,q0,Ngluons,randy);
		//if(traj->ptraj[Ngluons]==pf && traj->qtraj[Ngluons]==qf)
			traj->Print();
		
	}
	
	/*
	
	
	printf("Enter Ngluons: ");
	scanf("%d ",&Ngluons);
	Amax = ngluons + p0 + qf;
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
	*/
	return 0;
}
