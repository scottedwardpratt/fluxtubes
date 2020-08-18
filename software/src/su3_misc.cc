#ifndef __ADDMULTS_CC__
#define __ADDMULTS_CC__

#include "su3.h"
using namespace std;
using namespace NS_SU3;

void NS_SU3::ReadNpq(vector<vector<vector<double>>> &Npq,int Amax,vector<vector<double>> &biggestweight){
	int pmax,a,p,q,a1,a2,d;
	double weight;
	char filename[150];
	FILE *fptr;
	biggestweight.resize(Amax+1);
	for(a=0;a<=Amax;a++){
		biggestweight[a].resize(Amax+1);
	}
	pmax=GetPmax(Amax);
	Npq.resize(Amax+1);
	for(a=0;a<=Amax;a++){
		Npq[a].resize(pmax+1);
		for(p=0;p<=pmax;p++){
			Npq[a][p].resize(p+1);
			for(q=0;q<=p;q++){
				Npq[a][p][q]=0;
			}
		}
		if(a>0){
			sprintf(filename,"opentrajectories/Npq_A%d.dat",a);
			fptr=fopen(filename,"r");
			fscanf(fptr,"%d %d",&p,&q);
			while(!feof(fptr)){
				fscanf(fptr,"%lf",&Npq[a][p][q]);
				fscanf(fptr,"%d %d",&p,&q);
			}
			fclose(fptr);
		}
	}
	for(a1=0;a1<=Amax;a1++){
		for(a2=0;a2<=a1;a2++){
			biggestweight[a1][a2]=0.0;
			if(a1!=a2)
				biggestweight[a2][a1]=0.0;
			if(a1<a2)
				pmax=GetPmax(a1);
			else
				pmax=GetPmax(a2);
			for(p=0;p<=pmax;p++){
				for(q=0;q<=p;q++){
					d=degen(p,q);
					if(a1==a2){
						//weight=0.5*Npq[a1][p][q]*(Npq[a2][p][q]-1)/(d*d);
						weight=0.5*Npq[a1][p][q]*(Npq[a2][p][q]-1.0)/(d*d);
						if(p==q)
							weight*=4.0;
					}
					else{
						weight=0.5*Npq[a1][p][q]*Npq[a2][p][q]/(d*d);
						if(p==q)
							weight*=4.0;
					}
					
					if(weight>biggestweight[a1][a2]){
						biggestweight[a1][a2]=weight;
						if(a1!=a2)
							biggestweight[a2][a1]=biggestweight[a1][a2];
					}
				}
			}
		}
	}
}

void NS_SU3::WriteOpenTrajectories(int A,long long int ntraj,CRandy *randy){
	long long int itraj;
	int pmax=GetPmax(A),p,q,a;	
	double w;
	char filename[150];

	vector<int> ptraj;
	vector<int> qtraj;
	ptraj.resize(A+1);
	qtraj.resize(A+1);

	// Npq counts the number of trajectories for given pq
	vector<vector<double>> Npq;
	Npq.resize(pmax+1);
	for(p=0;p<=pmax;p++){
		Npq[p].resize(p+1);
		for(q=0;q<=p;q++)
			Npq[p][q]=0;
	}
	// Open Files to write trajectoies
	vector<vector<FILE *>> opentraj;
	opentraj.resize(pmax+1);
	for(p=0;p<=pmax;p++){
		opentraj[p].resize(p+1);
		for(q=0;q<=p;q++){
				sprintf(filename,"opentrajectories/A%d/p%d_q%d.dat",A,p,q);
				opentraj[p][q]=fopen(filename,"w");
		}
	}
	
	// NOW LET'S MAKE TRAJECTORIES
	for(itraj=0;itraj<ntraj;itraj++){
		ptraj[0]=qtraj[0]=0;
		p=q=0;
		for(a=1;a<=A;a++){
			RanStep(p,q,randy->ran(),w);
			if(p>GetPmax(a) || q>GetPmax(a)){
				printf("pq is too big=%d,%d, a=%d",p,q,a);
			}
			ptraj[a]=p;
			qtraj[a]=q;
		}
		if(p>pmax || q>pmax){
			printf("%d:%d, pmax=%d\n",p,q,pmax);
			exit(1);
		}
		for(a=0;a<=A;a++){
			//printf("A=%d, p=%d, q=%d, ptraj=%d,qtraj=%d\n",A,p,q,ptraj[A],qtraj[A]);
			if(p>=q){
				fprintf(opentraj[p][q],"%d %d ",ptraj[a],qtraj[a]);
			}
			else{
				fprintf(opentraj[q][p],"%d %d ",ptraj[a],qtraj[a]);
			}
			//printf("--------\n");
		}
		if(p>=q){
			fprintf(opentraj[p][q],"\n");
			Npq[p][q]+=1;
		}
		else{
			fprintf(opentraj[q][p],"\n");
			Npq[q][p]+=1;
		}
		if(((itraj+1)*10)%ntraj==0)
			printf("finished %g percent of trajectories\n",double(itraj+1)*100.0/double(ntraj));
	}
	for(p=0;p<=pmax;p++){
		for(q=0;q<=p;q++){
			fclose(opentraj[p][q]);
		}
	}
	sprintf(filename,"opentrajectories/Npq_A%d.dat",A);
	FILE *fptr=fopen(filename,"w");
	for(p=0;p<=pmax;p++){
		for(q=0;q<=p;q++){
			fprintf(fptr,"%d %d %g\n",p,q,Npq[p][q]);
		}
	}
	fclose(fptr);
}

void NS_SU3::CalcPQCount(int Amax,vector<vector<vector<double>>> &pqcount){
  int A,pmax,p,q;
  int nquarks=0,nanti=0,ngluons,casimir;
	double dtot,Mtot,nsinglets=0;
	double dtrue;
  Tpq *pq;
  Tpqlist **pqlist;
	CRandy *randy=new CRandy(-1234);
	ngluons=Amax;
  //Amax=ngluons+nquarks+nanti;
  pqlist=new Tpqlist *[Amax+1];
	pmax=GetPmax(ngluons)+nquarks+nanti;
  for(A=0;A<=Amax;A++)
		pqlist[A]=new Tpqlist(GetPmax(ngluons)+nquarks+nanti);
  pqlist[0]->add(0,0,1);
	
  for(A=1;A<=Amax;A++){
    pqlist[A]->clear();
    for(pq=pqlist[A-1]->first;pq!=NULL;pq=pq->next){
      if(A<=ngluons)
				addmults_list(pq->p,pq->q,pq->n,1,1,1,pqlist[A]);
      if(A>ngluons && A<=ngluons+nquarks)
				addmults_list(pq->p,pq->q,pq->n,1,0,1,pqlist[A]);
      if(A>ngluons+nquarks)
				addmults_list(pq->p,pq->q,pq->n,0,1,1,pqlist[A]);
    }
		//pqlist[A-1]->clear();
		//pqlist[A]->compress();
  }
	pqcount.resize(Amax+1);
	for(A=0;A<=Amax;A++){
		pqcount[A].resize(pmax+1);
		for(p=0;p<=pmax;p++){
			pqcount[A][p].resize(pmax+1,0);
		}
	}	
	for(A=0;A<=Amax;A++){
		dtot=Mtot=0.0;
		for(pq=pqlist[A]->first;pq!=NULL;pq=pq->next){
			p=pq->p; q=pq->q;
			pqcount[A][p][q]+=degen(p,q)*pq->n;
			if(A==Amax){
				dtot+=degen(p,q)*pq->n;
				Mtot+=pq->n;
				if(A==Amax && pq->p==0 && pq->q==0)
					nsinglets+=pq->n;
			}
		}	
		if(A==Amax){
			dtrue=pow(double(3),nquarks+nanti)*pow(8.0,ngluons);
			printf("---- A=%d, dtot=%g =? %g ----\n",A,dtot,dtrue);
			printf("Nsinglets=%g = %g fraction of multiplets and %g fraction of states\n",nsinglets,nsinglets/Mtot,nsinglets/dtot);
		}
	}
	for(A=0;A<=Amax;A++){
		pqlist[A]->clear();
	}
}

void NS_SU3::ClearPQCount(int Amax,vector<vector<vector<double>>> &pqcount){
	int pmax=GetPmax(Amax);
	for(int A=0;A<=Amax;A++){
		for(int p=0;p<=pmax;p++){
			for(int q=0;q<=pmax;q++){
				pqcount[A][p][q]=0;
			}
		}
	}
}

double NS_SU3::omega_massless(double T,double V){
	const double prefactor=(1.0/(2.0*PI*PI*HBARC*HBARC*HBARC));
	return 4.0*V*prefactor*T*T*T;
}

double NS_SU3::omegae_massless(double T,double V){
	const double prefactor=(1.0/(2.0*PI*PI*HBARC*HBARC*HBARC));
	return 12.0*V*prefactor*T*T*T*T;
}

double NS_SU3::omega_bessel(double mass,double T,double V){
	const double prefactor=(1.0/(PI*PI*HBARC*HBARC*HBARC));
	double z,answer;
	z=mass/T;
	answer=V*pow(mass,3)*prefactor*(z*cyl_bessel_k(0,z)+2.0*cyl_bessel_k(1,z))/(z*z);
	return answer;
}

double NS_SU3::omegae_bessel(double mass,double T,double V){
	const double prefactor=(1.0/(PI*PI*HBARC*HBARC*HBARC));
	double z,answer;
	z=mass/T;
	answer=V*pow(mass,4)*prefactor*(3.0*z*cyl_bessel_k(0,z)+(z*z+6.0)*cyl_bessel_k(1,z))/(z*z*z);
	return answer;
}

double NS_SU3::omegap_bessel(double mass,double T,double V){
	const double prefactor=(1.0/(PI*PI*HBARC*HBARC*HBARC));
	double z,answer;
	z=mass/T;
	answer=V*mass*mass*T*T*prefactor*(cyl_bessel_k(0,z)+2*cyl_bessel_k(1,z)/z);
	return answer;
}

int NS_SU3::GetPmax(int A){
	return lrint(0.1+1.5*A);
}

int NS_SU3::degen(int p,int q){
	int answer;
	answer=(p+1)*(q+1)*(p+q+2);
	return answer/2;
}

int NS_SU3::Casimir(int p,int q){
	return (p*p+q*q+3*p+3*q+p*q)/3;
}

#endif
