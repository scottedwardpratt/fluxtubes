#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cmath>
#include "subs/randy.h"

using namespace std;

#include "subs/misc.cc"
#include "subs/pq.cc"
#include "subs/tableaux.cc"
#include "subs/addmults.cc"
#include "subs/randy.cc"

int main(){
  int A,Amax,pmax,p,q;
  int nquarks=0,nanti=0,ngluons=6,casimir;
  double dtrue,dtot,Mtot,nsinglets,casbar=0.0,zcas=0.0;
  Tpq *pq;
  Tpqlist **pqlist;
	vector<double> ncasimir;
	CRandy *randy=new CRandy(-1234);
	char filename[100];
		
  printf("Enter Amax: ");
  scanf("%d",&Amax);
	ngluons=Amax;
  //Amax=ngluons+nquarks+nanti;
	ncasimir.resize((Amax+1)*(Amax+1));
  pqlist=new Tpqlist *[Amax+1];
	pmax=int(1.5*ngluons)+nquarks+nanti;
	//pmax=int(1.5*Amax);
  for(A=0;A<=Amax;A++)
		pqlist[A]=new Tpqlist(int(1.5*ngluons)+nquarks+nanti);
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
			//pqlist[A]->compress();
  }
	
  nsinglets=0;
	vector<vector<vector<double>>> pqcount;
	pqcount.resize(Amax+1);
	for(A=0;A<=Amax;A++){
		pqcount[A].resize(pmax+1);
		for(p=0;p<=pmax;p++){
			pqcount[A][p].resize(pmax+1,0.0);
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
				casimir=Casimir(p,q);
				if(ncasimir.size()<=casimir){
					printf("ncasimir length too short, casimir=%d\n",casimir);
					exit(1);
				}
				ncasimir[casimir]+=pq->n*dtot;		
				if(A==Amax && pq->p==0 && pq->q==0)
					nsinglets+=pq->n;
			}
		}
	  dtrue=pow(double(3),nquarks+nanti)*pow(8.0,ngluons);
	  printf("---- A=%d, dtot=%g =? %g ----\n",A,dtot,dtrue);
	  printf("Nsinglets=%g = %g percent of multiplets and %g percent of states\n",nsinglets,nsinglets/Mtot,nsinglets/dtot);
	}
	
 
		//for(A=0;A<=Amax;A++) pqlist[A]->print();
	pqlist[Amax]->print();
	printf("-------------------------\n");
	for(int icasimir=0;icasimir<ncasimir.size();icasimir++){
		zcas+=ncasimir[icasimir];
		casbar+=ncasimir[icasimir]*icasimir;
		//printf("ncasimir[%d]=%g\n",icasimir,ncasimir[icasimir]);
	}
	casbar=casbar/zcas;
	printf("<Casimir>=%g\n",casbar);
	
	int A1,A2;
	long long int NA=0,N12=0,N1,N2;
	double Cbar,CZbar,w;
	sprintf(filename,"Cbar_Amax%d.dat",Amax);
	FILE *fptr=fopen(filename,"w");
	NA=pqcount[Amax][0][0];
	for(A1=0;A1<=Amax;A1++){
		A2=Amax-A1;
		N12=N1=N2=0;
		Cbar=CZbar=0.0;
		for(p=0;p<=pmax;p++){
			for(q=0;q<=pmax;q++){
				N1+=pqcount[A1][p][q];
				N2+=pqcount[A2][p][q];
				w=pqcount[A1][p][q]*pqcount[A2][p][q]/(degen(p,q)*degen(p,q));
				N12+=w;
				CZbar+=w;
				Cbar+=w*Casimir(p,q);
			}
		}
		fprintf(fptr,"%3d %g\n",A1,Cbar/CZbar);
		printf("A1=%3d  Cbar=%g\n",A1,Cbar/CZbar);
		//printf("nsinglets for A=%d=%lld\n",A,NA);
		//printf("N1=%lld=?%lld, N2=%lld=?%lld\n",N1,(long long int)pow(8,A1),N2,(long long int)pow(8,A2));
		//printf("N12=%lld, NA/N12=%g\n",N12,double(NA)/double(N12));
	}
	fclose(fptr);

	return 0;
}
