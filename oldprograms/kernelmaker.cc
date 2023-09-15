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
  int A,Amax,pmax,p,q,pqmax,p1,p2,q1,q2;
  double dtrue,dtot,Mtot,nsinglets,casbar=0.0,zcas=0.0;
	vector<vector<vector<vector<vector<int>>>>> pqcount;
  Tpq *pq;
  Tpqlist **pqlist;
	vector<double> ncasimir;
	CRandy *randy=new CRandy(-1234);
	char filename[100];
  printf("Enter Amax: ");
  scanf("%d",&Amax);
	pqmax=lrint(2.5*Amax);
	pqcount.resize(Amax+1);
	for(A=0;A<=Amax;A++){
		pqcount[A].resize(pqmax+1);
		for(p1=0;p1<=pqmax;p1++){
			pqcount[A][p1].resize(pqmax+1);
			for(q1=0;q1<=pqmax;q1++){
				pqcount[A][p1][q1].resize(pqmax+1);
				for(p2=0;p2<=pqmax;p2++){
					pqcount[A][p1][q1][p2].resize(pqmax+1);
					for(q2=0;q2<=pqmax;q2++)
						pqcount[A][p1][q1][p2][q2]=0;
				}
			}
		}
	}
	printf("check a\n");
	
	pqlist=new Tpqlist *[Amax+1];
	for(A=0;A<=Amax;A++){
		pqlist[A]=new Tpqlist(pqmax);
	}
	ncasimir.resize((Amax+1)*(Amax+1));
  
	for(p1=0;p1<=pqmax;p1++){
		for(q1=0;q1<=pqmax;q1++){
			printf("check b: p1=%d, q1=%d\n",p1,q1);
			for(A=0;A<=Amax;A++)
				pqlist[A]->clear();
			printf("check bb\n");
			pqlist[0]->add(p1,q1,1);
			printf("check bbb\n");
			for(A=1;A<=Amax;A++){
				for(pq=pqlist[A-1]->first;pq!=NULL;pq=pq->next){
					printf("hmmmmm, p=%d, q=%d, n=%g\n",pq->p,pq->q,pq->n);
					addmults_list(pq->p,pq->q,pq->n,1,1,1,pqlist[A]);
					printf("OK\n");
				}
			}
			printf("check c, pqmax=%d\n",pqmax);
			for(A=0;A<=Amax;A++){
				for(pq=pqlist[A]->first;pq!=NULL;pq=pq->next){
					p=pq->p; q=pq->q;
					pqcount[A][p1][q1][p][q]+=degen(p,q)*pq->n;
				}
			}
			printf("check d\n");
		}
	}
	

	return 0;
}
