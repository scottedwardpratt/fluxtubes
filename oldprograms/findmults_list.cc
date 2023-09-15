#include <cstdlib>
#include <cmath>
#include <cstdio>

const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

using namespace std;

#include "subs/misc.cc"
#include "subs/pq.cc"
#include "subs/tableaux.cc"
#include "subs/addmults.cc"

int main(){
  int A,Amax;
  int nquarks,nanti,ngluons;
  double dtrue,dtot,nsinglets;
  Tpq *pq;
  Tpqlist **pqlist;
	
  printf("Enter Ngluons, Nquarks and Nantiquarks : ");
  scanf("%d %d %d",&ngluons,&nquarks,&nanti);
  Amax=ngluons+nquarks+nanti;
  pqlist=new Tpqlist *[Amax+1];
  for(A=0;A<=Amax;A++) pqlist[A]=new Tpqlist(int(1.5*ngluons)+nquarks+nanti);
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
	
  dtot=0;
  nsinglets=0;
  for(pq=pqlist[Amax]->first;pq!=NULL;pq=pq->next){
			//printf("(%d,%d), n=%d\n",pq->p,pq->q,pq->n);
    dtot+=degen(pq->p,pq->q)*pq->n;
    if(pq->p==0 && pq->q==0) nsinglets+=pq->n;
  }
  dtrue=pow(double(3),nquarks+nanti)*pow(8.0,ngluons);
  printf("dtot=%g =? %g\n",dtot,dtrue);
  printf("Nsinglets=%g\n",nsinglets);
		//for(A=0;A<=Amax;A++) pqlist[A]->print();
	pqlist[Amax]->print();
  return 0;
}
