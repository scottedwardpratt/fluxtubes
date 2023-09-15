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
  double dtrue,dtot,nsinglets;
	int p1,q1,p2,q2,P,Q;
	double C1,C2;
	Tpq *pq;
  Tpqlist *pqlist;
	
  printf("Enter p1,q1: ");
	scanf("%d %d",&p1,&q1);
  printf("Enter p2,q2: ");
	scanf("%d %d",&p2,&q2);
	
	Casimir12(p1,q1,C1,C2);
	printf("For (p1=%d, q1=%d) C1=%lf, C2=%lf\n",p1,q1,C1,C2);
	Casimir12(p2,q2,C1,C2);
	printf("For (p2=%d, q2=%d) C1=%lf, C2=%lf\n",p2,q2,C1,C2);

  pqlist=new Tpqlist(p1+p2+q1+q2+5);
	pqlist->clear();
	addmults_list(p1,q1,1,p2,q2,1,pqlist);
	pqlist->compress();

	
	double C2bar=0.0,C1f,C2f;
  dtot=0;
  nsinglets=0;
  for(pq=pqlist->first;pq!=NULL;pq=pq->next){
		Casimir12(pq->p,pq->q,C1f,C2f);
    dtot+=degen(pq->p,pq->q)*pq->n;
		C2bar+=degen(pq->p,pq->q)*pq->n*C2f;
    if(pq->p==0 && pq->q==0)
			nsinglets+=pq->n;
  }
  dtrue=degen(p1,q1)*degen(p2,q2);
	C2bar=C2bar/double(dtot);
	
  printf("Total states = %g =? %g\n",dtot,dtrue);
  printf("Nsinglets = %g, C2bar=%g\n",nsinglets,C2bar);
	pqlist->print();
  return 0;
}
