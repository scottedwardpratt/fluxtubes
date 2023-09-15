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
	int pa,qa,pb,qb,P,Q;
	double C1a,C2a,C1b,C2b;
	double C1bar,C2bar,C1f,C2f;
	Tpq *pq;
  Tpqlist *pqlist;
	
  printf("Enter pa,qa: ");
	scanf("%d %d",&pa,&qa);
  //printf("Enter pb,qb: ");
	//scanf("%d %d",&pb,&qb);
	pb=0; qb=1;
	
	Casimir12(pa,qa,C1a,C2a);
	printf("For (pa=%d, qa=%d) C1=%lf, C2=%lf\n",pa,qa,C1a,C2a);
	Casimir12(pb,qb,C1b,C2b);
	printf("For (pb=%d, qb=%d) C1=%lf, C2=%lf\n",pb,qb,C1b,C2b);

  pqlist=new Tpqlist(pa+pb+qa+qb+5);
	pqlist->clear();
	addmults_list(pa,qa,1,pb,qb,1,pqlist);
	pqlist->compress();

	C1bar=C2bar=0.0;
  dtot=0;
  nsinglets=0;
  for(pq=pqlist->first;pq!=NULL;pq=pq->next){
		Casimir12(pq->p,pq->q,C1f,C2f);
    dtot+=degen(pq->p,pq->q)*pq->n;
		C1bar+=degen(pq->p,pq->q)*pq->n*C1f;
		C2bar+=degen(pq->p,pq->q)*pq->n*C2f;
    if(pq->p==0 && pq->q==0)
			nsinglets+=pq->n;
  }
  dtrue=degen(pa,qa)*degen(pb,qb);
	C1bar=C1bar/double(dtot);
	C2bar=C2bar/double(dtot);
	
  printf("Total states = %g =? %g\n",dtot,dtrue);
  printf("Nsinglets = %g, <C1>=%g, Delta<C1>=%g,  <C2>=%g, Delta C2=%g\n",nsinglets,C1bar,C1bar-C1a,C2bar,C2bar-C2a);
	pqlist->print();
  return 0;
}
