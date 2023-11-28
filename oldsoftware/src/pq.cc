#ifndef __PQ_CC__
#define __PQ_CC__
#include "su3.h"
using namespace NS_SU3;
using namespace std;

Tpq::Tpq(int pset,int qset,double nset){
	p=pset;
	q=qset;
	n=nset;
	next=NULL;
}

Tpqlist::Tpqlist(int pqmax_set){
	int p,q;
	pqmax=pqmax_set;
	first=NULL;
	last=NULL;
	pqptr_array=new Tpq **[pqmax+1];
	for(p=0;p<=pqmax;p++){
		pqptr_array[p]=new Tpq *[pqmax+1];
		for(q=0;q<=pqmax;q++)
			pqptr_array[p][q]=NULL;
	}
}

void Tpqlist::print(){
	Tpq *pqptr;
	pqptr=first;
	do{
		printf("(%d,%d), degen=%d, number=%g\n",pqptr->p,pqptr->q,
		degen(pqptr->p,pqptr->q),pqptr->n);
		pqptr=pqptr->next;
	} while(pqptr!=NULL);
}

void Tpqlist::compress(){
	Tpq *pq1,*pq2,*oldpq1,*pqtemp,*nextpq1;

	for(pq2=first->next;pq2!=NULL;pq2=pq2->next){
		oldpq1=NULL;

		for(pq1=first;pq1!=pq2;pq1=nextpq1){

			nextpq1=pq1->next;
			if(pq1->p==pq2->p && pq1->q==pq2->q){
				pqptr_array[pq2->p][pq2->q]=pq2;
				if(pq1==first){
					pqtemp=first;
					first=first->next;
					pq1=first;
				}
				else{
					pqtemp=pq1;
					if(oldpq1==NULL){
						printf("oldpq1=NULL?? error!\n");
						exit(1);
					}
					oldpq1->next=pq1->next;
					pq1=oldpq1;
				}
				pq2->n+=pqtemp->n;
				delete pqtemp;
			}
			else{
				oldpq1=pq1;
			}
		}
	}

}

void Tpqlist::add(int p,int q,double n){
	Tpq *pq;
	if(p>pqmax || q>pqmax){
		printf("Tpqlist::add -- Trying to add to list something that is out of bounds!\n");
		printf("p=%d, q=%d, pqmax=%d\n",p,q,pqmax);
		exit(1);
	}
	if(pqptr_array[p][q]==NULL){
		pq=new Tpq(p,q,n);
		if(first==NULL){
			first=pq;
			last=pq;
		}
		else{
			last->next=pq;
			last=pq;
		}
		pqptr_array[p][q]=pq;
	}
	else{
		pq=pqptr_array[p][q];
		pq->n+=n;
	}
	if(p>pqmax||q>pqmax){
		printf("Tpqlist::add  p,q too large=%d,%d\n",p,q);
		exit(1);
	}
}

void Tpqlist::clear(){
	Tpq *next;
	next=first;
	while(first!=NULL){
		next=first->next;
		delete first;
		first=next;
	}
	first=NULL;
	last=NULL;
	int p,q;
	for(p=0;p<=pqmax;p++){
		for(q=0;q<=pqmax;q++){
			pqptr_array[p][q]=NULL;
		}
	}
}

void Tpqlist::remove(){
	Tpq *next;
	next=first;
	while(first!=NULL){
		next=first->next;
		delete first;
		first=next;
	}
	first=NULL;
	last=NULL;
	int p,q;
	for(p=0;p<=pqmax;p++){
		for(q=0;q<=pqmax;q++){
			pqptr_array[p][q]=NULL;
		}
		delete pqptr_array[p];
	}
	delete pqptr_array;
}

void Tpqlist::cgluon(int ell){
	clear();
	if(ell==1) add(1,1,1.0);
	if(ell==2){
		add(2,2,1.0);
		add(0,3,-1.0);
		add(3,0,-1.0);
		add(0,0,1.0);
	}
	if(ell>=3){
		add(ell,ell,1.0);
		add(ell-2,ell+1,-1.0);
		add(ell+1,ell-2,-1.0);
		add(ell-2,ell-2,-1.0);
		add(ell,ell-3,1.0);
		add(ell-3,ell,1.0);
		add(0,0,2.0);
	}
}

void Tpqlist::cquark(int ell){
	clear();
	if(ell==1) add(1,0,1.0);
	if(ell==2){
		add(2,0,1.0);
		add(0,1,-1.0);
	}
	if(ell>=3){
		add(ell,0,1.0);
		add(ell-2,1,-1.0);
		add(ell-3,0,1.0);
	}
}

#endif
