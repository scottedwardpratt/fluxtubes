#ifndef __TRAJECTORY_CC__
#define __TRAJECTORY_CC__

#include "su3.h"

using namespace NS_SU3;
using namespace std;

void Ctraj::FindTrajectory(int p0,int q0,int A,vector<vector<vector<double>>> pqcount,CRandy* &randy){
	Resize(A);
	int pmax=GetPmax(A+2*p0+2*q0);
	Tpqlist *goodlist=new Tpqlist(pmax+2);
	Tpq *pq;
	int p,q,pprime,qprime,a,ipq,count;
	int dp[7]={2,1,0,-1,-2,1,-1},dq[7]={-1,1,0,2,1,-2,-1};
	double r,wsum,wtot,w[7];
	p=p0;
	q=q0;
	ptraj[0]=p0;
	qtraj[0]=q0;
	for(a=1;a<=A;a++){
		wtot=0.0;
		//goodlist=new Tpqlist(pmax+2);
		goodlist->clear();
		//for(pprime=ppmin;pprime<=ppmax;pprime++){
		//for(qprime=qpmin;qprime<=qpmax;qprime++){
		for(ipq=0;ipq<7;ipq++){
			pprime=p-dp[ipq];
			qprime=q-dq[ipq];
			if(pprime>=0 && qprime>=0 && pprime<=pmax && qprime<=pmax){
				if(pprime!=0 || qprime!=0 || p!=0 || q!=0){
					count=1;
					if(ipq==2)
						count=2;
					goodlist->add(pprime,qprime,count);
				}
			}
		}
		wtot=0.0;
		ipq=0;
		for(pq=goodlist->first;pq!=NULL;pq=pq->next){
			pprime=pq->p; qprime=pq->q;
			w[ipq]=pqcount[A-a][pprime][qprime]*pq->n/degen(pprime,qprime);
			wtot+=w[ipq];
			ipq+=1;
		}
		//printf("a=%d, wtot=%g\n",a,wtot);
		r=randy->ran();
		pq=goodlist->first;
		wsum=0.0;
		ipq=0;
		do{
			pprime=pq->p; qprime=pq->q;
			wsum+=w[ipq]/wtot;
			if(wsum>1.000001){
				printf("wsum=%g\n",wsum);
				exit(1);
			}
			pq=pq->next;
			ipq+=1;
		}while(r>wsum);
		ptraj[a]=pprime; qtraj[a]=qprime;
		p=pprime; q=qprime;
		goodlist->clear();
		//delete goodlist;

	}
	goodlist->remove();
	delete goodlist;
	//Print();
}

void Ctraj::PrintCasimirs(){
	for(int a=0;a<=A;a++){
		printf("%9.3f  %9.3f",Cquad[a],Ccubic[a]);
	}
	printf("\n");
}

void Ctraj::CalcCasimirs(){
	for(int a=0;a<=A;a++){
		NS_SU3::CalcCasimirs(ptraj[a],qtraj[a],Cquad[a],Ccubic[a]);
	}
}

void Ctraj::Print(){
	int a;
	for(a=0;a<=A;a++){
		printf("(%2d,%2d)",ptraj[a],qtraj[a]);
	}
	printf("\n -----------------------------------------------------------\n");
}

void Ctraj::Write(FILE *fptr){
	int a;
	for(a=0;a<=A;a++){
		fprintf(fptr,"%4d %4d %4d %9.3f %9.3f\n",a,ptraj[a],qtraj[a],Cquad[a],Ccubic[a]);
	}
}


#endif
