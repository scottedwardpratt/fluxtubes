#ifndef __TRAJECTORY_CC__
#define __TRAJECTORY_CC__

#include "su3.h"

using namespace NS_SU3;
using namespace std;

void Ctraj::FindOpenTrajectory(int Nquarks0,int Ngluons,int Nquarksf,CRandy* &randy){
	int dpgluon[7]={2,1,0,-1,-2,1,-1},dqgluon[7]={-1,1,0,2,1,-2,-1};
	int dpquark[3]={1,-1,0},dqquark[3]={0,1,-1};
	//int dantip[3]={0,1,-1},dqanti[3]={1,-1,0};
	int ipq,igluon,iquark,pnew,qnew,count;
	double wsum,rcheck,weight;
	Resize(Nquarks0+Ngluons+Nquarksf);
	ptraj[0]=0;
	qtraj[0]=0;
		
	for(iquark=1;iquark<=Nquarks0;iquark++){
		rcheck=randy->ran();
		wsum=0.0;
		ipq=0;
		do{
			if(ipq==3){
				printf("ipq=%d\n",ipq);
				exit(1);
			}
			pnew=ptraj[iquark-1]+dpquark[ipq];
			qnew=qtraj[iquark-1]+dqquark[ipq];
			if(pnew>=0 && qnew>=0){
				weight=degen(pnew,qnew)/double(degen(ptraj[iquark-1],qtraj[iquark-1])*3);
				wsum+=weight;
				if(wsum<rcheck)
					ipq+=1;
				else{
					ptraj[iquark]=ptraj[iquark-1]+dpquark[ipq];
					qtraj[iquark]=qtraj[iquark-1]+dqquark[ipq];
					//printf("iquark=%d, ptraj=%d, qtraj=%d\n",iquark, ptraj[iquark],qtraj[iquark]);
				}
			}
			else
				ipq+=1;
		}while(wsum<rcheck && wsum<=1.0 && ipq<3);
		if(wsum>1.0){
			printf("wsum=%g, must be <=1\n",wsum);
			exit(1);
		}
	}
	
	for(igluon=Nquarks0+1;igluon<=Nquarks0+Ngluons;igluon++){
		rcheck=randy->ran();
		wsum=0.0;
		ipq=0;
		do{
			if(ipq==7){
				printf("ipq=%d\n",ipq);
				exit(1);
			}
			pnew=ptraj[igluon-1]+dpgluon[ipq];
			qnew=qtraj[igluon-1]+dqgluon[ipq];
			if(pnew>=0 && qnew>=0){
				count=1;
				if(ipq==2 && ptraj[igluon-1]!=0 && qtraj[igluon-1]!=0)
					count=2;
				weight=count*degen(pnew,qnew)/double(degen(ptraj[igluon-1],qtraj[igluon-1])*8);
				wsum+=weight;
				if(wsum<rcheck)
					ipq+=1;
				else{
					ptraj[igluon]=pnew;
					qtraj[igluon]=qnew;
					//printf("igluon=%d, ptraj=%d, qtraj=%d\n",igluon, ptraj[igluon],qtraj[igluon]);
				}
			}
			else
				ipq+=1;
		}while(wsum<rcheck && wsum<=1.0 && ipq<7);
		if(wsum>1.000001 || ipq==7){
			printf("wsum=%g,ipq=%d\n",wsum,ipq);
			exit(1);
		}
	}
	
	wsum=0.0;	
	for(iquark=Nquarks0+Ngluons+1;iquark<=Nquarks0+Ngluons+Nquarksf;iquark++){
		rcheck=randy->ran();
		wsum=0.0;
		ipq=0;
		do{
			if(ipq==3){
				printf("ipq=%d\n",ipq);
				exit(1);
			}
			pnew=ptraj[iquark-1]+dpquark[ipq];
			qnew=qtraj[iquark-1]+dqquark[ipq];
			if(pnew>=0 && qnew>=0){
				weight=degen(pnew,qnew)/double(degen(ptraj[iquark-1],qtraj[iquark-1])*3);
				wsum+=weight;
				if(wsum<rcheck)
					ipq+=1;
				else{
					ptraj[iquark]=ptraj[iquark-1]+dpquark[ipq];
					qtraj[iquark]=qtraj[iquark-1]+dqquark[ipq];
					//printf("iquark=%d, ptraj=%d, qtraj=%d\n",iquark, ptraj[iquark],qtraj[iquark]);
				}
			}
			else
				ipq+=1;
		}while(wsum<rcheck && wsum<=1.0 && ipq<3);
		if(wsum>1.0){
			printf("wsum=%g\n",wsum);
			exit(1);
		}
	}
	
}
	

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
	printf("\n");
}

void Ctraj::Write(FILE *fptr){
	int a;
	for(a=0;a<=A;a++){
		fprintf(fptr,"%4d %4d %4d %9.3f %9.3f\n",a,ptraj[a],qtraj[a],Cquad[a],Ccubic[a]);
	}
}


#endif
