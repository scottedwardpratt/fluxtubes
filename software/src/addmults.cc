#ifndef __SU3_MISC_CC__
#define __SU3_MISC_CC__

#include "su3.h"
using namespace std;
using namespace NS_SU3;

//#define PRINT_TABLEAUS

void NS_SU3::RanStepWeighted(int &p,int &q,double randy,double &weight){
	int p1,q1,p2,q2;
	p1=p; q1=q;
	p2=q2=1; // gluons
	int degen1=degen(p1,q1);
	double ransum=0.0,maxsum=0.0;
	 
	int ip2,iq2,i,na0min,na1min,na1max,nb1min,nb1max,iw;
	int na[3],nb[3],nx[3];
	Ttableauxlist alist;
	Ttableaux *atableaux;
	Ttableaux *abtableaux;
	Ttableaux *testtableaux;

	for(i=0;i<3;i++) nb[i]=0;
	nx[0]=p1+q1;
	nx[1]=q1;
	nx[2]=0;

	na0min=p2+q2-(p1+q1);
	if(na0min<0) na0min=0;
	for(na[0]=na0min;na[0]<=p2+q2;na[0]++){
		na1min=(p2+q2-na[0])-q1;
		if(na1min<0) na1min=0;
		na1max=p1;
		if(na1max>p2+q2-na[0]) na1max=p2+q2-na[0];
		for(na[1]=na1min;na[1]<=na1max;na[1]++){
			na[2]=(p2+q2-na[0]-na[1]);
			atableaux=new Ttableaux(na,nb,nx);
			alist.add(atableaux);
		}
	}

	for(atableaux=alist.first;atableaux!=NULL;
	atableaux=atableaux->next){
		nb1min=q2-atableaux->length[1]+atableaux->length[2];
		if(nb1min<0) nb1min=0;
		nb1max=q2;
		if(nb1max>atableaux->length[0]-atableaux->length[1])
			nb1max=atableaux->length[0]-atableaux->length[1];
		if(nb1max>atableaux->na[0]) nb1max=atableaux->na[0];
		for(nb[1]=nb1min;nb[1]<=nb1max;nb[1]++){
			nb[2]=q2-nb[1];
			if(nb[2]<=atableaux->na[0]+atableaux->na[1]-nb[1]){
				p=atableaux->na[0]+nx[0]-atableaux->na[1]-nb[1]-nx[1];
				q=atableaux->na[1]+nb[1]+nx[1]-atableaux->na[2]-nb[2];
				maxsum+=1.0;
			}
		}
	}
	weight=maxsum/8.0;
	//printf("maxsum=%g, weight=%g\n",maxsum,weight);

	randy*=maxsum;
	for(atableaux=alist.first;atableaux!=NULL;
	atableaux=atableaux->next){
		nb1min=q2-atableaux->length[1]+atableaux->length[2];
		if(nb1min<0) nb1min=0;
		nb1max=q2;
		if(nb1max>atableaux->length[0]-atableaux->length[1])
			nb1max=atableaux->length[0]-atableaux->length[1];
		if(nb1max>atableaux->na[0]) nb1max=atableaux->na[0];
		for(nb[1]=nb1min;nb[1]<=nb1max;nb[1]++){
			nb[2]=q2-nb[1];
			if(nb[2]<=atableaux->na[0]+atableaux->na[1]-nb[1]){
				p=atableaux->na[0]+nx[0]-atableaux->na[1]-nb[1]-nx[1];
				q=atableaux->na[1]+nb[1]+nx[1]-atableaux->na[2]-nb[2];
				ransum+=1.0;
				if(ransum>maxsum+0.5){
					printf("ransum=%g, maxsum=%g\n",ransum,maxsum);
				}
				if(ransum>randy){
					alist.clear();
					return;
				}
			}
		}
	}
	alist.clear();
	printf("DIDN'T FIND RANSTEP!!!!\n");
	exit(1);
}

void NS_SU3::RanStep(int &p,int &q,double randy,double &weight){
	weight=1.0;
	int p1,q1,p2,q2;
	p1=p; q1=q;
	p2=q2=1; // gluons
	int degen1=degen(p1,q1);
	double ransum=0.0;
	 
	int ip2,iq2,i,na0min,na1min,na1max,nb1min,nb1max,iw;
	int na[3],nb[3],nx[3];
	Ttableauxlist alist;
	Ttableaux *atableaux;
	Ttableaux *abtableaux;
	Ttableaux *testtableaux;
	for(i=0;i<3;i++) nb[i]=0;
	nx[0]=p1+q1;
	nx[1]=q1;
	nx[2]=0;

	na0min=p2+q2-(p1+q1);
	if(na0min<0) na0min=0;
	for(na[0]=na0min;na[0]<=p2+q2;na[0]++){
		na1min=(p2+q2-na[0])-q1;
		if(na1min<0) na1min=0;
		na1max=p1;
		if(na1max>p2+q2-na[0]) na1max=p2+q2-na[0];
		for(na[1]=na1min;na[1]<=na1max;na[1]++){
			na[2]=(p2+q2-na[0]-na[1]);
			atableaux=new Ttableaux(na,nb,nx);
			alist.add(atableaux);
		}
	}

	for(atableaux=alist.first;atableaux!=NULL;
	atableaux=atableaux->next){
		nb1min=q2-atableaux->length[1]+atableaux->length[2];
		if(nb1min<0) nb1min=0;
		nb1max=q2;
		if(nb1max>atableaux->length[0]-atableaux->length[1])
			nb1max=atableaux->length[0]-atableaux->length[1];
		if(nb1max>atableaux->na[0]) nb1max=atableaux->na[0];
		for(nb[1]=nb1min;nb[1]<=nb1max;nb[1]++){
			nb[2]=q2-nb[1];
			if(nb[2]<=atableaux->na[0]+atableaux->na[1]-nb[1]){
				p=atableaux->na[0]+nx[0]-atableaux->na[1]-nb[1]-nx[1];
				q=atableaux->na[1]+nb[1]+nx[1]-atableaux->na[2]-nb[2];
				ransum+=degen(p,q)/(8.0*degen1);
				if(ransum>randy){
					alist.clear();
					return;
				}
			}
		}

	}
	alist.clear();
	printf("DIDN'T FIND RANSTEP!!!!\n");
	exit(1);
}

void NS_SU3::addmults(int p1,int q1,int p2,int q2,Tweightlist *weightlist){

	if(fabs(weightlist->inweight[0])>1.0E-10){
		int p,q,ip2,iq2,i,na0min,na1min,na1max,nb1min,nb1max,iw;
		int na[3],nb[3],nx[3];
		Ttableauxlist *alist;
		Ttableaux *atableaux;
		Ttableaux *abtableaux;
		Ttableaux *testtableaux;

		alist=new Ttableauxlist;

		for(i=0;i<3;i++) nb[i]=0;
		nx[0]=p1+q1;
		nx[1]=q1;
		nx[2]=0;

		na0min=p2+q2-(p1+q1);
		if(na0min<0) na0min=0;
		for(na[0]=na0min;na[0]<=p2+q2;na[0]++){
			na1min=(p2+q2-na[0])-q1;
			if(na1min<0) na1min=0;
			na1max=p1;
			if(na1max>p2+q2-na[0]) na1max=p2+q2-na[0];
			for(na[1]=na1min;na[1]<=na1max;na[1]++){
				na[2]=(p2+q2-na[0]-na[1]);
				atableaux=new Ttableaux(na,nb,nx);
				alist->add(atableaux);
			}
		}

		for(atableaux=alist->first;atableaux!=NULL;
		atableaux=atableaux->next){
			nb1min=q2-atableaux->length[1]+atableaux->length[2];
			if(nb1min<0) nb1min=0;
			nb1max=q2;
			if(nb1max>atableaux->length[0]-atableaux->length[1])
				nb1max=atableaux->length[0]-atableaux->length[1];
			if(nb1max>atableaux->na[0]) nb1max=atableaux->na[0];
			for(nb[1]=nb1min;nb[1]<=nb1max;nb[1]++){
				nb[2]=q2-nb[1];
				if(nb[2]<=atableaux->na[0]+atableaux->na[1]-nb[1]){
					p=atableaux->na[0]+nx[0]-atableaux->na[1]-nb[1]-nx[1];
					q=atableaux->na[1]+nb[1]+nx[1]-atableaux->na[2]-nb[2];
					for(iw=0;iw<weightlist->length;iw++){
						weightlist->outweight[iw][p][q]+=weightlist->inweight[iw];
					}
				}
			}

		}
		alist->clear();
	}
}

void NS_SU3::addgluon_list(int p,int q,int n0,Tpqlist *pqlist){
	int dp[7]={2,1,0,-1,-2,1,-1},dq[7]={-1,1,0,2,1,-2,-1};
	int pprime,qprime,ipq,count,ncheck=0;
	if(p==0 && q==0){
		pqlist->add(1,1,n0);
		ncheck+=8;
	}
	else{
		for(ipq=0;ipq<7;ipq++){
			pprime=p+dp[ipq]; qprime=q+dq[ipq];
			if(pprime>pqlist->pqmax || qprime>pqlist->pqmax){
				printf("NS_SU3::addgluon_list -- pprime, qprime too big, =%d,%d\n",pprime,qprime);
				exit(1);
			}
			if(pprime>=0 && qprime>=0){
				count=1;
				if(ipq==2 && p!=0 && q!=0)
					count=2;
				pqlist->add(pprime,qprime,n0*count);
				ncheck+=count*degen(pprime,qprime);
			}
		}
	}
	if(ncheck!=8*degen(p,q)){
		printf("ncheck=%d =? %d\n",ncheck,8*degen(p,q));
		exit(1);
	}
}

void NS_SU3::addquark_list(int p,int q,int n0,Tpqlist *pqlist){
	int dp[3]={1,-1,0},dq[3]={0,1,-1};
	int pprime,qprime,ipq,count,ncheck=0;
	if(p==0 && q==0){
		pqlist->add(1,0,n0);
		ncheck=3;
	}
	else{
		for(ipq=0;ipq<3;ipq++){
			pprime=p+dp[ipq]; qprime=q+dq[ipq];
			if(pprime>pqlist->pqmax || qprime>pqlist->pqmax){
				printf("NS_SU3::addgluon_list -- pprime, qprime too big, =%d,%d\n",pprime,qprime);
				exit(1);
			}
			if(pprime>=0 && qprime>=0){
				count=1;
				//if(ipq==2 && pprime!=0 && qprime!=0)
					//count=2;
				pqlist->add(pprime,qprime,n0*count);
				ncheck+=count*degen(pprime,qprime);
			}
		}
	}
	if(ncheck!=3*degen(p,q)){
		printf("ncheck=%d =? %d\n",ncheck,3*degen(p,q));
		printf("pq=(%d,%d), p'q'=(%d,%d)\n",p,q,pprime,qprime);
		exit(1);
	}
}

void NS_SU3::addantiquark_list(int p,int q,int n0,Tpqlist *pqlist){
	int dp[3]={0,1,-1},dq[3]={1,-1,0};
	int pprime,qprime,ipq,count,ncheck=0;
	if(p==0 && q==0){
		pqlist->add(1,0,n0);
		ncheck=3;
	}
	else{
		for(ipq=0;ipq<3;ipq++){
			pprime=p+dp[ipq]; qprime=q+dq[ipq];
			if(pprime>pqlist->pqmax || qprime>pqlist->pqmax){
				printf("NS_SU3::addgluon_list -- pprime, qprime too big, =%d,%d\n",pprime,qprime);
				exit(1);
			}
			if(pprime>=0 && qprime>=0){
				count=1;
				//if(ipq==2 && pprime!=0 && qprime!=0)
					//count=2;
				pqlist->add(pprime,qprime,n0*count);
				ncheck+=count*degen(pprime,qprime);
			}
		}
	}
	if(ncheck!=3*degen(p,q)){
		printf("ncheck=%d =? %d\n",ncheck,3*degen(p,q));
		printf("pq=(%d,%d), p'q'=(%d,%d)\n",p,q,pprime,qprime);
		exit(1);
	}
}

void NS_SU3::addmults_list(int p1,int q1,double n1,
int p2,int q2,long long n2,Tpqlist *pqlist){
	int p,q,ip2,iq2,i,na0min,na1min,na1max,nb1min,nb1max;
	int delp,delpq;
	int na[3],nb[3],nx[3];
	Ttableauxlist *alist;
	Ttableauxlist *ablist;
	Ttableaux *atableaux;
	Ttableaux *abtableaux;
	alist=new Ttableauxlist;
	for(i=0;i<3;i++)
		nb[i]=0;
	nx[0]=p1+q1;
	nx[1]=q1;
	nx[2]=0;

	na0min=p2+q2-(p1+q1);
	if(na0min<0) na0min=0;
	for(na[0]=na0min;na[0]<=p2+q2;na[0]++){
		na1min=(p2+q2-na[0])-q1;
		if(na1min<0) na1min=0;
		na1max=p1;
		if(na1max>p2+q2-na[0]) na1max=p2+q2-na[0];
		for(na[1]=na1min;na[1]<=na1max;na[1]++){
			na[2]=(p2+q2-na[0]-na[1]);
			atableaux=new Ttableaux(na,nb,nx);
			alist->add(atableaux);
		}
	}

	ablist=new Ttableauxlist;
	for(atableaux=alist->first;atableaux!=NULL;
	atableaux=atableaux->next){
#ifdef PRINT_TABLEAUS    
		printf("___________________________________________\n");
#endif
		nb1min=q2-atableaux->length[1]+atableaux->length[2];
		if(nb1min<0) nb1min=0;
		nb1max=q2;
		if(nb1max>atableaux->length[0]-atableaux->length[1])
			nb1max=atableaux->length[0]-atableaux->length[1];
		if(nb1max>atableaux->na[0]) nb1max=atableaux->na[0];
		for(nb[1]=nb1min;nb[1]<=nb1max;nb[1]++){
			nb[2]=q2-nb[1];
			if(nb[2]<=atableaux->na[0]+atableaux->na[1]-nb[1]){
				abtableaux=new  Ttableaux(atableaux->na,nb,nx);
				ablist->add(abtableaux);
#ifdef PRINT_TABLEAUS    
				abtableaux->print();
#endif
				q=abtableaux->length[1]-abtableaux->length[2];
				p=abtableaux->length[0]-abtableaux->length[1];

				if(p>pqlist->pqmax || q>pqlist->pqmax){
					printf("addmults_list, p,q too big, %d,%d, pqmax=%d\n",p,q,pqlist->pqmax);
				}
				pqlist->add(p,q,n1*n2);
				delp=p-p1;
				delpq=p+q-p1-q1;
				if(abs(delpq)>2)
					printf("delp=%d, p1,q1=%d,%d, p2,q2=%d,%d, N=%g, p,q=%d,%d\n",delp,p1,q1,p2,q2,n1*n2,p,q);
			}
		}
	}
	alist->clear();
	ablist->clear();
	delete alist;
	delete ablist;
}

#endif
