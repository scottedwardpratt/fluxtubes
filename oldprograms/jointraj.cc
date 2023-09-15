#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <cmath>
#include "subs/randy.h"

using namespace std;

#include "subs/misc.cc"
#include "subs/pq.cc"
#include "subs/tableaux.cc"
#include "subs/addmults.cc"
#include "subs/randy.cc"

class Ctraj{
public:
	vector<int> ptraj,qtraj,casimir;
	int A;
	void Resize(int Aset){
		A=Aset;
		ptraj.resize(A+1);
		qtraj.resize(A+1);
		casimir.resize(A+1);
	}
	void CalcCasimirs();
	void PrintCasimirs();
	Ctraj(int Aset){
		Resize(Aset);
	}
	void MakeTrajectory(vector<vector<vector<int>>> &Npq,CRandy* &randy,int A);
	void Print();
};

void Ctraj::MakeTrajectory(vector<vector<vector<int>>> &Npq,CRandy* &randy,int A){
	Resize(A);
	double d,weight,CZ;
	int Amax,a,a1,a2,pmax1,pmax2,p,q,pp,qq;
	int n1,n2,n,flip;
	FILE *fptr,*fptr1,*fptr2;
	char filename1[100],filename2[100];
	char dummy[400];
	
	Amax=Npq.size();
	ptraj.resize(2*Amax+1);
	qtraj.resize(2*Amax+1);
	if(A>2*Amax){
		printf("A=%d is too big, >2Amax, Amax=%d\n",A,Amax);
		exit(1);
	}
	
	if(A%2==0){
		a1=a2=A/2;
	}
	else{
		a1=(A+1)/2;
		a2=A-a1;
	}
	CZ=0.0;
	pmax1=lrint(1.5*a1);
	for(p=0;p<pmax1;p++){
		for(q=0;q<p;q++){
			d=degen(p,q);
			weight=0.5*Npq[a1][p][q]*Npq[a2][p][q]/(d*d);
			if(p==q)
				weight*=2.0;
			CZ+=weight;
		}
	}
	do{
		p=floorl(pmax1*randy->ran());
		q=floorl(pmax1*randy->ran());
		d=degen(p,q);
		weight=0.5*Npq[a1][p][q]*Npq[a2][p][q]/(d*d);
		if(p==q)
			weight*=2.0;
	}while(weight/CZ<randy->ran());  // use better method
	//
	pp=p; qq=q;
	if(p<q){
		pp=q; qq=p;
	}
	sprintf(filename1,"opentrajdata/A%d/p%d_q%d.dat",a1,pp,qq);
	fptr1=fopen(filename1,"r");
	if(a1!=a2){
		sprintf(filename2,"opentrajdata/A%d/p%d_q%d.dat",a2,pp,qq);
		fptr2=fopen(filename2,"r");
	}
	else
		fptr2=fptr1;
	
	n1=floorl(Npq[a1][pp][qq]*randy->ran());
	n2=floorl(Npq[a2][pp][qq]*randy->ran());
	

	
	if(a1==a2){
		if(n1>n2){
			flip=n1; n1=n2; n2=flip;
		}
	}
	for(n=0;n<n1;n++)
		fgets(dummy,400,fptr1);
	for(a=0;a<=a1;a++){
		fscanf(fptr1,"%d %d",&ptraj[a],&qtraj[a]);
	}
	if(a1==a2){
		for(n=n1;n<n2;n++)
			fgets(dummy,400,fptr2);
	}
	else{
		for(n=0;n<n2;n++)
			fgets(dummy,400,fptr2);
	}
	for(a=0;a<a2;a++){
		fscanf(fptr2,"%d %d",&ptraj[A-a],&qtraj[A-a]);
	}
	pp=ptraj[a1]; qq=qtraj[a1];
	fscanf(fptr2,"%d %d",&ptraj[a1],&qtraj[a1]);
	fclose(fptr1);
	if(a1!=a2)
		fclose(fptr2);
	
	//Flip if p and q are flipped
	if(pp!=ptraj[a1] || qq!=qtraj[a1]){
		if(pp==qtraj[a1] && qq==ptraj[a1]){
			for(a=a1+1;a<=A;a++){
				flip=ptraj[a]; ptraj[a]=qtraj[a]; qtraj[a]=flip;
			}
		}
		else{
			printf("DISASTER, pp=%d, qq=%d, ptraj[a1]=%d, qtraj[a1]=%d\n",pp,qq,ptraj[a1],qtraj[a1]);
			exit(1);
		}
	}
	CalcCasimirs();
}

void Ctraj::PrintCasimirs(){
	for(int a=0;a<=A;a++){
		printf("%6d ",casimir[a]);
	}
	printf("\n");
}

void Ctraj::CalcCasimirs(){
	for(int a=0;a<=A;a++){
		casimir[a]=Casimir(ptraj[a],qtraj[a]);
	}
}

void Ctraj::Print(){
	int a;
	for(a=0;a<=A;a++){
		printf("(%3d,%3d) ",ptraj[a],qtraj[a]);
		if(a==A/2)
			printf("|\n|");
	}
	printf("\n");
}

void MakeNpq(vector<vector<vector<int>>> &Npq,int Amax){
	int pmax,a,p,q;
	char filename[150];
	FILE *fptr;
	pmax=int(Amax*1.5);
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
			sprintf(filename,"opentrajdata/Npq_A%d.dat",a);
			fptr=fopen(filename,"r");
			fscanf(fptr,"%d %d",&p,&q);
			while(!feof(fptr)){
				fscanf(fptr,"%d",&Npq[a][p][q]);
				fscanf(fptr,"%d %d",&p,&q);
			}
			fclose(fptr);
		}
	}
}



int main(){
	vector<vector<vector<int>>> Npq;
	int A,Ntraj;
	printf("Enter A and Ntraj: ");
	scanf("%d %d",&A,&Ntraj);
	MakeNpq(Npq,40);
	CRandy *randy=new CRandy(-time(NULL));
	Ctraj traj(0);
	for(int itraj=0;itraj<Ntraj;itraj++){
		traj.MakeTrajectory(Npq,randy,A);
		//traj.PrintCasimirs();		
	}
}
