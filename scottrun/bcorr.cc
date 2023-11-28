#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	long long int itraj,Ntraj;
	int A,Amax,Btarget,Bproj,Ngluons,Nquarks0,Nquarksf;
	CRandy *randy=new CRandy(-time(NULL));
	vector<int> plist;
	vector<int> qlist;
	vector<double> balance;
	Ctraj *traj;
	int Nsinglets=0,p,q,oldp,oldq,b,NTrajWrite=100,ia,ib,dely;
	long long int Nb=0;
	double q2,q3;
	vector<double> Q2bar,Q3bar,pminusqbar;
 
	//printf("Enter B of target: ");
	//scanf("%d",&Btarget);
	//printf("Enter B of projectile: ");
	//scanf("%d",&Bproj);
	Btarget=Bproj=0;
	printf("Enter Ngluons: ");
	scanf("%d",&Ngluons);
	printf("Enter Ntraj: ");
	scanf("%lld",&Ntraj);
	
	Nquarks0=3*Btarget;
	Nquarksf=3*Bproj;
	Amax=Nquarks0+Ngluons+Nquarksf;
	Q2bar.resize(Amax+1);
	Q3bar.resize(Amax+1);
	pminusqbar.resize(Amax+1);
	balance.resize(Amax+1);
	for(A=0;A<=Amax;A++){
		Q2bar[A]=Q3bar[A]=pminusqbar[A]=balance[A]=0.0;
	}
	
	traj=new Ctraj(Nquarks0+Ngluons+Nquarksf);
	string filename,dirname="trajectories/B0"+to_string(Btarget)+"Bf"+to_string(Bproj)+"Ng"+to_string(Ngluons);
	string command="mkdir -p "+dirname;
	system(command.c_str());
	command="rm -r -f "+dirname+"/*.txt";
	system(command.c_str());
	FILE *fptr;
	fptr=fopen(filename.c_str(),"w");
	
	int N3=0;
	double pq3bar=0.0;
	
	for(itraj=0;itraj<Ntraj;itraj++){
		traj->FindOpenTrajectory(Nquarks0,Ngluons,Nquarksf,randy);
		if(traj->ptraj[Amax]==0 && traj->qtraj[Amax]==0){
			plist.clear();
			qlist.clear();
			//traj->Print();
			N3+=1;
			pq3bar+=traj->ptraj[3]-traj->qtraj[3];
			if(Nsinglets<NTrajWrite){
				filename=dirname+"/traj"+to_string(Nsinglets)+".txt";
				//printf("filename=%s\n",filename.c_str());
				fptr=fopen(filename.c_str(),"w");
			}
			Nsinglets+=1;
			for(A=0;A<Amax;A++){
				p=traj->ptraj[A];
				q=traj->qtraj[A];
				pminusqbar[A]+=p-q;
				CalcCasimirs(p,q,q2,q3);
				Q2bar[A]+=q2;
				Q3bar[A]+=q3;
				if(Nsinglets<=NTrajWrite){
					if(A==Nquarks0-1)
						fprintf(fptr,"0 0 0 0 0\n");
					if(A>=Nquarks0 && A<Amax+1-Nquarksf)
						fprintf(fptr,"%5d %5d %5d %10g %10g\n",A+1-Nquarks0,p,q,q2,q3);
					if(A==Amax-Nquarksf)
						fprintf(fptr,"%5d 0 0 0 0\n",Ngluons+2);
				}
				if(A>0){
					b=(p-q-oldp+oldq)/3;
					if(b>0){
						plist.push_back(A);
						Nb+=1;
					}
					if(b<0){
						qlist.push_back(A);
						Nb+=1;
					}
					if(abs(b)>1){
						printf("b=%d\n",b);
						Nb+=1;
					}
				}
				oldp=p;
				oldq=q;
			}
			for(ia=0;ia<plist.size();ia++){
				for(ib=ia;ib<plist.size();ib++){
					if(ia!=ib){
						dely=abs(lrint(plist[ia]-plist[ib]));
						balance[dely]-=1.0;
						dely=abs(lrint(qlist[ia]-qlist[ib]));
						balance[dely]-=1.0;
					}
					dely=abs(lrint(plist[ia]-qlist[ib]));
					balance[dely]+=1.0;
					dely=abs(lrint(qlist[ia]-plist[ib]));
					balance[dely]+=1.0;	
				}
			}
			if(!feof(fptr))
				fclose(fptr);
			if(plist.size()!=qlist.size()){
				printf("plist/qlist sizes do not match!!!1 = %lu,%lu\n",plist.size(),qlist.size());
				exit(1);
			}
			
		}
		if(10*(itraj+1)%Ntraj==0)
			printf("finished %g percent\n",100*double(itraj+1)/double(Ntraj));
	}
	
	printf("pq3bar=%g\n",pq3bar/double(N3));
	printf("Nsinglets=%d\n",Nsinglets);
	
	filename="results/B0"+to_string(Btarget)+"Bf"+to_string(Bproj)+"Ng"+to_string(Ngluons)+".txt";
	fptr=fopen(filename.c_str(),"w");
	FILE *bfile=fopen("bcorr.txt","w");
	double bnorm=0.0;
	for(A=0;A<=Amax;A++){
		balance[A]=balance[A]/Nb;
		bnorm+=balance[A];
		fprintf(bfile,"%d  %g\n",A,balance[A]);
		pminusqbar[A]=pminusqbar[A]/double(Nsinglets);
		Q2bar[A]=Q2bar[A]/double(Nsinglets);
		Q3bar[A]=Q3bar[A]/double(Nsinglets);
		//printf("%5d %9.4f %9.4f %9.4f\n",A,Q2bar[A],Q3bar[A],pminusqbar[A]);
		fprintf(fptr,"%5d %9.4f %9.4f %9.4f\n",A,Q2bar[A],Q3bar[A],pminusqbar[A]);
	}
	fclose(fptr);
	fclose(bfile);
	printf("bnorm=%g\n",bnorm);
	
	return 0;
}
