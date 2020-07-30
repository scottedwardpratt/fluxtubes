#include "su3.h"
using namespace NS_SU3;
using namespace std;

int main(){
	int Amax,p,q,A1,A2,pmax;
	double N12=0.0;
	double CZbar,Cbar,w;
	CRandy *randy=new CRandy(-time(NULL));
	char filename[150];
	vector<vector<vector<double>>> pqcount;
  printf("Enter Amax: ");
  scanf("%d",&Amax);
	pmax=GetPmax(Amax);
	CalcPQCount(Amax,pqcount);
	
	sprintf(filename,"Cbar_Amax%d.dat",Amax);
	FILE *fptr=fopen(filename,"w");
	for(A1=0;A1<=Amax;A1++){
		A2=Amax-A1;
		N12=0.0;
		Cbar=CZbar=0.0;
		for(p=0;p<=pmax;p++){
			for(q=0;q<=pmax;q++){
				w=(pqcount[A1][p][q]/degen(p,q))*(pqcount[A2][p][q]/degen(p,q));
				N12+=w;
				CZbar+=w;
				Cbar+=w*Casimir(p,q);
			}
		}
		fprintf(fptr,"%3d %g\n",A1,Cbar/CZbar);
		printf("A1=%3d  Cbar=%g\n",A1,Cbar/CZbar);
	}
	fclose(fptr);
	
	return 0;
}
