#ifndef __TABLEAUX_CC__
#define __TABLEAUX_CC__

//#define PRINT_TABLEAUS  
#define EPSILON 1.0E-60

#include "su3.h"
using namespace NS_SU3;
using namespace std;

void Ttableaux::getpq(int &p,int &q){
  p=length[0]-length[1];
  q=length[1]-length[2];
}

bool Ttableaux::testequal(Ttableaux *other){
  bool test;
  int i;
  test=1;
  for(i=0;i<3;i++){
    if(na[i]!=other->na[i] || nb[i]!=other->nb[i] || nx[i]!=other->nb[i]){
      test=0;
      goto TestFailure;
    }
  }

 TestFailure:
  return test;
}

void Ttableaux::copy(Ttableaux *tableaux_ptr){
  int i,j;
  //next=tableaux_ptr->next;                                                    
  next=NULL;
  for(i=0;i<3;i++){
    length[i]=tableaux_ptr->length[i];
    na[i]=tableaux_ptr->na[i];
    nb[i]=tableaux_ptr->nb[i];
    nx[i]=tableaux_ptr->nx[i];
  }
}

Ttableaux::Ttableaux(int *naset,int *nbset,int *nxset){
  int i;
  next=NULL;
  for(i=0;i<3;i++){
    na[i]=naset[i];
    nb[i]=nbset[i];
    nx[i]=nxset[i];
    length[i]=na[i]+nb[i]+nx[i];
  }
}

void Ttableaux::print(){
  int i,j,p,q;
  q=length[1]-length[2];
  p=length[0]-length[1];
  printf("___ P=%d ___ Q=%d ___ degen=%d\n",p,q,degen(p,q));
  for(i=0;i<3;i++){
    for(j=0;j<nx[i];j++) printf("x ");
    for(j=0;j<na[i];j++) printf("a ");
    for(j=0;j<nb[i];j++) printf("b ");
    printf("\n");
  }
}

Ttableauxlist::Ttableauxlist(){
  first=NULL;
  last=NULL;
  void clear();
}

void Ttableauxlist::clear(){
  Ttableaux *tptr;
  while(first!=NULL){
    tptr=first->next;
    delete first;
    first=tptr;
  }
  first=NULL;
  last=NULL;
}

void Ttableauxlist::add(Ttableaux *tableaux){
  if(last==NULL){
    first=tableaux;
    last=tableaux;
  }
  else{
    last->next=tableaux;
    last=tableaux;
  }
}

#endif
