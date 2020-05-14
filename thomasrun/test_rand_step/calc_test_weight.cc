#include "su3.h"

using namespace NS_SU3;
using namespace std;

#include <stdio.h>

//point of this program is to
//calculate the probability of going to (p', q', n') from
// (p,q,n) at a time a

void
get_good_list (int p_given, int q_given, int p_new, int q_new, int delta_p[7],
               int delta_q[7], int maxval, int num_multiplet,
               Tpqlist * good_list);

int
main (int argc, char **argv)
{
  int p1,q1,n1,pprime,qprime,a,Amax,pmax,ipq; 
  int dp[7] = { 2, 1, 0, -1, -2, 1, -1 }, dq[7] = { -1, 1, 0, 2, 1, -2, -1 };
  vector<vector<vector<double>>> pqcount;
  double w[7];
  double wtot;

  a = atoi(argv[1]);
  p1 = atoi(argv[2]);			//grab a,p,q,n and the Amax used.
  q1 = atoi(argv[3]);
  n1 = atoi(argv[4]);
  Amax = atoi(argv[5]);
  CalcPQCount(Amax, pqcount); //object for calculating weights
  pmax = Amax;


  Tpq *pq;
  Tpqlist *goodlist= new Tpqlist(pmax+1);
  
  get_good_list(p1,q1,pprime,qprime,dp,dq,pmax,n1,goodlist);	  
  
  ipq = 0;  
  wtot = 0.;

  for(pq=goodlist->first;pq!=NULL;pq=pq->next){
      pprime=pq->p; qprime=pq->q;
      w[ipq]=pqcount[Amax-a][pprime][qprime]*pq->n/degen(pprime,qprime);
      wtot+=w[ipq];
      ipq+=1;
    }
  
  char filename[100];
  sprintf(filename, "calculated_weights/weights_%d_%d_%d_%d_%d_.txt",a,p1,q1,n1,Amax);

  FILE *fptr = fopen (filename, "w");
  fprintf (fptr, "p'    q'    w\n"); 

  ipq = 0;
  for(pq=goodlist->first;pq!=NULL;pq=pq->next)
  { 
      pprime=pq->p; qprime=pq->q;
      double weight = w[ipq]/wtot;
      char output_string[150];
      sprintf(output_string,"%i    %i    %f\n",pprime,qprime,weight);
      fprintf(fptr,"%s",output_string);
      ipq += 1;
  }

  fclose(fptr);
  delete pq;
  goodlist -> clear();
 }


void
get_good_list (int p_given, int q_given, int p_new, int q_new, int delta_p[7],
              int delta_q[7], int maxval, int num_multiplet,
              Tpqlist * good_list)
{
  for (int iter_pq = 0; iter_pq < 7; iter_pq++)
    {                           //using the rules of young tableaux diagrams, get the possible multiplets we may jump too

      if ((p_given == 0) && (q_given == 0))
        {
          good_list->add (1, 1, 1);
          break;
        }

      p_new = p_given - delta_p[iter_pq];
      q_new = q_given - delta_q[iter_pq];

      if (p_new >= 0 && q_new >= 0 && p_new <= maxval && q_new <= maxval)
        {
          num_multiplet = 1;
          if (iter_pq == 2)
            num_multiplet = 2;
          good_list->add (p_new, q_new, num_multiplet);
        }

    }
}


