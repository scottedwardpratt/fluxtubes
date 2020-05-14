#include "su3.h"
using namespace NS_SU3;
using namespace std;

#include <stdio.h>

//Purpose of this program is to investigate if 
//a truly random walk through p-q space and see if the 
//same weights appear as specified in the more guided 
//random walk done in traj.cc

void
get_good_list (int p_given, int q_given, int p_new, int q_new, int delta_p[7],
	       int delta_q[7], int maxval, int num_multiplet,
	       Tpqlist * good_list);

int
main ()
{

  int Amax, p1, q1, n1, pmax, pprime, qprime, a, ipq, count, nprime;
  double r, wsum, wtot, w[7];
  CRandy *randy = new CRandy (-time (NULL));
  char filename[150];

  Amax = 9;			//Seting Amax to a relatively low number so that we can generate a lot of trajectories quickly
  pmax = GetPmax (Amax);

  long Ntraj, itraj;
  Ntraj = 5e9;

  sprintf (filename, "rand_traj/good_random_traj_%d.txt", Amax);	//open a file for our trajectories that come back to zero

  FILE *fptr = fopen (filename, "w");
  fprintf (fptr, "t    p    q    n\n");

  int dp[7] = { 2, 1, 0, -1, -2, 1, -1 }, dq[7] = { -1, 1, 0, 2, 1, -2, -1 };	//initalize tableaux "machinery"
  Tpqlist *goodlist = new Tpqlist (pmax);
  Tpq *pq;

  bool good_traj;		// boolean for if our trajectory is good or not
  vector < int >traj_vectors[Amax + 1];	//storing our trajectory with vectors
  vector < int >curr_traj {0, 0, 0, 0};

  for (itraj = 0; itraj < Ntraj; itraj++)	//generate Ntraj total random walks
    {

      good_traj = false;
      p1 = q1 = curr_traj[0] = curr_traj[1] = curr_traj[2] = curr_traj[3] = 0;

      traj_vectors[0] = curr_traj;

      for (a = 1; a <= Amax; a++)	// adding 'A' total gluons to our trajectory
	{
	  wtot = 0.0;
	  goodlist->clear ();
	  get_good_list (p1, q1, pprime, qprime, dp, dq, pmax, count, goodlist);	//grabs the allowed states we can jump too and stores them in goodlist

	  wtot = 0.0;		//we're gonna weight which multiplet to jump to by the number of states in each multiplet
	  ipq = 0;
	  for (pq = goodlist->first; pq != NULL; pq = pq->next)	//calculate the weight
	    {
	      pprime = pq->p;
	      qprime = pq->q;
	      w[ipq] = degen (pprime, qprime) * pq->n;
	      wtot += w[ipq];
	      ipq += 1;
	    }

	  r = randy->ran ();	//grab a number on the interval [0,1]
	  pq = goodlist->first;
	  wsum = 0.0;
	  ipq = 0;

	  do			//add the weights, and the first one to pass 'r' is the multiplet we select
	    {
	      pprime = pq->p;
	      qprime = pq->q;
	      nprime = pq->n;
	      wsum += w[ipq] / wtot;
	      if (wsum > 1.000001)
		{
		  printf ("wsum=%g\n", wsum);
		  exit (1);
		}
	      pq = pq->next;
	      ipq += 1;
	    }

	  while (r > wsum);

	  curr_traj[0] = a;
	  p1 = curr_traj[1] = pprime;
	  q1 = curr_traj[2] = qprime;
	  n1 = curr_traj[3] = nprime;

	  if ((a == Amax) && (p1 == 0) & (q1 == 0))	//check if we landed back on the origin.
	    {
	      good_traj = true;
	    }

	  traj_vectors[a] = curr_traj;	//add this value of (p,q) to our array
	  goodlist->clear ();
	}

      int len = sizeof (traj_vectors) / sizeof (traj_vectors[0]);

      if (good_traj == false)
	{
	  //if it's a 'bad' trajectory, start back from the origin.
	  good_traj = false;
	  continue;
	}

      else
	{
	  //if it's a good trajectory, we want to save it.

	  int p_iter, q_iter, n_iter;
	  char output_string[100];

	  for (int i = 0; i < len; i++)	//output the step, and p and q to a file
	    {

	      p_iter = traj_vectors[i][1];
	      q_iter = traj_vectors[i][2];
	      n_iter = traj_vectors[i][3];

	      sprintf (output_string, "%i    %i    %i    %i\n", i, p_iter,
		       q_iter, n_iter);

	      fprintf (fptr, "%s", output_string);

	    }

	  good_traj = false;	//reset good_traj and make another trajectory
	  continue;

	}

    }

  fclose (fptr);
  return 0;
}


void
get_good_list (int p_given, int q_given, int p_new, int q_new, int delta_p[7],
	      int delta_q[7], int maxval, int num_multiplet,
	      Tpqlist * good_list)
{
  for (int iter_pq = 0; iter_pq < 7; iter_pq++)
    {				//using the rules of young tableaux diagrams, get the possible multiplets we may jump too

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
