#ifndef __SU3_H__
#define __SU3_H__

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include "constants.h"
#include "defs.h"
#include "randy.h"
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

class Tpq{
public:
	int p,q;
	double n;
	Tpq *next;
	Tpq(int pset,int qset,double nset);
};

class Tpqlist{
public:
	int pqmax;
	Tpq ***pqptr_array;
	Tpq *first;
	Tpq *last;
	Tpqlist(int Amax);
	void add(int p,int q,double n);
	void clear();
	void remove();
	void cgluon(int ell);
	void cquark(int ell);
	void compress();
	void print();
};

class Ttableaux{
public:
  Ttableaux *next;
  int length[3],na[3],nb[3],nx[3];
  Ttableaux(int *naset,int *nbset,int *nxset);
  void copy(Ttableaux *tableaux_ptr);
  void print();
  bool testequal(Ttableaux *other);
  void getpq(int &p,int &q);
};

class Tweightlist{
public:
	int length;
	double inweight[5];
	double **outweight[5];
};

class Ttableauxlist{
public:
  Ttableaux *first;
  Ttableaux *last;
  Ttableauxlist();
  void add(Ttableaux *tableaux);
  void clear();
};

class Ctraj{
public:
	vector<int> ptraj,qtraj;
	vector<double> Cquad,Ccubic;
	int A;
	void Resize(int Aset){
		A=Aset;
		ptraj.resize(A+1);
		qtraj.resize(A+1);
		Cquad.resize(A+1);
		Ccubic.resize(A+1);
	}
	void CalcCasimirs();
	void PrintCasimirs();
	Ctraj(int Aset){
		Resize(Aset);
	}
	void JoinTrajectories(vector<vector<vector<int>>> &Npq,CRandy* &randy,int A,vector<vector<double>> &biggestweight);
	void FindTrajectory(int p0,int q0,int A,vector<vector<vector<double>>> pqcount,CRandy* &randy);
	void Print();
	void Write(FILE *fptr);
};

namespace NS_SU3{
	void RanStepWeighted(int &p,int &q,double randy,double &weight);
	void RanStep(int &p,int &q,double randy,double &weight);
	void addmults(int p1,int q1,int p2,int q2,Tweightlist *weightlist);
	void addmults_list(int p1,int q1,double n1,int p2,int q2,int long long n2,Tpqlist *pqlist);
	void addgluon_list(int p,int q,int n0,Tpqlist *pqlist);
	void addquark_list(int p,int q,int n0,Tpqlist *pqlist);
	void addantiquark_list(int p,int q,int n0,Tpqlist *pqlist);
	void ReadNpq(vector<vector<vector<double>>> &Npq,int Amax,vector<vector<double>> &biggestweight);
	void WriteOpenTrajectories(int Amax,long long int ntraj,CRandy *randy);
	void CalcPQCount(int p0,int q0,int Amax,vector<vector<vector<double>>> &pqcount);
	void ClearPQCount(int Amax,vector<vector<vector<double>>> &pqcount);
	int GetPmax(int A);
	int degen(int p,int q);
	void CalcCasimirs(int p,int q,double &Cquad,double &Ccubic);
	double omega_massless(double T,double V);
	double omegae_massless(double T,double V);
	double omega_bessel(double mass,double T,double V);
	double omegae_bessel(double mass,double T,double V);
	double omegap_bessel(double mass,double T,double V);
}

#endif
