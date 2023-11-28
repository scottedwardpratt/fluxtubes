#ifndef __BALANCE_DEFS_H__
#define __BALANCE_DEFS_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstring>
#include <complex>
#include <sys/stat.h>
#include <ctime>
#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <cmath>
#include <vector>
#include <map>
#include <unordered_map>
#include <Eigen/Dense>
#include "constants.h"

using namespace std;

class CHydroBalance;
class CHydroMesh;
class CHyperMesh;
class CEoS;
class CCharge;
class CPart;
class CResList;
class CResInfo;
class CBranchInfo;
class CHyperElement;
class CBalance;
class CAcceptance;
class CAction;
class CB3D;
class CB3DCell;
class CRandom;
class CParameterMap;
class CBalanceArrays;
class CSampler;
class CLocalInfo;
class CAction;
class CRegenerate;
class CSEInfo;
class CHYDROtoB3D;

typedef unordered_map<long int,CResInfo *> CResInfoMap;
typedef pair<long int, CResInfo*> CResInfoPair;
typedef vector<CBranchInfo *> CBranchList; //gives branchlist name
typedef multimap<int,CCharge* > CChargeMap;
typedef pair<int,CCharge* > CChargePair;
typedef multimap<int,CPart* > CPartMap;
typedef pair<int,CPart* > CPartPair;
//typedef array<double,4> FourVector;
typedef double FourVector[4];

typedef multimap<double,CAction *> CActionMap;
typedef pair<double,CAction*> CActionPair;
#endif
