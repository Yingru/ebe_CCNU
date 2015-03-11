#ifndef __CEOSH
#define __CEOSH
#include<cmath>
#include<iostream>
#include<deque>
#include<algorithm>
#include"PhyConst.h"
#include<fstream>

using namespace std;
using namespace PhyConst;
typedef const double CD;

/*Tabular Eos s95p-v1 */
struct EPTable{
  double e0, de;
  int Ne;
  deque<double> e;
  deque<double> p;
  deque<double> s;
  deque<double> T;
};

class CEOSL{
 public:
  EPTable t[4];   //The 4 tables for different energy region
  CEOSL();
  void init();
  double lin_int(CD &x1, CD &x2, CD &y1, CD &y2, CD &x);
  double P(CD &eps);
  double T(CD &eps);
  double S(CD &eps);
};



class CEos{
	public:
	int IEOS;
	CEOSL eosl;

	CEos();
	CEos(int ieos);
	void SetEos(int ieos);
	double P(CD &eps, CD &bd );//pressure density
	double T(CD &eps, CD &bd );//temperature
	double S(CD &eps, CD &bd );//entropy density
};


#endif
