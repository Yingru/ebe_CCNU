#ifndef SOLVEU_MU_H
#define SOLVEU_MU_H
#include<iostream>
#include<cassert>
#include<cmath>
using namespace std;

bool SolveUmu(double Tmn[4][4], double g[4], double Umu[4], double &e);

/*Solve 4 velocity U^{mu} and energy density e from EM tensor T^{\mu\nu} */
/*Tmn = T^{\mu\nu},  g=diag{g_{00},g_{11},g_{22},g_{33}}, Umu=U^{\mu}, */
bool SolveUmu1(double Tmn[4][4], double g[4], double Umu[4], double &e);
bool SolveUmu2(double Tmn[4][4], double g[4], double Umu[4], double &e);

bool SolveUmu3(double Tmn[4][4], double g[4], double Umu[4], double &e);


class CT44{
    public:
	int i, j, k;
	double T[4][4];
	CT44():i(0),j(0),k(0){
	    for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
		    T[i][j] = 0;
	}
	CT44 & operator+=(const CT44 &rhs){
	    for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
		    this->T[i][j] += rhs.T[i][j];
	    return *this;
	}
	bool InSameCell(const CT44 &rhs){
	    return ( (i==rhs.i)&&(j==rhs.j)&&(k==rhs.k) );
	}
};



#endif
