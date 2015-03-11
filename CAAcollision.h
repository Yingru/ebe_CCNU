#ifndef CAACOLLISION_H
#define CAACOLLISION_H

#include<cmath>

using namespace std;

typedef const double CD;
class CNucleus{
	public:
		double A;          // mass of atom
		double Ro0;        // mass density
		double R;          //nuclear radii (unit: fm)
		double Eta;        //diffussion parameter of nuclear radii(unit: fm)
		CNucleus(void);	   
		CNucleus(const CNucleus& rhs);
		CNucleus(CD& A1, CD& ro0, CD& r, CD& eta);
		double Thickness(CD& x, CD& y);//Thickness at (x, y)
};

class CAAcollision
{
	public:
		CNucleus Ap;        //Projectile nuclear
		CNucleus At;        //Target nuclear
		double Si0;         //Cross section 
		double b;           //Impact parameter
		double kNw;         //Contribution of wounded collisions dafault 0.75
		double kNb;         //Contribution of binary collisions  dafault 0.25
		double Eta_flat;    //rapidity flat
		double Eta_gw;      //rapidity gaussion width
		double Nw(CD& x, CD& y);   //Number of wounded collisions
		double Nb(CD& x, CD& y);   //Number of binary collisions
		double NCollision(CD& x, CD& y, CD& z);
		//constructor 
		CAAcollision(void);
		CAAcollision(const CNucleus& A1, const CNucleus &A2, CD& sig,CD& b0, CD& fNw, CD& fNb, CD& Eflat, CD& Egw);
};//sig-Si0  b0-b  fNw-kNw  fNb-kNb



#endif
