#include "CAAcollision.h"



//constructor
//overloading construct function
CAAcollision::CAAcollision(const CNucleus& A1, const CNucleus & A2, CD& sig, CD& b0, CD& fNw, CD& fNb, CD& Eflat, CD& Egw)\
		:Ap(A1),At(A2)
{
	Si0=sig;   //cross section
	b=b0;    //impact parameter
	kNw=fNw;
	kNb=fNb;
	Eta_flat=Eflat;
	Eta_gw=Egw;
}

double CAAcollision::Nw(CD& x, CD& y)
{
	double Ta=Ap.Thickness(x-0.5*b,y);
	double Tb=At.Thickness(x+0.5*b,y);
	return Ta*(1.0-pow(1.0-Si0*Tb/Ap.A, Ap.A))\
		+Tb*(1.0-pow(1.0-Si0*Ta/At.A, At.A));
}//Wounded collisions

double CAAcollision::Nb(CD& x, CD& y)
{
	double Ta=Ap.Thickness(x-0.5*b,y);
	double Tb=At.Thickness(x+0.5*b,y);
	return Si0*Ta*Tb;
}//Binary collisions

double CAAcollision::NCollision(CD& x,CD& y, CD& z)
{
	double heta;

    if(abs(z)>Eta_flat) heta=exp(-pow(abs(z)-Eta_flat,2.0)/(2.0*Eta_gw*Eta_gw));
    //if(abs(z)>Eta_flat) heta=2.0/(exp((abs(z)-Eta_flat)/Eta_gw)+1.0);
	else heta=1.0;

	return (kNw*Nw(x,y) + kNb*Nb(x,y))*heta;
}//Total collisions


CNucleus::CNucleus(const CNucleus& rhs)\
		:A(rhs.A), Ro0(rhs.Ro0), R(rhs.R), Eta(rhs.Eta)
{
}// Copy construct function

CNucleus::CNucleus(CD& A1, CD& ro0, CD& r, CD& eta)\
		:A(A1),Ro0(ro0),R(r),Eta(eta)
{
}//Nucleus class construct function

//Calc Nucleus thickness for Glauber model
double CNucleus::Thickness(CD& x, CD& y)
{
//	double r, RhoA;
//	double thickness=0.0;
//	double ddz = R/100.0;
//	double z = -5.0*R;
//	while(z<=5.0*R){
//		r=sqrt(x*x+y*y+z*z);
//		RhoA=Ro0/(exp((r-R)/Eta)+1.0);
//		thickness = thickness + RhoA*ddz;
//		z=z+ddz;
//	}//This method is slow. The two methods generate the same result
//	return thickness;

	double r,thickness,f,z ,cut,a[10],b[10];
	double zk[5]={-0.9061798, -0.5386492, 0.0, 0.5386493, 0.9061798};//range(-1,1)
	double Ak[5]={0.2369269,0.4786287,0.5688889,0.4786287,0.2369269};//weight
	thickness = 0.0;
	cut  = 5.0*R;
	for (int j=0; j!=10; ++j){
		a[j] = j*cut/10.0;
		b[j] = (j+1)*cut/10.0;
		for(int i=0; i!=5; ++i){
			z=(a[j]+b[j])/2.0 + (b[j]-a[j])/2.0* zk[i];
			r=sqrt( x*x+y*y+z*z );
			f=Ro0/(1.0+exp((r-R)/Eta));
			thickness= thickness +Ak[i]*f*(b[j]-a[j])/2.0;
		}
	}
	return 2.0*thickness;
}
/*There is a gauss intergration .
 * Change integration variable z to zk, z=(a+b)/2.0+(b-a)/2.0*zk;
 * Then \int_{a}^{b} dz 
 * =(b-a)/2.0 *\int_{-1}^{1} dzk 
 * =\sum (b-a)/2.0 *A_{i}*g(zk)*/


