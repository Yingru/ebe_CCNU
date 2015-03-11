#include"ran2.h"
#include"PhyConst.h"
#include"HijingIni.h"
#include<cmath>
#include<iostream>

using namespace std;

double fws(double r){
    double ws=1.0/(exp((r-6.37)/0.54)+1.0);
    return ws*ws;
}//Give the Hiing initial particles a squared woods-saxon distribution

void wsxy(double &x, double &y){
/*Generate a woods-saxon distribution x and y */
    double r, phi, fr;
    double fws_max=fws(0.0);
    
    do{
        r=15.0*ran2(&NSEED);        // x
        fr=ran2(&NSEED);            // f(x)
      }while(fr>fws(r)/fws_max);

    phi=ran2(&NSEED)*2.0*PhyConst::pi;

    x=r*cos(phi);
    y=r*sin(phi);
}

void generate_xy(double &x, double &y, CD &E, CD &PX, CD &PY, CD &PZ){
     double pT2 = PX*PX + PY*PY;
     double pTjet=3.5*3.5;  //jet energy cut 1.5GeV/c
     double r,phi;
     if(pT2<pTjet)wsxy(x,y);
     else{
       r=6.0+ran2(&NSEED);  //jet_r= 6-7 fm
       phi=ran2(&NSEED)*2.0*PhyConst::pi;
       x=r*cos(phi);
       y=r*sin(phi);
     }
}//generate x,y for hijing initial condition


//int main(int argc, char**argv)
//{
//    double x,y;
//    generate_xy(x,y,1.6,1.1,1.1,0.2);
//    cout<<"x="<<x<<" y="<<y<<endl;
//    return 0;
//}
