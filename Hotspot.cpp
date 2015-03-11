    #include<iostream>
    #include<cmath>
    #include"ran2.h"
    long NSEED=5111111;
    double pi=3.1415926;
    using namespace std;
    typedef const double CD;

    double Ed_hotspot(CD &x, CD &y, CD &h){
        CD emax=5.0;// Gev/fm^3
        CD sigm=0.1;//Gaussian width
        CD coef=emax/sqrt(2.0*pi)/sigm;
        CD r=7.0; // radius position of the hot spot
        double phi=ran2(&NSEED);
        cout<<"phi="<<phi<<endl;
        double result;

        if(abs(h)<3.0){
           double dx= x - r*cos(phi);
           double dy= y - r*sin(phi);
           double dr2= dx*dx + dy*dy;
           result = coef*exp(-dr2/(2.0*0.1*0.1));
        }
        else result=0.0;
        return result;
    }

    int main(){
      cout<<Ed_hotspot(7,0,0)<<endl;
      return 0;
    }
