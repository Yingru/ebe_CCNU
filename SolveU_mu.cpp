#include<iostream>
#include<fstream>
#include<cmath>
#include"SolveU_mu.h"
#include"PhyConst.h"
//#include"nr_lusolve.h"

using namespace std;

inline double Det(double A[4][4]){
   return A[1][1]*A[2][2]*A[3][3] - A[1][1]*A[2][3]*A[3][2] \
         -A[1][2]*A[2][1]*A[3][3] + A[1][2]*A[2][3]*A[3][1] \
         +A[1][3]*A[2][1]*A[3][2] - A[1][3]*A[2][2]*A[3][1];
}//Calculate |A_{ij}| i,j from 1 to 3

inline void RestA(double AI[4][4],double A[4][4], double B[4], int I){
         for(int r=0; r<=3; r++)
            for(int c=0; c<=3; c++){
               if(c==I)AI[r][c]=B[r];
               else AI[r][c]=A[r][c];
         }
}


bool LinearSolve(double A[4][4], double B[4], double U[4]){
     double A1[4][4], A2[4][4], A3[4][4];     
     RestA(A1, A, B, 1);      
     RestA(A2, A, B, 2);      
     RestA(A3, A, B, 3);      
     double detA=Det(A);
     assert(detA != 0);
     U[0]=0.0;
     U[1]=Det(A1)/detA;
     U[2]=Det(A2)/detA;
     U[3]=Det(A3)/detA;
     return true;
}

bool SolveUmu3(double Tmn[4][4], double g[4], double Umu[4], double &e){
    const double NCount=1.0E7;
    const double acu=1.0E-5;
    double e_old = 0.0;
    double U_mu[4];   //U_{mu} = g_{\mu\nu}U^{\nu}

    assert(Tmn[0][0]<0 || Tmn[0][0]>0 || Tmn[0][0]==0);
    if(Tmn[0][0]<PhyConst::acu){
	e=0.0;
	Umu[0]=1.0;
	Umu[1]=0.0;
	Umu[2]=0.0;
	Umu[3]=0.0;
    }
    else{
	assert(Tmn[0][0]>PhyConst::acu);
	for(int i=1; i<=3; i++){
	    Umu[i] =Tmn[0][i]/Tmn[0][0];
	    U_mu[i]=g[i]*Umu[i];
	}
	Umu[0] = sqrt(max(1.0-Umu[1]*U_mu[1]-Umu[2]*U_mu[2]-Umu[3]*U_mu[3],1.0E-30));  //normalization
	U_mu[0]= g[0]*Umu[0];
	for(int k=0; k<=3; k++)  e += Tmn[0][k]*U_mu[k]/U_mu[0];

	/*****Solve the equation 
	  (A0*U1+A1*U2)=B1
	  (A2*U1+A3*U2)=B2
	 ************************/
	double A0, A1, A2, A3, B1, B2;
	int iter=0;
	while( fabs(e-e_old) > acu ){
	    iter ++; 
	    //assert(iter<NCount); /*If more than NCount loops, stop program*/
        if(iter>NCount){
        cout<<"spend too much time!"<<endl;
        cout<<"e="<<e<<endl;
        cout<<"e_old="<<e_old<<endl;
        break;
        }
	    e_old = e;
        double A[4][4],B[4],d;
        int n=3;
        int *indx;

        for(int i=0; i<4; i++){
          B[i]=-Tmn[i][0]*g[0]*Umu[0];
          for(int j=0; j<4; j++){
             if(i!=j)A[i][j]=g[j]*Tmn[i][j];
             else A[i][j]= g[j]*Tmn[i][j]-e;
          }
        }

        
        LinearSolve(A,B,Umu);
        for(int i=1; i<=3; i++){
            U_mu[i] = g[i]*Umu[i];
        }


        //ludcmp(A,n,indx,&d);
        //lubksb(A,n,indx,B );   
        //for(int i=1; i<=3; i++){
        //    Umu[i] = B[i];
        //    U_mu[i] = g[i]*Umu[i];
        //}

	    Umu[0] = sqrt(max(1.0-Umu[1]*U_mu[1]-Umu[2]*U_mu[2]-Umu[3]*U_mu[3],1.0E-30));  //normalization
        U_mu[0]= g[0]*Umu[0];

	    e = 0.0;
	    for(int k=0; k<=3; k++)  e += Tmn[0][k]*U_mu[k]/U_mu[0];
	}
    cout<<"iter="<<iter<<"\n";
    }
    return true;
}




bool SolveUmu2(double Tmn[4][4], double g[4], double Umu[4], double &e){
    const double NCount=1.0E7;
    const double acu=1.0E-5;
    double e_old = 0.0;
    double U_mu[4];   //U_{mu} = g_{\mu\nu}U^{\nu}
    Umu[3]=0.0;
    U_mu[3]=0.0;

    assert(Tmn[0][0]<0 || Tmn[0][0]>0 || Tmn[0][0]==0);
    if(Tmn[0][0]<PhyConst::acu){
	e=0.0;
	Umu[0]=1.0;
	Umu[1]=0.0;
	Umu[2]=0.0;
	Umu[3]=0.0;
    }
    else{
	assert(Tmn[0][0]>PhyConst::acu);
	for(int i=1; i<=2; i++){
	    Umu[i] =Tmn[0][i]/Tmn[0][0];
	    U_mu[i]=g[i]*Umu[i];
	}
	Umu[0] = sqrt(max(1.0-Umu[1]*U_mu[1]-Umu[2]*U_mu[2]-Umu[3]*U_mu[3],1.0E-30));  //normalization
	U_mu[0]= g[0]*Umu[0];
	for(int k=0; k<=2; k++)  e += Tmn[0][k]*U_mu[k]/U_mu[0];

	/*****Solve the equation 
	  (A0*U1+A1*U2)=B1
	  (A2*U1+A3*U2)=B2
	 ************************/
	double A0, A1, A2, A3, B1, B2;
	int iter=0;
	while( fabs(e-e_old) > acu ){
	    iter ++; 
	    assert(iter<NCount); /*If more than NCount loops, stop program*/
	    e_old = e;

	    A0 = g[1]*Tmn[1][1]-e;
	    A1 = g[2]*Tmn[1][2];
	    A2 = g[1]*Tmn[2][1];
	    A3 = g[2]*Tmn[2][2]-e;
	    B1 =-g[0]*Tmn[0][1]*Umu[0];
	    B2 =-g[0]*Tmn[0][2]*Umu[0];

	    Umu[2] = (B2-B1*A2/A0)/(A3-A1*A2/A0);
	    Umu[1] = (B1-A1*Umu[2])/A0;
	    Umu[0] = sqrt(max(1.0-Umu[1]*U_mu[1]-Umu[2]*U_mu[2]-Umu[3]*U_mu[3],1.0E-30));  //normalization
	    for(int i=0; i<=2; i++)U_mu[i] = g[i]*Umu[i];

	    e = 0.0;
	    for(int k=0; k<=2; k++)  e += Tmn[0][k]*U_mu[k]/U_mu[0];
	}
    cout<<"iter="<<iter<<"\n";
    }
    return true;
}



bool SolveUmu1(double Tmn[4][4], double g[4], double Umu[4], double &e){
    const double NCount=1.0E6;
    const double acu=1.0E-7;
    double e_old = -1.0;
    double V[4];
    double U_mu[4];   //U_{mu} = g_{\mu\nu}U^{\nu}
    Umu[3]=0.0;
    U_mu[3]=0.0;

    assert(Tmn[0][0]>0.0);//T^{\tau\tau} should be larger than 0;
    if(Tmn[0][0]<=0.0){
	e=0.0;
	Umu[0]=1.0;
	Umu[1]=0.0;
	Umu[2]=0.0;
	Umu[3]=0.0;
    }
    else{

	e=Tmn[0][0];
	for(int i=1; i<=2; i++){
	    Umu[i] =Tmn[0][i]/Tmn[0][0];
	    U_mu[i]=g[i]*Umu[i];
	}
	Umu[0] = sqrt(1.0-Umu[1]*U_mu[1]-Umu[2]*U_mu[2]-Umu[3]*U_mu[3]);  //normalization
	U_mu[0]= g[0]*Umu[0];

	int iter=0; 
	while( fabs(e-e_old) > acu ){
	    iter ++; 
	    //assert(iter<NCount); /*If more than NCount loops, stop program*/
	    if(iter>NCount)return false;
	    e_old = e;

	    e = 0;
	    for(int k=0; k<=2; k++)  e += Tmn[0][k]*U_mu[k]/U_mu[0];
	    for(int i=1; i<=2; i++) {
		Umu[i] = 0;
		for(int k=0; k<3; k++) Umu[i] +=  Tmn[i][k]/e * U_mu[k];
		U_mu[i] = g[i]*Umu[i];
	    }
	    //	Umu[1] = (Tmn[1][0]*U_mu[0] + Tmn[1][2]*U_mu[2])/e/(1.0-Tmn[1][1]/e*g[1]);
	    //
	    //
	    //	Umu[2] = (Tmn[2][0]*U_mu[0] + Tmn[2][1]*U_mu[1])/e/(1.0-Tmn[2][2]/e*g[2]);
	    //
	    //
	    //	U_mu[1]=Umu[1]*g[1];
	    //	U_mu[2]=Umu[2]*g[2];

	    Umu[0] = sqrt(1.0-Umu[1]*U_mu[1]-Umu[2]*U_mu[2]-Umu[3]*U_mu[3]);  //normalization
	    U_mu[0]= Umu[0]*g[0];
	}
    }
    //cout<<"iter="<<iter<<endl;
    return true;
}

bool SolveUmu(double Tmn[4][4], double g[4], double Umu[4], double &e)
{
    const double NCount=1.0E6;
    const double acu=1.0E-6;
    double e_old = -1.0;
    double U_mu[4];   //U_{mu} = g_{\mu\nu}U^{\nu}

    assert(Tmn[0][0]>0.0);//T^{\tau\tau} should be larger than 0;

    e=Tmn[0][0];
    Umu[3]=0.0;
    U_mu[3]=0.0;
    for(int i=1; i<3; i++){
	Umu[i] =Tmn[0][i]/Tmn[0][0];
	U_mu[i]=g[i]*Umu[i];
    }

    /* The comment below are another way of initialization */	
    Umu[0] = sqrt(1.0-Umu[1]*U_mu[1]-Umu[2]*U_mu[2]-Umu[3]*U_mu[3]);  //normalization
    U_mu[0]= g[0]*Umu[0];

    //	e=0;
    //	for(int i=0; i<4; i++)
    //		for(int j=0; j<4; j++){
    //			e += U_mu[i]*Tmn[i][j]*U_mu[j];
    //		}
    //	/*Calc energy density e from e=u_{\mu}T^{\mu\nu}u_{\nu} */

    int iter=0; 
    while( fabs(e-e_old) > acu ){
	iter ++; assert(iter<NCount); /*If more than NCount loops, stop program*/
	e_old = e;

	for(int i=0; i<3; i++){
	    double A=0.0;                       //A=T^{\rho\mu}U_{\rho}
	    for(int j=0; j<3; j++){
		A += Tmn[i][j] * U_mu[j];
	    }
	    Umu[i] = A / e_old;
	    U_mu[i]= g[i]*Umu[i];
	}//Calc U^{\mu} by iteration
	Umu[0] = sqrt(1.0-Umu[1]*U_mu[1]-Umu[2]*U_mu[2]-Umu[3]*U_mu[3]);  //normalization
	U_mu[0]= g[0]*Umu[0];

	e=0;
	for(int i=0; i<3; i++)
	    for(int j=0; j<3; j++){
		e += U_mu[i]*Tmn[i][j]*U_mu[j];
	    }
	/*Calc energy density e from e=u_{\mu}T^{\mu\nu}u_{\nu} */
    }
    //cout<<"Iterat "<<iter<<" times"<<endl;
    return true;
}


//#define TEST
#ifdef TEST
int main(int argc, char **argv)
{
    double tau0=0.4;  
    double m=0.14;
    double vx=0.5;
    double vy=0.4;
    double vz=0.3;

    double P[4]={0.0 , m*vx, m*vy, m*vz};
    double X[4]={0.0};
    int NTOTAL;
    //ifstream fin_ampt("ampt_Ini/P0.dat");
    ifstream fin_ampt("../ampt-v1.21-v2.21/ana/ampt_b0-10_Ini/P0.dat");
    //ifstream fin_ampt("results/ED_XY0.dat");
    fin_ampt>>NTOTAL;
    //while(!fin_ampt.eof()){
    for(int i=0; i<4; i++)fin_ampt>>P[i];
    for(int i=0; i<4; i++)fin_ampt>>X[i];

    double MT = sqrt(P[0]*P[0] - P[3]*P[3]);
    double Y  = 0.5*log((P[0]+P[3])/(P[0]-P[3]));
    double etas=0.5*log((X[0]+X[3])/(X[0]-X[3]));
    //tau0=sqrt(X[0]*X[0]-X[3]*X[3]);
    double v_etas=tanh(Y-etas)/tau0;
    cout<<"tau0="<<tau0<<" Setting v_etas="<<v_etas<<endl;
    cout<<"mass="<<sqrt(MT*MT-P[1]*P[1]-P[2]*P[2])<<endl;
    P[0]= MT*cosh(Y-etas);
    P[3]= MT*sinh(Y-etas)/tau0;

    
    double T[4][4];
    double U[4];
    /*
    P[0] = sqrt( m*m+P[1]*P[1]+P[2]*P[2]+P[3]*P[3] );
    double MT = sqrt(P[0]*P[0] - P[3]*P[3]);
    double Y  = 0.5*log((P[0]+P[3])/(P[0]-P[3]));
    double etas=1.0;
    double v_etas=tanh(Y-etas)/tau0;
    cout<<"Setting v_etas="<<v_etas<<endl;
    P[0]= MT*cosh(Y-etas);
    P[3]= MT*sinh(Y-etas)/tau0;
    */


    double g[4]={1.0,-1.0, -1.0, -(tau0*tau0)};
    for(int i=0; i<4; i++){
	for(int j=0; j<4; j++){
	T[i][j] = P[i]*P[j]/P[0];  // /(tau0*dx*dy*detas);
	}
    }
    cout<<"The solution is:\n";
    cout<<"U^{mu}= ";
    double e;
	SolveUmu3(T,g,U,e);
	for(int i=0; i<4; i++)
	    cout<<U[i]/U[0]<<"  ";
	cout<<endl;
	cout<<"energy density ="<<e<<endl;
    cout<<endl;
    //}
}
#endif
