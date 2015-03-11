#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include"mpi.h"
using namespace std;
int myid, size;
MPI_Status status;
const int Num0 = 100;
const int Num1 = 100;
const int Num2 = 100;
const int VELOCITY=1;
const int ANGULAR=2;

const double etamax = 12.0;
const double x1max  = 12.0;
const double x2max  = 12.0;
const double deta = etamax/Num0;
const double dx1  = x1max/Num1;
const double dx2  = x2max/Num2;
double U1[4][100][100][100];
double U2[4][100][100][100];
double mystack[2][100][100][100];
double eta[Num0],x1[Num1],x2[Num2];
double tau,a0, cut, RR, c0;
//((((((((((((((((((((((((((((((((((Function Declear((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
void checkType(int &steptype);//check steptype to see wheather we should use (lnT,y,v1,v2) or (lnT,y,alpha,theta) as the step variable.

double derivative(const int& N1,const int& N2,const int& N3,const int& I,const int& J,const int& K,const int& steptype);//calc d_{\tau},d_{\eta},d_{v1},d_{v2//} or d_{\alpha}, d_{\theta}.

void matrixes(const int &i, const int &j, const int &k, const double &cs2, const double &v1_or_alpha,const double &v2_or_theta, const double &eta,const double &y, double AD0[4][4], double AD1[4][4], double AD2[4][4], const int &steptype);
//construct matrix AD0, AD1, AD2, or M^{-}.AD0.M, M^{-}.AD1.M, M^{-}.AD2.M

void popFromStack(double mystack[2][Num0][Num1][Num2],const int& steptype);//temperaly change alpha theta to v1,v2 in order to write to file.
void pushToStack(double mystack[2][Num0][Num1][Num2],const int& steptype);//change back to alpha theta in U1(2,i,j,k) and U1(3,i,j,k) for ANGULAR step.

void step(const double &dtau, const int &steptype);//step according to the steptype=VELOCITY or ANGULAR.

double WoodSaxonZ(const double &x1, const double &x2, const double &b0);//calc the density profile in xy plane.
void Initial(const double &b,const double &tau0, const int &Num0, const int &Num1, const int &Num2);//give an initial value of U1.
void normalizedWS(const int &AAA);//CALC the constand c0 to normalize the woodsaxon distribution. 
//((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((





//
///////////////////////////////Main program/////////////////////////////////////////
int main(int argc, char** argv)
{
double tau0, b, tauend, ymax, Ifind, Jfind, Kfind,dtau;
int key, kskip;
int steptype = VELOCITY ;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&size);
MPI_Comm_rank(MPI_COMM_WORLD,&myid);
for(int i=myid; i<Num0; i+=size) eta[i] =(i+1)*deta;
for(int j=myid; j<Num1; j+=size) x1[j]  =(j+1)*dx1;
for(int k=myid; k<Num2; k+=size) x2[k]  =(k+1)*dx2;
if(myid==0){
normalizedWS(197);
b = 9.0;
tau0 = 0.6;
dtau = 0.001;
tauend= 12.6;
kskip = 100;
Initial(b, tau0, 100, 100, 100);

tau = tau0;
key = kskip;
ofstream f11("xyeta0_LnTYv1v2.dat");
ofstream f12("xyeta1_LnTYv1v2.dat");
ofstream f13("xyeta2_LnTYv1v2.dat");
ofstream f14("xyeta3_LnTYv1v2.dat");

ofstream f21("xetay0_LnTYv1v2.dat");
ofstream f22("xetay1_LnTYv1v2.dat");
ofstream f23("xetay2_LnTYv1v2.dat");
ofstream f24("xetay3_LnTYv1v2.dat");

ofstream f31("yetax0_LnTYv1v2.dat");
ofstream f32("yetax1_LnTYv1v2.dat");
ofstream f33("yetax2_LnTYv1v2.dat");
ofstream f34("yetax3_LnTYv1v2.dat");

while(tau <= tauend){
  cout<<tau<<"   "<<dtau<<"   "<<tauend<<endl;
	if(key != kskip) goto st;
  pushToStack(mystack,steptype);
  for(int j=0; j!=Num1; ++j)
	for(int k=0; k!=Num2; ++k){
	  int i=0;  //near eta=0 plane
	  f11<<left<<setw(15)<<tau<<setw(15) <<x1[j]<<setw(15) <<x2[k]<<setprecision(8)<<setw(15) <<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  i=29;  //near eta=0.9 plane
	  f12<<left<<setw(15)<<tau<<setw(15) <<x1[j]<<setw(15) <<x2[k]<<setw(15) <<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  i=59;  //near eta= 1.8 plane
	  f13<<left<<setw(15)<<tau<<setw(15) <<x1[j]<<setw(15) <<x2[k]<<setw(15) <<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  i=99;  //near eta=3.0 plane
	  f14<<left<<setw(15)<<tau<<setw(15) <<x1[j]<<setw(15) <<x2[k]<<setw(15) <<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	}

  for(int i=0; i!=Num0; ++i)
	for(int j=0; j!=Num1; ++j){
	  int k=0;  //near x2=0 plane
	  f21<<left<<setw(15)<<tau<<setw(15)<<x1[j]<<setw(15)<<eta[i]<<setw(15)<<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  k=29;  //near x2=4.5 plane
	  f22<<left<<setw(15)<<tau<<setw(15)<<x1[j]<<setw(15)<<eta[i]<<setw(15)<<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  k=59;  //near x2=9.0 plane
	  f23<<left<<setw(15)<<tau<<setw(15)<<x1[j]<<setw(15)<<eta[i]<<setw(15)<<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  k=99;  //near x2= 15 plane
	  f24<<left<<setw(15)<<tau<<setw(15)<<x1[j]<<setw(15)<<eta[i]<<setw(15)<<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;
	}

  for(int i=0; i!=Num0; ++i)
	for(int k=0; k!=Num2; ++k){
	  int j=0;
	  f31<<left<<setw(15)<<tau<<setw(15)<<x2[k]<<setw(15)<<eta[i]<<setw(15)<<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  j=29;
	  f32<<left<<setw(15)<<tau<<setw(15)<<x2[k]<<setw(15)<<eta[i]<<setw(15)<<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  j=59;
	  f33<<left<<setw(15)<<tau<<setw(15)<<x2[k]<<setw(15)<<eta[i]<<setw(15)<<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;

	  j=99;
	  f34<<left<<setw(15)<<tau<<setw(15)<<x2[k]<<setw(15)<<eta[i]<<setw(15)<<U1[0][i][j][k]<<setw(15)
	  <<U1[1][i][j][k]<<setw(15)<<U1[2][i][j][k]<<setw(15)<<U1[3][i][j][k]<<endl;
	}
	key = 0;

  popFromStack(mystack,steptype);
st:	  key += 1;
	  tau += dtau;
	  checkType(steptype);
	  step(dtau,steptype);
	  ymax = 0.0;
	  for( int i=0; i!=Num0; ++i)
	    for(int j=0; j!=Num1; ++j)
		for(int k=0; k!=Num2; ++k)
		  for(int l=0; l!=4; ++l){
			U1[l][i][j][k] = U2[l][i][j][k];
			if(U2[1][i][j][k] > ymax){
			  ymax = U2[1][i][j][k];
			  Ifind = i;
			  Jfind = j;
			  Kfind = k;
			}
		  }
	   cout<<"ymax= "<<setprecision(8)<< ymax<<"   "<<Ifind<<"   "<<Jfind<<"   "<<Kfind<<endl;

}
   		
		f11.close();
		f12.close();
		f13.close();
		f14.close();
		f21.close();
		f22.close();
		f23.close();
		f24.close();
		f31.close();
		f32.close();
		f33.close();
		f34.close();
}
MPI_Finalize();
return 0;					
}
///////////////////////////////End Main/////////////////////////////////////////////////////////

//*******************WoodSaxonZ(x1, x2, b00)************************************
double WoodSaxonZ(const double &x1, const double &x2, const double &b0){
double r,temp,f,a[10],b[10], z ;
double zk[5]={-0.9061798, -0.5386492, 0.0, 0.5386493, 0.9061798};//range(-1,1)
double Ak[5]={0.2369269,0.4786287,0.5688889,0.4786287,0.2369269};//weight
temp = 0.0;
for (int j=0; j!=10; ++j){
  a[j] = j*cut/10.0;
  b[j] = (j+1)*cut/10.0;
  for(int i=0; i!=5; ++i){
	z=(a[j]+b[j])/2.0 + (b[j]-a[j])/2.0* zk[i];
	r=sqrt(x1*x1+(x2-b0/2.0)*(x2-b0/2.0)+z*z);
	f=c0/(1.0+exp((r-RR)/a0));
	temp= temp +Ak[i]*f*(b[j]-a[j])/2.0;
	}
  }
return 2.0*temp;
/*There is a gauss intergration .
Change integration variable z to zk, z=(a+b)/2.0+(b-a)/2.0*zk;
Then \int_{a}^{b} dz 
=(b-a)/2.0 *\int_{-1}^{1} dzk 
=\sum (b-a)/2.0 *A_{i}*g(zk)
*/
}
//******************************************************************************
//&&&&&&&&&&&&&&&&&&& Initial() &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
void Initial(const double &b,const double &tau0, const int &Num0, const int &Num1, const int &Num2){
double a1, nxyz, eta, x1, x2, Ta, Tb, nbc, nwn, nmix,sig, fbc, fwn, sigmaz, eta0;
sigmaz = 0.3 ;
sig = 4.0;
fbc = 0.25; 
fwn = 0.75;
eta0 = 5.0;
a1 = 0.5;//
for(int j=0; j!=Num1; ++j)
	for(int k=0; k!=Num2; ++k){
	  x1=dx1*(j+1.0);
	  x2=dx2*(k+1.0);
	  Ta=WoodSaxonZ(x1, x2, b);
	  Tb=WoodSaxonZ(x1, x2, -b);
	  nbc = sig*Ta*Tb;
	  nwn= Ta*(1.0-exp(-sig*Tb))+Tb*(1.0-exp(-sig*Ta));
	  nmix = fbc*nbc +fwn*nwn;
		for(int i=0; i!=Num0; ++i){
// z=tau0*sinh(eta)
// nxyz=nmix/(sqrt(2*pi)*sigmaz)*exp(-z*z/(2*sigmaz*sigmaz))
		  eta = deta*(i+1.0);
		  nxyz=nmix*1.0/(1.0+exp(eta- eta0)/a1);
		  U1[0][i][j][k]= log(pow(nxyz,0.25));
		  U1[1][i][j][k]= eta;
	  	  U1[2][i][j][k]=0.0;
		  U1[3][i][j][k]=0.0;}
	}
		  

}
//&&&&&&&&&&&&&&&&&&&Initial end&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//%%%%%%%%%%%%%%%%%%normalizedWS(AAA) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void normalizedWS(const int &AAA){
double temp, f, pi, a[10], b[10], r, z;
double A= double(AAA);
double zk[5]={-0.9061798, -0.5386492, 0.0, 0.5386493, 0.9061798};
double Ak[5]={0.2369269,0.4786287,0.5688889,0.4786287,0.2369269};
temp = 0.0;
a0 = 0.54;
RR = 1.12*pow(A, 1.0/3.0) - 0.86*pow(A,-1.0/3.0);
pi=4.0*atan(1.0);

cut = 1.5*RR;
for (int j=0; j!=10; ++j){
  a[j] = j*cut/10.0;
  b[j] = (j+1)*cut/10.0;
  for(int i=0; i!=5; ++i){
	z=(a[j]+b[j])/2.0 + (b[j]-a[j])/2.0* zk[i];
	f=1.0/(1.0+exp((z-RR)/a0))*z*z;
	temp= temp +Ak[i]*f*(b[j]-a[j])/2.0;
	}
  }
c0 = A/(4.0 *pi *temp);

}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//checkType( steptype ) steptype == VELOCITY or ANGULAR;
void checkType(int &steptype)
{	if(steptype == VELOCITY && tau>=5.0){
	   double v, vmax=0.0;
		for(int i=0; i!=Num0; ++i)
		  for(int j=0; j!=Num1; ++j)
		    for(int k=0; k!=Num2; ++k){
		      v = sqrt(U1[2][i][j][k]*U1[2][i][j][k]+U1[3][i][j][k]*U1[3][i][j][k]);
		      if(v>vmax) vmax = v;
		    }
	   if(vmax > 0.99){
	       double v1, v2;
		for(int i=0; i!=Num0; ++i)
		  for(int j=0; j!=Num1; ++j)
		    for(int k=0; k!=Num2; ++k){
                    v1 = U1[2][i][j][k];
		    v2 = U1[3][i][j][k];
		    v=sqrt(v1*v1+v2*v2);
		U1[2][i][j][k] = 0.5*log((1+v)/(1-v));
		U1[3][i][j][k] = atan(v2/v1);    
		}
		steptype = ANGULAR;
	   }
       }
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//pushToStack(U1,mystack) if(steptype==ANGULAR) push alpha,theta to mystack, change U1(2,i,j,k) and U1(3,i,j,k) to v1,v2
//after print them out to files, change them back using popFromStack(U1,mystack);
void pushToStack(double mystack[2][Num0][Num1][Num2],const int& steptype)
{
  if(steptype == ANGULAR){
	for(int i=0; i!=Num0; ++i)
	  for(int j=0; j!=Num1; ++j)
	    for(int k=0; k!=Num2; ++k){
	      mystack[0][i][j][k] = U1[2][i][j][k];
	      mystack[1][i][j][k] = U1[3][i][j][k];
	      U1[2][i][j][k] = tanh( mystack[0][i][j][k] )*cos(mystack[1][i][j][k]);
	      U1[3][i][j][k] = tanh( mystack[0][i][j][k] )*sin(mystack[1][i][j][k]);
	    }
  }
}	


void popFromStack(double mystack[2][Num0][Num1][Num2],const int& steptype)
{
  if(steptype == ANGULAR){
	for(int i=0; i!=Num0; ++i)
	  for(int j=0; j!=Num1; ++j)
	    for(int k=0; k!=Num2; ++k){
	    U1[2][i][j][k]=mystack[0][i][j][k]  ;
	    U1[3][i][j][k]=mystack[1][i][j][k]  ;
	    }
  }


}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!!!!


//#####################matrixes()#############################################################################
void matrixes(const int &i, const int &j, const int &k, const double &cs2, const double &v1_or_alpha,const double &v2_or_theta, const double &eta,
	const double &y, double AD0[4][4], double AD1[4][4], double AD2[4][4], const int &steptype){
double v, alpha,theta,v1,v2, coshalpha2, crad, srad, trad, comm;
double sinT, cosT,sinhA, coshA,cothA, tanhA, coshA2, sechA2, temp[3][4][4];
	if(steptype == VELOCITY){
	  v1 = v1_or_alpha;
 	  v2 = v2_or_theta;
 	  v=sqrt(v1*v1+v2*v2);
	  alpha = 0.5*log((1.0+v)/(1.0-v));
	  
	}
	else if(steptype== ANGULAR){
	  alpha=v1_or_alpha;
	  theta=v2_or_theta;
	  v = tanh(alpha);
	  v1 = v*cos(theta);
	  v2 = v*sin(theta);
sinT=sin(theta);   cosT=cos(theta); sinhA=sinh(alpha); coshA=cosh(alpha); tanhA=sinhA/coshA; cothA=1.0/tanhA;
coshA2=coshA*coshA; sechA2=1.0/coshA2;	  
	}


if(v>=1) cout<<left<<"Alarm ,tachyonic "<<v1<<setw(15) <<v2<<setw(15) <<eta<<setw(15) <<y<<setw(15) <<i<<setw(15) <<j<<setw(15) <<k<<endl;
coshalpha2= 1.0/(1.0 - v*v);
crad= cosh(y-eta);
srad= sinh(y-eta);
trad= tanh(y-eta);
comm= 1.0+(1.0/cs2-1.0)*coshalpha2*crad*crad;
AD0[0][0]= (1.0/cs2-1.0)*coshalpha2*srad*crad/comm;
AD0[0][1]= coshalpha2/comm;
AD0[0][2]= 0.0;
AD0[0][3]= 0.0;

AD0[1][0]= (1.0/cs2-v*v)/comm;
AD0[1][1]= AD0[0][0];
AD0[1][2]= 0.0;
AD0[1][3]= 0.0;

AD0[2][0]=v1*(1.0-v*v)*trad/comm;
AD0[2][1]=-v1/comm;
AD0[2][2]=trad;
AD0[2][3]=0.0;

AD0[3][0]=v2*(1.0-v*v)*trad/comm;
AD0[3][1]=-v2/comm;
AD0[3][2]=0.0;
AD0[3][3]=trad;


AD1[0][0]=v1*(1.0/cs2-1.0)*coshalpha2*crad/comm;
AD1[0][1]=-v1*coshalpha2*srad/comm;
AD1[0][2]=coshalpha2*crad/comm;
AD1[0][3]=0.0;

AD1[1][0]=-v1*(1.0/cs2-1.0)*srad/comm;
AD1[1][1]=v1*(1.0/cs2-v*v)*coshalpha2*crad/comm;
AD1[1][2]=-srad/comm;
AD1[1][3]=0.0;

AD1[2][0]=(crad*(1.0/cs2+v2*v2*(1.0/cs2-1.0))-srad*trad*(1.0-v*v))/comm;
AD1[2][1]=v1*v1*srad/comm;
AD1[2][2]=v1*(1.0/crad+crad*(-1.0+coshalpha2*(1.0/cs2-1.0)))/comm;
AD1[2][3]=0.0;

AD1[3][0]=-(1.0/cs2-1.0)*v1*v2*crad/comm;
AD1[3][1]=v1*v2*srad/comm;
AD1[3][2]=-v2*crad/comm;
AD1[3][3]=v1/crad;

AD2[0][0]=v2*(1.0/cs2-1.0)*coshalpha2*crad/comm;
AD2[0][1]=-v2*coshalpha2*srad/comm;
AD2[0][2]=0.0;
AD2[0][3]=coshalpha2*crad/comm;

AD2[1][0]=-v2*(1.0/cs2-1.0)*srad/comm;
AD2[1][1]=v2*(1.0/cs2-v*v)*coshalpha2*crad/comm;
AD2[1][2]=0.0;
AD2[1][3]=-srad/comm;

AD2[2][0]=-v1*v2*(1.0/cs2-1.0)*crad/comm;
AD2[2][1]=v1*v2*srad/comm;
AD2[2][2]=v2/crad;
AD2[2][3]=-v1*crad/comm;

AD2[3][0]=(crad*(1.0/cs2+v1*v1*(1.0/cs2-1.0))-srad*trad*(1.0-v*v))/comm;
AD2[3][1]=v2*v2*srad/comm;
AD2[3][2]=0.0;
AD2[3][3]=v2*(1.0/crad+crad*(-1.0+coshalpha2*(1.0/cs2-1.0)))/comm;



if(steptype == ANGULAR)
{//change matrix AD to M^{-}.AD.M
    for(int i=0; i!=4; ++i)
	for(int j=0; j!=4; ++j)
	{temp[0][i][j] = AD0[i][j];
	 temp[1][i][j] = AD1[i][j];
	 temp[2][i][j] = AD2[i][j];
	}
double s[3][4][4];
 for(int id=0; id!=3; ++id)
 {
  s[id][0][0]=temp[id][0][0];
  s[id][0][1]=temp[id][0][1];
  s[id][0][2]=sechA2*(temp[id][0][2]*cosT + temp[id][0][3]*sinT);
  s[id][0][3]=tanhA *(temp[id][0][3]*cosT - temp[id][0][2]*sinT);

  s[id][1][0]=temp[id][1][0];
  s[id][1][1]=temp[id][1][1];
  s[id][1][2]=sechA2*(temp[id][1][2]*cosT + temp[id][1][3]*sinT);
  s[id][1][3]=tanhA *(temp[id][1][3]*cosT - temp[id][1][2]*sinT);

  s[id][2][0]=coshA2*(temp[id][2][0]*cosT + temp[id][3][0]*sinT);
  s[id][2][1]=coshA2*(temp[id][2][1]*cosT + temp[id][3][1]*sinT);  
  s[id][2][2]=temp[id][2][2]*cosT*cosT + sinT*cosT*(temp[id][3][2]+temp[id][2][3])+temp[id][3][3]*sinT*sinT;
  s[id][2][3]=sinhA*coshA*( temp[id][2][3]*cosT*cosT + sinT*cosT*(temp[id][3][3]-temp[id][2][2])-temp[id][3][2]*sinT*sinT );		
  
  s[id][3][0]=cothA*(temp[id][3][0]*cosT - temp[id][2][0]*sinT); 	
  s[id][3][1]=cothA*(temp[id][3][1]*cosT - temp[id][2][1]*sinT); 
  s[id][3][2]=1/(sinhA*coshA)*(temp[id][3][2]*cosT*cosT + sinT*cosT*(temp[id][3][3]-temp[id][2][2])-temp[id][2][3]*sinT*sinT);
  s[id][3][3]=temp[id][3][3]*cosT*cosT - sinT*cosT*(temp[id][3][2]+temp[id][2][3])+temp[id][2][2]*sinT*sinT;
  }
    for(int i=0; i!=4; ++i)
	for(int j=0; j!=4; ++j)
	{AD0[i][j]=s[0][i][j];
	 AD1[i][j]=s[1][i][j];
	 AD2[i][j]=s[2][i][j] ;
	}	
}  

}
//############################################################################################################


//@@@@@@@@@@@@@@@@@@@@@derivative@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
double derivative(const int& N1,const int& N2,const int& N3,const int& I,const int& J,const int& K,const int& steptype){
/*N1 = {1, 2} denotes 2 step of lax method
 *N2 = {eta, x1, x2} 
 *N3 = {lnT, y, v1, v2}   for steptype= VELOCITY
 *N3 = {lnT, y, alpha,theta} for steptype= ANGULAR
 */
//physical boundary y=0 at eta=0;  v1=0 at x1=0; v2=0 at x2=0; steptype= VELOCITY	
//physical boundary y=0 at eta=0;  theta = pi/2 at x1=0 ; theta =0 at x2=0; 
double deriv=0.0;
double pi=4.0*atan(1.0);
if(N1==0){ //for U1 or U2
  if(N2==0){  //*****respect to eta direction
    if(I>0 and I< (Num0-1))
	deriv = (U1[N3][I+1][J][K]-U1[N3][I-1][J][K])/(deta*2);
    else if(I==0){
	if(N3 !=1) deriv = (U1[N3][I+1][J][K]-U1[N3][I][J][K])/deta;
	else if(N3==1) //***** dy/deta near (eta=0, y=0)
	  deriv =(U1[N3][I+1][J][K]-0.0)/(deta*2); //why (2*deta)?
	}
    else if(I==(Num0-1))
	deriv = (U1[N3][I][J][K]-U1[N3][I-1][J][K])/deta;
	}
  
  else if(N2==1){//respect to x1 direction.
    if(J>0 and J<Num1-1) 
	deriv = (U1[N3][I][J+1][K]-U1[N3][I][J-1][K])/(dx1*2);
    else if(J==0){
      if(steptype == VELOCITY){
	if(N3 != 2) deriv=(U1[N3][I][J+1][K]-U1[N3][I][J][K])/dx1;
	else if(N3 ==2) deriv=(U1[N3][I][J+1][K]-0.0)/(dx1*2);//v1=0 at x1=0
	}
      else if(steptype== ANGULAR){
	if((N3!=3)&& !(N3==2 && K==0)) deriv=(U1[N3][I][J+1][K]-U1[N3][I][J][K])/dx1;
        else if(N3==3) deriv=(U1[N3][I][J][K]-pi/2.0)/(2*dx1);//theta = pi/2 at x1=0;
        else if(N3==2 && K==0) deriv=(U1[N3][I][J+1][K]-0.0)/(2*dx1); //alpha=0 at x1=0, x2=0;
      } 
     }
    else if(J==Num1-1)
	deriv=(U1[N3][I][J][K]-U1[N3][I][J-1][K])/dx1	;
  }

  else if(N2==2){//respect to x2 direction
    if(K>0 and K<Num2-1)
	deriv= (U1[N3][I][J][K+1]-U1[N3][I][J][K-1])/(2.0*dx2);
    else if(K==0){
	  if(steptype == VELOCITY){
	    if(N3!=3) deriv=(U1[N3][I][J][K+1]-U1[N3][I][J][K])/dx2;
	    else if(N3==3) deriv = (U1[N3][I][J][K+1]-0.0)/(2.0*dx2);
	  }
          else if(steptype == ANGULAR){
	   if((N3!=3)&& !(N3==2 && J==0)) deriv=(U1[N3][I][J][K+1]-U1[N3][I][J][K])/dx2;
           else if(N3==3) deriv=(U1[N3][I][J][K]-0.0)/(2*dx2);//theta=0 at x2=0;
	   else if(N3==2 && J==0) deriv=(U1[N3][I][J][K]-0.0)/(2*dx2);//alpha=0 at x1=0,x2=0;
          
	  }
	}

    else if(K==Num2-1)
	deriv= (U1[N3][I][J][K]-U1[N3][I][J][K-1])/dx2;
    }
  }

else if(N1==1){//two step method. this is the second step.
  if(N2==0){//correspond to eta direction
    if(I>0 and I<Num0-1)
	deriv = (U2[N3][I+1][J][K]-U2[N3][I-1][J][K])/(deta*2);
    else if(I==0){
	if(N3 !=1) deriv = (U2[N3][I+1][J][K]-U2[N3][I][J][K])/deta;
	else if(N3==1) deriv =(U2[N3][I+1][J][K]-0.0)/(deta*2); 
	}
    else if(I==Num0-1)
	deriv = (U2[N3][I][J][K]-U2[N3][I-1][J][K])/deta;
	}
  
  else if(N2==1){//correspond to x1 direction
    if(J>0 and J<Num0-1) 
	deriv = (U2[N3][I][J+1][K]-U2[N3][I][J-1][K])/(dx1*2);
    else if(J==0){
        if(steptype == VELOCITY){
	  if(N3 != 2) deriv=(U2[N3][I][J+1][K]-U2[N3][I][J][K])/dx1;
	  else if(N3 ==2) deriv=(U2[N3][I][J+1][K]-0.0)/(dx1*2);//
	 }
        else if(steptype== ANGULAR){
	  if((N3!=3)&& !(N3==2 && K==0)) deriv=(U2[N3][I][J+1][K]-U2[N3][I][J][K])/dx1;
          else if(N3==2 && K==0) deriv=(U2[N3][I][J+1][K]-0.0)/(2*dx1);
          else if(N3==3) deriv=(U2[N3][I][J][K]-pi/2.0)/(2*dx1);
        } 
     }

    else if(J==Num1-1)
	deriv=(U2[N3][I][J][K]-U2[N3][I][J-1][K])/dx1	;
  }

  else if(N2==2){//correspond to x2 direction
    if(K>0 and K<Num2-1)
	deriv= (U2[N3][I][J][K+1]-U2[N3][I][J][K-1])/(2.0*dx2);

    else if(K==0){
	  if(steptype == VELOCITY){
	    if(N3!=3) deriv=(U2[N3][I][J][K+1]-U2[N3][I][J][K])/dx2;
	    else if(N3==3) deriv = (U2[N3][I][J][K+1]-0.0)/(2.0*dx2);
	  }
          else if(steptype == ANGULAR){
	   if((N3!=3)&& !(N3==2 && J==0)) deriv=(U2[N3][I][J][K+1]-U2[N3][I][J][K])/dx2;
           else if(N3==3) deriv=(U2[N3][I][J][K]-0.0)/(2*dx2);
	   else if(N3==2 && J==0) deriv=(U2[N3][I][J][K]-0.0)/(2*dx2);          
	  }
	}
    
    else if(K==Num2-1)
	deriv= (U2[N3][I][J][K]-U2[N3][I][J][K-1])/dx2;
    }
  }

else cout<<"I catch nothing"<<endl;

return deriv;

}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


//^^^^^^^^^^^^^^^^^^^^^^^^step()^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void step(const double &dtau, const int &steptype){

double dydx[2][3][4],  change[4],cs2,v1_or_alpha, v2_or_theta, eta, y;
double AD0[4][4], AD1[4][4], AD2[4][4];
cs2=1.0/3.0;
for(int i=0; i!=Num0; ++i)
  for(int j=0; j!=Num1; ++j)
    for(int k=0; k!=Num2; ++k){
	  v1_or_alpha=U1[2][i][j][k];
	  v2_or_theta=U1[3][i][j][k];
	  eta=deta*(i+1);
	  y=U1[1][i][j][k];
	  matrixes(i,j,k,cs2,v1_or_alpha,v2_or_theta,eta,y,AD0,AD1,AD2,steptype);
		for(int jj=0; jj!=3; ++jj)
		  for(int kk=0; kk!=4; ++kk){
		  dydx[0][jj][kk]=derivative(0,jj,kk,i,j,k, steptype);
		  //jj={eta, x1, x2} and kk={lnT,y,v1_or_alpha,v2_or_theta}.
		  }
		for(int l=0; l!=4; ++l){
		  change[l]=0.0;
		  for(int m=0; m!=4; ++m){
		  change[l]=change[l]+AD0[l][m]*dydx[0][0][m]/tau+AD1[l][m]*dydx[0][1][m]+AD2[l][m]*dydx[0][2][m];
		  }
		U2[l][i][j][k]=U1[l][i][j][k]-change[l]*dtau/2.0;
		}
    }


for(int i=0; i!=Num0; ++i)
  for(int j=0; j!=Num1; ++j)
    for(int k=0; k!=Num2; ++k){
	  v1_or_alpha=U2[2][i][j][k];
	  v2_or_theta=U2[3][i][j][k];
	  eta=deta*(i+1);
	  y=U2[1][i][j][k];
	  matrixes(i,j,k,cs2,v1_or_alpha,v2_or_theta,eta,y,AD0,AD1,AD2,steptype);
		for(int jj=0; jj!=3; ++jj)
		  for(int kk=0; kk!=4; ++kk){
		  dydx[1][jj][kk]=derivative(1,jj,kk,i,j,k,steptype);
		  //jj={eta, x1, x2} and kk={lnT,y,v1_or_alpha,v2_or_theta}.
		  }
		for(int l=0; l!=4; ++l){
		  change[l]=0.0;
		  for(int m=0; m!=4; ++m){
		  change[l]=change[l]+AD0[l][m]*dydx[1][0][m]/tau+AD1[l][m]*dydx[1][1][m]+AD2[l][m]*dydx[1][2][m];
		  }
		U2[l][i][j][k]=U1[l][i][j][k]-change[l]*dtau;
		}
    }

}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

