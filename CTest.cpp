#include<iostream>
#include<iomanip>
#include<fstream>
#include"CEos.h"
#include"PhyConst.h"
#include<vector>
#include<algorithm>
#include<cassert>
#include<cmath>
#include<cstdlib>
int a[5][5]={0};
void f(int b[5][5]){

    for(int i=0; i<5; i++)
    	for(int j=0; j<5; j++){
	    b[i][j] += i*j;
	    }
}
using namespace PhyConst;

double e(double T){
  return 3.0*pi*pi/30.0*T*T*T*T/hbarc/hbarc/hbarc;
}

using namespace std;
int main(int argc, char **argv)
{
//CEos Eos(EOSL);
//
//ofstream f_eos("Eos.dat",ios_base::out);
//
//for(int i=0; i<=1000; i++){
//double epsi=i*1.0e-2;
//double epsip=(i+1)*1.0e-2;
//if(Eos.P(epsi,0.0)>Eos.P(epsip,0.0))cout<<"P(i)="<<Eos.P(epsi,0.0)<<" P(i+1)="<<Eos.P(epsip,0.0)<<endl;
////f_eos<<setw(15)<<epsi<<setw(15)<<Eos.P(epsi,0.0)<<endl;
//}
//
//f_eos.close();

//	int myints[]={2,3,1,4};
//	vector<int> a(myints, myints+4);
//	sort(a.begin(), a.end());
//	for(vector<int>::iterator it=a.begin(); it!=a.end(); ++it)
//		cout<<" "<<*it<<endl;

//  Test Eos

assert(argc==2);
double e=atof(argv[1]);
CEos Eos(EOSI);
cout<<"T(e="<<e<<")="<<Eos.T(e, 0.0)<<endl;
//
//cout.precision(10);
//cout<<e(0.130)<<endl;
//ofstream fout("testfilesize.dat");
//    for(int i=0; i<500; i++)
//    	for(int j=0; j<500; j++)
//	   fout<<double(i*j) <<" ";
//	   fout.close();
//
return 0;
}
