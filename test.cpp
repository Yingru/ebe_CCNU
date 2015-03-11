#include"CEos.h"
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cassert>
#include<sys/types.h>
#include<sys/stat.h>
#include<sys/wait.h>
#include<unistd.h>
#include"CAAcollision.h"
using namespace std;

int main(int argc, char** argv)
{

    /*
       CNucleus A1(197, 0.17, 6.38, 0.535);      //projectile nucleu
       CNucleus A2(197, 0.15, 6.38, 0.535);
       cout<<A1.Thickness(0,0)<<endl;
       cout<<A2.Thickness(0,0)<<endl;
       CAAcollision AACentral1(A1, A1, 4.0, 0.0, 0.95, 0.05, 2.95, 0.04);                //central collisions for kFactor calculation
       CAAcollision AACentral2(A2, A2, 4.0, 0.0, 0.95, 0.05, 2.95, 0.04);                //central collisions for kFactor calculation
       cout<<"A1 NCollision="<<AACentral1.NCollision(8.0,0.0,0.0)<<endl;
       cout<<"A2 NCollision="<<AACentral2.NCollision(8.0,0.0,0.0)<<endl;

       cout<<" int((0.68-0.6)/0.04)="<<int((0.68-0.6)/0.04)<<endl;
       return 0;
     */



    //      if(mkdir("mkdir_test",0777) ==-1)cout<<"mkdir erro!"<<endl; 


    pid_t pid;			
    int status;

    if( (pid=vfork())<0 ) cout<<"vfork error"<<endl;
    else if(pid==0){
        cout<<"vfork succeed!"<<endl;
        chdir("Reso3D/");
        if(execle("reso","reso",argv[1],(char *)0)==0)cout<<"exec succeed"<<endl;
        chdir("../");
        _exit(0);
    }
    else{
        if(waitpid(pid, &status, 0)==pid){
            cout<<"execle succeed"<<endl;
        }
        else{
            cout<<"Child not finished !"<<endl;
        }
    }
            cout<<"This is parent!"<<endl;

    /* 
       CEos Eos(EOSL2);
       double e;
       while(cin>>e){
       cout<<Eos.T(e,0);

       }
     */ 


    //ifstream fin0("results/dx0p2/dNdEta2.dat");
    //ifstream fin1("results/dx0p2/dNdEta1.dat");
    //
    //double temp;
    //while(!fin0.eof()){
    //    double N0, N1;
    //    fin0>>temp>>N0;
    //    fin1>>temp>>N1;
    //    cout<<N0/N1<<endl;
    //}
    //

    //Add a cvs test change here
    return 0;
}
