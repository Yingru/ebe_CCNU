#include"CHydro.h"
/*Default Construction function of hydro class*/
CHydro3D::CHydro3D() {
}//End Default Construction
////////////////////////////////////////////////////////////////////////////////
/* Calc the total entropy density at time step tau      */
inline double CHydro3D::STotal(){

    double stotal = 0.0;

    for(int i=NX0; i<=NX; i++)
        for(int j=NY0; j<=NY; j++)
            for(int k=NZ0; k<=NZ; k++){
                stotal += Eos.S(Ed(i,j,k),0.0)*tau*dx*dy*dz*Gamma(i,j,k);
            }
    return stotal;
}
//////////////////////////////////////////////////////////////////////

/*Hydrodynamics evolution */
void CHydro3D::Evolution(){

    /* used to calc Ed^{dt/2}, V^{dt/2}, P^{dt/2} */
    M3D TH00(NX0,NX, NY0,NY, NZ0,NZ);
    M3D TH01(NX0,NX, NY0,NY, NZ0,NZ);
    M3D TH02(NX0,NX, NY0,NY, NZ0,NZ);
    M3D TH03(NX0,NX, NY0,NY, NZ0,NZ);
    M3D THJT(NX0,NX, NY0,NY, NZ0,NZ);

    M3D Ed0(NX0,NX, NY0,NY, NZ0,NZ);
    M3D Vx0(NX0,NX, NY0,NY, NZ0,NZ);
    M3D Vy0(NX0,NX, NY0,NY, NZ0,NZ);
    M3D Vz0(NX0,NX, NY0,NY, NZ0,NZ);
    int NT=0;
    int DNT;//This value is different for two different hypersf calculation
    char file_prog[50];
    sprintf(file_prog,"time/time_%d.dat",IEVENT);
    ofstream fprog(file_prog);


    while( ED_MAX > setting.edec ){

        DNT=setting.NTD;


        if(NT%DNT==0){
            for(int i=NX0; i<=NX; i++)
                for(int j=NY0; j<=NY; j++)
                    for(int k=NZ0; k<=NZ; k++){
                        Ed0(i,j,k)=Ed(i,j,k);
                        Vx0(i,j,k)=Vx(i,j,k);
                        Vy0(i,j,k)=Vy(i,j,k);
                        Vz0(i,j,k)=Vz(i,j,k);
                    }//end for 
            double ST = STotal();
            fprog<<"time="<<tau<<" ED_MAX="<<ED_MAX<<" Stotal="<<ST<<endl;
        }




        if(NT%5==0){
/*            char file_EDXY_ETA0[50];
            sprintf(file_EDXY_ETA0,"results/event%d/ED_XY_%d.dat",IEVENT,NT);
            ofstream fout(file_EDXY_ETA0);
            fout.precision(3);

            for(int i=NX0; i<=NX; i++)
                for(int j=NY0; j<=NY; j++){
                	fout<<i<<' '<<j<<' '<<0<<' '<<Ed(i,j,0)<<endl;
            }
            fout.close();

            char file_EDXH_ETA0[50];
            sprintf(file_EDXH_ETA0,"results/event%d/ED_XH_%d.dat",IEVENT,NT);
            ofstream fout1(file_EDXH_ETA0);
            fout1.precision(3);

            for(int i=NX0; i<=NX; i++)
                for(int j=NZ0; j<=NZ; j++){
                	fout1<<i<<' '<<0<<' '<<j<<' '<<Ed(i,0,j)<<endl;
            }
            fout1.close();
*/

            Output(NT);
        }
        for(int i=NX0; i<=NX; i++)
            for(int j=NY0; j<=NY; j++)
                for(int k=NZ0; k<=NZ; k++){
                    TH00(i,j,k)=TT00(i,j,k);
                    TH01(i,j,k)=TT01(i,j,k);
                    TH02(i,j,k)=TT02(i,j,k);
                    TH03(i,j,k)=TT03(i,j,k);
                    THJT(i,j,k)=TJT(i,j,k);
                }//end for half time step variables

        //2nd order runge kutta method for source term 
        StepUpdate2(TH00,TH01,TH02,TH03,THJT,TT00,TT01,TT02,TT03,TJT,0.5*dt);
        StepUpdate2(TT00,TT01,TT02,TT03,TJT,TH00,TH01,TH02,TH03,THJT,    dt);

        tau += dt;  
        NT ++;

        //StepUpdate();

        //CalcHyperSf(Ed0, Vx0, Vy0, Vz0) ; 
        CalcHyperSf2(NT, Ed0, Vx0, Vy0, Vz0);


        bool outp=(abs(tau-1.60)<0.001) || 
            (abs(tau-3.60)<0.001) ||
            (abs(tau-5.60)<0.001) || 
            (abs(tau-7.60)<0.001) || 
            (abs(tau-9.60)<0.001) || 
            (abs(tau-11.60)<0.001) || 
            (abs(tau-13.60)<0.001);
        if(outp){
            char file_edx[128];
            char file_vxx[128];
            char file_vyy[128];
            char file_vzz[128];
            sprintf(file_edx,"results/event%d/edx_%f.dat",IEVENT,tau);
            sprintf(file_vxx,"results/event%d/vxx_%f.dat",IEVENT,tau);
            sprintf(file_vyy,"results/event%d/vyy_%f.dat",IEVENT,tau);
            sprintf(file_vzz,"results/event%d/vzz_%f.dat",IEVENT,tau);
            ofstream fedx(file_edx);
            ofstream fvxx(file_vxx);
            ofstream fvyy(file_vyy);
            ofstream fvzz(file_vzz);

            for(int i=NX0; i<=NX; i++){
                fedx<<i*dx<<' '<<Ed(i,0,0)<<endl;
                fvxx<<i*dx<<' '<<Vx(i,0,0)<<endl;
                fvyy<<i*dx<<' '<<Vy(0,i,0)<<endl;
            }
            for(int k=NZ0; k<=NZ; k++){
                fvzz<<k*dz<<' '<<Vz(0,0,k)*tau<<endl;
            }

            fedx.close();
            fvxx.close();
            fvyy.close();
            fvzz.close();
        }
    }

    fprog.close();

}//End Evolution


////////////////////////////////////////////////////////////////////////////////
/*Freeze out procedure*/
void CHydro3D::Freezeout2(){
    //calc dN/dEta by intergral over 
    ifstream f_hypsf("results/HyperSf.dat");
    ofstream f_progress("results/Frz_progress.dat");
    CD t_dec=0.137; //freeze out temperature
    CD m_pi=0.140; //pion mass Gev
    CD dphi=2.0*PhyConst::pi/40.0;
    CD dpt=0.1;
    CD deta=0.25;
    CD dof=3.0/pow(2.0*PhyConst::pi,3.0);
    CD Gevfm3=pow(1.0/PhyConst::hbarc,3.0);

    CI Neta=20;
    CI Npt=40;
    CI Nphi=40;
    double mt, pt, phi, Y, eta, gamma, vx, vy, vz;
    double e_dec,DA0, DA1, DA2, DA3, time, x, y, etas;
    double dNdEta[Neta]={0.0};
    double oldtime=0.0;

    while(!f_hypsf.eof()){
        f_hypsf>>e_dec>>DA0>>DA1>>DA2>>DA3>>vx>>vy>>vz>>time>>x>>y>>etas;
        if(time != oldtime) f_progress<<time<<"\n";
        oldtime=time;

        for(int k=0; k!=Neta; k++){
            eta=k*deta;
            for(int i=0; i!=Npt; i++){
                pt=i*0.1;
                mt=sqrt(pt*pt+m_pi*m_pi);
                double ptce=pt*cosh(eta);
                double dYdEta=ptce/sqrt(m_pi*m_pi+ptce*ptce);//   P/E  
                Y=atanh(dYdEta*tanh(eta));//  dY/deta=P/E
                double mtcy=mt*cosh(Y-etas);
                double mtsy=mt*sinh(Y-etas);
                for(int j=0; j!=Nphi; j++){
                    phi=j*dphi;
                    gamma=1.0/sqrt(max(1.0e-30,1.0-vx*vx-vy*vy-time*time*vz*vz));
                    double ptcp=pt*cos(phi);
                    double ptsp=pt*sin(phi);
                    double efkt=gamma*(mtcy-ptcp*vx-ptsp*vy-mtsy*time*vz)/t_dec;
                    double volum=mtcy*DA0-ptcp*DA1-ptsp*DA2-mtsy*DA3;
                    dNdEta[k]= dNdEta[k] +Gevfm3*dpt*dphi*dYdEta*volum/(exp(efkt)-1.0);
                }//end for phi
            }//end for pt
        }//end for eta 
    }//end while
    f_hypsf.close();
    f_progress.close();

    ofstream f_dNdY("results/dNdY.dat");
    for(int k=0; k!=Neta; k++){
        dNdEta[k]=dNdEta[k]*dof;
        eta=k*deta;
        f_dNdY<<setw(10)<<eta<<setw(20)<<dNdEta[k]<<"\n";
    }
    f_dNdY.close(); 

}//End Freezeout()
////////////////////////////////////////////////////////////////////////////////
/*Freeze out procedure using Gauss integral over phi direction and 
  Gauss-Laguerre integral over pt direction*/
void CHydro3D::Freezeout_xyz_half(){
    //calc dN/dEta by intergral over 
    char FHyperSF[60];
    sprintf(FHyperSF,"results/HyperSf%d.dat",IEVENT);
    ifstream f_hypsf(FHyperSF);

    char FProgress[60];
    sprintf(FProgress,"results/Frz_progress%d.dat",IEVENT);
    ofstream f_progress(FProgress);


    char FdNdPtdPhi[60];
    sprintf(FdNdPtdPhi,"results/dNdPtdPhi%d.dat",IEVENT);
    ofstream f_dNdPtdPhi(FdNdPtdPhi);

    //CD t_dec=0.130; //freeze out temperature
    CD t_dec=Eos.T(setting.edec,0.0);

    //CD m_pi=0.140; //pion mass Gev
    CD m_pi=0.0;
    CD dof=1.0/pow(2.0*PhyConst::pi,3.0);

    CD etalo=-5.0;
    CD etahi=5.0;
    CD ptlo =0.0;
    CD philo=0.0;
    CD phihi=0.5*PhyConst::pi;
    CD Gevfm3=pow(1.0/PhyConst::hbarc,3.0);

    CI Neta=41;
    CI Npt=15;
    CI Nphi=48;

    double *pgauss, *wgauss;
    double *pgaula, *wgaula;
    pgauss = gaulep48; wgauss=gaulew48;
    pgaula = gala15x;  wgaula=gala15w;

    double deta, dpt, dphi;
    double mt, pt, phi, Y, eta, gamma;
    double DA0, DA1, DA2, DA3, vx, vy, vz, time, x, y, etas;
    double dNdEdPtdPhi[Neta][Npt][Nphi]={0.0};
    double dNdYdPt[Npt][Nphi]={0.0};  //dNdYdPtdPhi at Y=0
    double dNdEta[Neta]={0.0};
    double oldtime=0.0;
    double invslope=0.25;
    deta=(etahi-etalo)/(Neta-1.0);

    while(!f_hypsf.eof()){
        f_hypsf>>DA0>>DA1>>DA2>>DA3>>vx>>vy>>vz>>etas;
        //	if(floor(time) != floor(oldtime) ) f_progress<<time<<"\n";
        //	oldtime=time;
        gamma=1.0/sqrt(max(1.0e-8,1.0-vx*vx-vy*vy-vz*vz));

        for(int k=0; k!=Neta; k++){
            eta=(k-Neta/2)*deta;                //eta from -5 to 5
            for(int i=0; i!=Npt; i++){
                pt=invslope*pgaula[i]+ptlo;
                mt=sqrt(pt*pt+m_pi*m_pi);
                double ptce=pt*cosh(eta);
                double dYdEta=ptce/sqrt(m_pi*m_pi+ptce*ptce);//   P/E  
                Y=atanh(dYdEta*tanh(eta));//  dY/deta=P/E
                double mtcy=mt*cosh(Y-etas);
                double mtsy=mt*sinh(Y-etas);
                for(int j=0; j!=Nphi/2; j++){
                    phi=(phihi-philo)*pgauss[j];
                    double ptcp=pt*cos(phi);
                    double ptsp=pt*sin(phi);
                    double efkt=gamma*(mtcy-ptcp*vx-ptsp*vy-mtsy*vz)/t_dec;
                    double volum=mtcy*DA0-ptcp*DA1-ptsp*DA2-mtsy*DA3;
                    double pvf=(efkt==0)?0.0:volum/(exp(efkt)-1.0);
                    if(k==20)dNdYdPt[i][j] += pvf;   //for Y==0
                    dNdEdPtdPhi[k][i][j] += dYdEta*pvf;
                }//end for phi
            }//end for pt
        }//end for Eta
    }//end while
    f_hypsf.close();
    f_progress.close();

    double dNdPt[Npt]={0.0};
    double dNdPhi[Nphi]={0.0};
    for(int k=0; k!=Neta; k++){
        //	eta=(etahi-etalo)*pgauss[k];
        eta=(k-Neta/2)*deta;
        for(int i=0; i!=Npt; i++){
            double sumphi=0.0;
            pt=invslope*pgaula[i]+ptlo;
            for(int j=0; j!=Nphi/2; j++){
                phi=(phihi-philo)*pgauss[j];
                sumphi = sumphi + wgauss[j]*dNdEdPtdPhi[k][i][j];//Integral over phi
                if(k==20){
                    dNdPt[i] += wgauss[j]*dNdYdPt[i][j];                     //dNdPt at y=0
                    dNdPhi[j]+= Gevfm3*pt*wgaula[i]*invslope*dNdYdPt[i][j];  //dNdPhi at y=0
                }
                f_dNdPtdPhi<<setw(15)<<eta<<setw(15)<<pt<<setw(15)<<phi<<setw(15)<<dNdEdPtdPhi[k][i][j]<<"\n";
            }
            sumphi = sumphi*(phihi-philo);
            if(k==20)dNdPt[i]=dNdPt[i]*(phihi-philo);
            dNdEta[k] += Gevfm3*sumphi*pt*wgaula[i]*invslope; //Integral over pt 
        }
    }
    f_dNdPtdPhi.close();

    char FdNdPt[60];
    sprintf(FdNdPt,"results/dNdPt%d.dat",IEVENT);
    ofstream f_dNdPt(FdNdPt);
    for(int i=0; i!=Npt; i++){
        pt=invslope*pgaula[i]+ptlo;
        dNdPt[i]=4.0*Gevfm3*dof*dNdPt[i]/(2.0*PhyConst::pi);
        f_dNdPt<<setw(10)<<pt<<setw(20)<<dNdPt[i]<<"\n";
    }
    f_dNdPt.close();

    char FdNdPhi[60];
    sprintf(FdNdPhi,"results/dNdPhi%d.dat",IEVENT);
    ofstream f_dNdPhi(FdNdPhi);
    for(int j=0; j!=Nphi; j++){
        phi=(phihi-philo)*pgauss[j];
        dNdPhi[j]=4.0*dof*dNdPhi[j];
        f_dNdPhi<<setw(10)<<phi<<setw(10)<<dNdPhi[j]<<"\n";
    }
    f_dNdPhi.close();

    char FdNdEta[60];
    sprintf(FdNdEta,"results/dNdEta%d.dat",IEVENT);
    ofstream f_dNdY(FdNdEta);
    for(int k=0; k!=Neta; k++){
        eta=(k-Neta/2)*deta;
        dNdEta[k]=4.0*dNdEta[k]*dof;
        f_dNdY<<setw(10)<<eta<<setw(20)<<dNdEta[k]<<"\n";
    }
    f_dNdY.close(); 
}//End Freezeout()


////////////////////////////////////////////////////////////////////////////////
/*Freeze out procedure using Gauss integral over phi direction and 
  Gauss-Laguerre integral over pt direction*/
void CHydro3D::Freezeout(){
    //calc dN/dEta by intergral over 
    char FHyperSF[60];
    sprintf(FHyperSF,"results/HyperSf%d.dat",IEVENT);
    ifstream f_hypsf(FHyperSF);

    char FProgress[60];
    sprintf(FProgress,"results/Frz_progress%d.dat",IEVENT);
    ofstream f_progress(FProgress);

    char FileOut[60];
    sprintf(FileOut,"results/dNdEtadPhi%d.dat",IEVENT);
    ofstream f_dNdEdPhi(FileOut);	

    char FdNdPtdPhi[60];
    sprintf(FdNdPtdPhi,"results/dNdPtdPhi%d.dat",IEVENT);
    ofstream f_dNdPtdPhi(FdNdPtdPhi);

    //CD t_dec=0.130; //freeze out temperature
    CD t_dec=Eos.T(setting.edec,0.0);

    CD m_pi=0.140; //pion mass Gev
    //CD m_pi=0.0;
    CD dof=1.0/pow(2.0*PhyConst::pi,3.0);

    CD etalo=-6.0;
    CD etahi=6.0;
    CD ptlo =0.0;
    CD philo=0.0;
    CD phihi=2.0*PhyConst::pi;
    CD Gevfm3=pow(1.0/PhyConst::hbarc,3.0);

    CI Neta=41;
    CI Npt=15;
    CI Nphi=48;

    double *pgauss, *wgauss;
    double *pgaula, *wgaula;
    pgauss = gaulep48; wgauss=gaulew48;
    pgaula = gala15x;  wgaula=gala15w;

    double deta, dpt, dphi;
    double mt, pt, phi, Y, eta, gamma;
    double DA0, DA1, DA2, DA3, vx, vy, vz, time, x, y, etas;
    double dNdEdPtdPhi[Neta][Npt][Nphi]={0.0};
    double dNdYdPt[Npt][Nphi]={0.0};  //dNdYdPtdPhi at Y=0
    double dNdEdPhi[Neta][Nphi]={0.0};
    double dNdEta[Neta]={0.0};
    double oldtime=0.0;
    double invslope=0.25;
    deta=(etahi-etalo)/(Neta-1.0);

    while(!f_hypsf.eof()){
        f_hypsf>>DA0>>DA1>>DA2>>DA3>>vx>>vy>>vz>>etas;
        //	if(floor(time) != floor(oldtime) ) f_progress<<time<<"\n";
        //	oldtime=time;
        gamma=1.0/sqrt(max(1.0e-8,1.0-vx*vx-vy*vy-vz*vz));

        for(int k=0; k!=Neta; k++){
            eta=(k-Neta/2)*deta;                //eta from -5 to 5
            for(int i=0; i!=Npt; i++){
                pt=invslope*pgaula[i]+ptlo;
                mt=sqrt(pt*pt+m_pi*m_pi);
                double ptce=pt*cosh(eta);
                double dYdEta=ptce/sqrt(m_pi*m_pi+ptce*ptce);//   P/E  
                Y=atanh(dYdEta*tanh(eta));//  dY/deta=P/E
                double mtcy=mt*cosh(Y-etas);
                double mtsy=mt*sinh(Y-etas);
                for(int j=0; j!=Nphi/2; j++){
                    phi=(phihi-philo)*pgauss[j];
                    double ptcp=pt*cos(phi);
                    double ptsp=pt*sin(phi);
                    double efkt=gamma*(mtcy-ptcp*vx-ptsp*vy-mtsy*vz)/t_dec;
                    double volum=mtcy*DA0-ptcp*DA1-ptsp*DA2-mtsy*DA3;
                    double pvf=(efkt==0)?0.0:volum/(exp(efkt)-1.0);
                    if(k==20)dNdYdPt[i][j] += pvf;   //for Y==0
                    dNdEdPtdPhi[k][i][j] += dYdEta*pvf;
                }//end for phi
            }//end for pt
        }//end for Eta
    }//end while
    f_hypsf.close();
    f_progress.close();

    double dNdPt[Npt]={0.0};
    double dNdPhi[Nphi]={0.0};
    for(int k=0; k!=Neta; k++){
        //	eta=(etahi-etalo)*pgauss[k];
        eta=(k-Neta/2)*deta;
        for(int i=0; i!=Npt; i++){
            double sumphi=0.0;
            pt=invslope*pgaula[i]+ptlo;
            for(int j=0; j!=Nphi/2; j++){
                phi=(phihi-philo)*pgauss[j];
                sumphi = sumphi + wgauss[j]*dNdEdPtdPhi[k][i][j];//Integral over phi
                if(k==20){
                    dNdPt[i] += wgauss[j]*dNdYdPt[i][j];                     //dNdPt at y=0
                    dNdPhi[j]+= Gevfm3*pt*wgaula[i]*invslope*dNdYdPt[i][j];  //dNdPhi at y=0
                }
                f_dNdPtdPhi<<setw(15)<<eta<<setw(15)<<pt<<setw(15)<<phi<<setw(15)<<dNdEdPtdPhi[k][i][j]<<"\n";
            }
            sumphi = sumphi*(phihi-philo);
            if(k==20)dNdPt[i]=dNdPt[i]*(phihi-philo);
            dNdEta[k] += Gevfm3*sumphi*pt*wgaula[i]*invslope; //Integral over pt 
        }
    }
    f_dNdPtdPhi.close();

    char FdNdPt[60];
    sprintf(FdNdPt,"results/dNdPt%d.dat",IEVENT);
    ofstream f_dNdPt(FdNdPt);
    for(int i=0; i!=Npt; i++){
        pt=invslope*pgaula[i]+ptlo;
        dNdPt[i]=Gevfm3*dof*dNdPt[i]/(2.0*PhyConst::pi);
        f_dNdPt<<setw(10)<<pt<<setw(20)<<dNdPt[i]<<"\n";
    }
    f_dNdPt.close();

    char FdNdPhi[60];
    sprintf(FdNdPhi,"results/dNdPhi%d.dat",IEVENT);
    ofstream f_dNdPhi(FdNdPhi);
    for(int j=0; j!=Nphi/2; j++){
        phi=(phihi-philo)*pgauss[j];
        dNdPhi[j]=dof*dNdPhi[j];
        f_dNdPhi<<setw(10)<<phi<<setw(10)<<dNdPhi[j]<<"\n";
    }
    f_dNdPhi.close();

    char FdNdEta[60];
    sprintf(FdNdEta,"results/dNdEta%d.dat",IEVENT);
    ofstream f_dNdY(FdNdEta);
    for(int k=0; k!=Neta; k++){
        eta=(k-Neta/2)*deta;
        dNdEta[k]=dNdEta[k]*dof;
        f_dNdY<<setw(10)<<eta<<setw(20)<<dNdEta[k]<<"\n";
    }
    f_dNdY.close(); 

    for(int k=0; k!=Neta; k++){
        for(int j=0; j!=Nphi; j++){
            for(int i=0; i!=Npt; i++){
                pt=invslope*pgaula[i]+ptlo;
                dNdEdPhi[k][j] += Gevfm3*pt*wgaula[i]*invslope*dNdEdPtdPhi[k][i][j];
            }
            f_dNdEdPhi<<" "<<dof*dNdEdPhi[k][j];
        }
        f_dNdEdPhi<<"\n";
    }

}//End Freezeout()

////////////////////////////////////////////////////////////////////////////////
/*Freeze out procedure*/
void CHydro3D::Freezeout3(){
    char FileOut[60];
    sprintf(FileOut,"results/dNdEtadPhi%d.dat",IEVENT);
    char FdNdEta[60];
    sprintf(FdNdEta,"results/dNdEta%d.dat",IEVENT);

    char FHyperSF[60];
    sprintf(FHyperSF,"results/HyperSf%d.dat",IEVENT);
    ifstream f_hypsf(FHyperSF);

    char FProgress[60];
    sprintf(FProgress,"results/Frz_progress%d.dat",IEVENT);
    //    ofstream f_progress(FProgress);

    ofstream f_dNdEdPhi(FileOut);	
    ofstream f_dNdEta(FdNdEta);
    //CD t_dec=0.137;//freeze out temperature
    CD t_dec=Eos.T(setting.edec,0.0);
    CD m_pi=0.139; //pion mass Gev
    //CD m_pi=0.0;
    CD dof=1.0/pow(2.0*PhyConst::pi,3.0);

    CD etalo=-6.0;
    CD etahi=6.0;
    CD ptlo =0.0;
    CD philo=0.0;
    CD phihi=2.0*PhyConst::pi;
    CD Gevfm3=pow(1.0/PhyConst::hbarc,3.0);

    CI Neta=41;
    CI Npt=15;
    CI Nphi=41;

    double invslope=1.0/12.0;

    double *pgauss, *wgauss;
    double *pgaula, *wgaula;
    pgauss = gaulep48; wgauss=gaulew48;
    pgaula = gala15x;  wgaula=gala15w;

    double deta, dpt, dphi;
    double mt, pt, phi, Y, eta, gamma;
    double DA0, DA1, DA2, DA3, vx, vy, vz, time, x, y, etas;
    double dNdEdPtdPhi[Neta][Npt][Nphi]={0.0};
    double dNdYdPt[Npt][Nphi]={0.0};  //dNdYdPtdPhi at Y=0
    double dNdEdPhi[Neta][Nphi]={0.0};
    double dNdEta[Neta]={0.0};
    deta = (etahi-etalo)/(Neta-1.0);
    dphi = (phihi-philo)/(Nphi-1.0);
    double oldtime;

    double ptce, dYdEta, mtcy, mtsy, ptcp, ptsp, efkt, volum, pvf;
    while(!f_hypsf.eof()){
        f_hypsf>>DA0>>DA1>>DA2>>DA3>>vx>>vy>>vz>>etas;
        //for(deque<SF>::iterator it=hypsf.begin(); it!=hypsf.end(); it ++){
        //cout<<"progress"<<it<<' '<<endl;
        //DA0 = (*it).DA0;  
        //DA1 = (*it).DA1;  
        //DA2 = (*it).DA2;  
        //DA3 = (*it).DA3;  
        //vx  = (*it).vx;
        //vy  = (*it).vy;
        //vz  = (*it).vz;
        //etas= (*it).etas;

        gamma=1.0/sqrt(max(1.0e-8,1.0-vx*vx-vy*vy-vz*vz));

        for(int k=0; k!=Neta; k++){
            eta=etalo+k*deta;                //eta from -5 to 5
            for(int i=0; i!=Npt; i++){
                pt=invslope*pgaula[i]+ptlo;
                mt=sqrt(pt*pt+m_pi*m_pi);
                ptce=pt*cosh(eta);
                dYdEta=ptce/sqrt(m_pi*m_pi+ptce*ptce);//   P/E  
                Y=atanh(dYdEta*tanh(eta));//  dY/deta=P/E

                mtcy=mt*cosh(Y-etas);
                mtsy=mt*sinh(Y-etas);
                for(int j=0; j!=Nphi-1; j++){
                    phi=philo + j*dphi;
                    ptcp=pt*cos(phi);
                    ptsp=pt*sin(phi);
                    efkt=gamma*(mtcy-ptcp*vx-ptsp*vy-mtsy*vz)/t_dec;
                    volum=mtcy*DA0-ptcp*DA1-ptsp*DA2-mtsy*DA3;
                    pvf=volum/(exp(efkt)-1.0);   //Gev.fm^3
                    if(abs(eta)<=1.0)dNdYdPt[i][j] += pvf*dYdEta*deta;   //for |eta|<=1.0
                    dNdEdPtdPhi[k][i][j] += dYdEta*pvf;
                }//end for phi
            }//end for pt
        }//end for Eta
    }//end while

    double dNdPt[Npt]={0.0};
    double dNdPt_dot_v2[Npt]={0.0};
    double dNdPhi[Nphi]={0.0};

    for(int k=0; k!=Neta; k++){
        eta = (k-Neta/2.0)*deta;
        for(int j=0; j!=Nphi-1; j++){
            for(int i=0; i!=Npt; i++){
                pt=invslope*pgaula[i]+ptlo;
                dNdEdPhi[k][j] += Gevfm3*pt*wgaula[i]*invslope*dNdEdPtdPhi[k][i][j];
                if(abs(eta)<=1.0){
                    dNdPt[i] += dphi*dNdYdPt[i][j];                     //dNdPt at y=0
                    phi=philo+j*dphi;
                    dNdPt_dot_v2[i] += cos(2.0*phi)*dphi*dNdYdPt[i][j];                     //dNdPt at y=0
                    dNdPhi[j]+= Gevfm3*pt*wgaula[i]*invslope*dNdYdPt[i][j];  //dNdPhi at y=0
                }
            }
            f_dNdEdPhi<<" "<<dof*dNdEdPhi[k][j];
            dNdEta[k] += dof*dNdEdPhi[k][j]*dphi;
        }
        f_dNdEdPhi<<"\n";
        f_dNdEta<<eta<<" "<<dNdEta[k]<<endl;
    }


    char FdNdPt[60];
    sprintf(FdNdPt,"results/dNdPt%d.dat",IEVENT);
    ofstream f_dNdPt(FdNdPt);

    char Fv2[60];
    sprintf(Fv2,"results/v2%d.dat",IEVENT);
    ofstream f_v2(Fv2);

    for(int i=0; i!=Npt; i++){
        pt=invslope*pgaula[i]+ptlo;
        double v2=dNdPt_dot_v2[i]/max(dNdPt[i],1.0E-33);
        dNdPt[i]=Gevfm3*dof*dNdPt[i]/(2.0*PhyConst::pi);
        f_dNdPt<<setw(10)<<pt<<setw(20)<<dNdPt[i]<<"\n";
        f_v2<<setw(10)<<pt<<setw(20)<<v2<<"\n";
    }
    f_dNdPt.close();
    f_v2.close();

    char FdNdPhi[60];
    sprintf(FdNdPhi,"results/dNdPhi%d.dat",IEVENT);
    ofstream f_dNdPhi(FdNdPhi);
    for(int j=0; j!=Nphi-1; j++){
        phi=philo + j*dphi;
        dNdPhi[j]=dof*dNdPhi[j];
        f_dNdPhi<<setw(10)<<phi<<setw(10)<<dNdPhi[j]<<"\n";
    }
    f_dNdPhi.close();


    }

    void CHydro3D::ReadParticles(){//Read particles' information
        double buff;
        char oneline[256];
        //ifstream fin("Reso3D/resoweak.dat");
        ifstream fin("s95p-v1/pdg05.dat");
        CParticle p;

        while(!fin.eof()){
            p.stable = 0;
            fin>>p.monval>>p.name>>p.mass>>p.width         \
                >>p.gspin>>p.baryon>>p.strange>>p.charm     \
                >>p.bottom>>p.gisospin>>p.charge>>p.decays;

        //    fin>>buff>>buff>>buff; //Not used in s95p-v1/pdg05.dat
            fin.getline(oneline,256);//Skip the current line
            int numpar;
            for(int k=0; k<p.decays; k++){
                fin>>buff>>numpar;
                fin.getline(oneline,256);//Skip the daughters' data
            }
            //if(numpar==1)p.stable=1;
            if(p.width < 1.0E-8)p.stable=1;

            particles.push_back(p);
        }
        particles.pop_back();  //The last line of the file is empty
        fin.close();


        CParticle antiB;//anti-baryon
        int N=particles.size();
        for(deque<CParticle>::size_type i=0; i!=N; i++){
            if(particles[i].baryon){
                antiB.monval = -particles[i].monval;
                antiB.name   = "A";
                antiB.name.append(particles[i].name);
                antiB.mass = particles[i].mass;
                antiB.width = particles[i].width;
                antiB.gspin = particles[i].gspin;
                antiB.baryon = -particles[i].baryon;
                antiB.strange= -particles[i].strange;
                antiB.charm = -particles[i].charm;
                antiB.bottom = -particles[i].bottom;
                antiB.gisospin=particles[i].gisospin;
                antiB.charge=-particles[i].charge;
                antiB.decays=particles[i].decays;
                antiB.stable=particles[i].stable;
                particles.push_back(antiB);
            }
        }
    }

    void CHydro3D::dNdYdptdphi(CParticle par){
        if(par.monval > 0){
        //Only cal positive particles, the anti-particles has the same spectra

        CD mass=par.mass;
        CD fermi_boson = (abs(par.baryon)) ? +1.0 : -1.0; // Fermi or boson
        //CD mu_b=(abs(par.baryon)) ? 0.130 : 0.0;              // To fix the pi/p ratio
        CD mu_b=0.0;

        int MC=abs(par.monval);
        //CD t_dec=0.137;//freeze out temperature
        CD t_dec=Eos.T(setting.edec,0.0);
        cout<<"t_dec="<<t_dec<<endl;

        CD g = par.gspin;
        CD dof = 1/pow(2.0*PhyConst::pi,3.0);
        char filename[128];
        if(par.monval>=0)sprintf(filename,"results/event%d/dNdYPtdPtdPhi_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename,"results/event%d/dNdYPtdPtdPhi_A%d.dat",IEVENT,MC);

        ofstream f_dNdYdPtdPhi(filename);

        CD rylo=-5.5;
        CD ryhi= 5.5;
        CD ptlo =0.0;
        CD philo=0.0;
        CD phihi=2.0*PhyConst::pi;
        CD Gevfm3=pow(1.0/PhyConst::hbarc,3.0);


        CI Nry=NETA;   //rapidity
        CI Npt=NPT;
        CI Nphi=NPHI;
        double PHI[Nphi];
        double WPHI[Nphi];//Gaussian weight for integral
        double COSPHI[Nphi];
        double SINPHI[Nphi];


        CD invslope=1.0/12.0;

        double *pgauss, *wgauss;
        double *pgaula, *wgaula;
        pgauss = gaulep48; wgauss=gaulew48;
        pgaula = gala15x;  wgaula=gala15w;

        double dry, dpt, dphi;
        double mt, pt, phi, Y, gamma;
        double DA0, DA1, DA2, DA3, vx, vy, vz, time, x, y, etas;
        double dNdYPtdPtdPhi[Nry][Npt][Nphi]={0.0};

        double ptce, dYdEta, mtcy, mtsy, ptcp, ptsp, efkt, volum, pvf;

        dry = (ryhi-rylo)/double(Nry-1.0);
        for(int i=0;i<24;i++){
            PHI[47-i] = PhyConst::pi*(pgauss[i] + 1.0);
            PHI[i] = PhyConst::pi*(1.0 - pgauss[i]);
        }//PHI from 0 to 2pi
        for(int j=0;j<Nphi;j++){
            COSPHI[j] = cos(PHI[j]);
            SINPHI[j] = sin(PHI[j]);
        }//It will be faster 


        cout<<"The size of the deque is:"<<hypsf.size()<<endl;
        for(deque<SF>::iterator it=hypsf.begin(); it!=hypsf.end(); it ++){
            int I_PRG=(it-hypsf.begin());

            if(I_PRG%10000==0)cout<<I_PRG<<endl;
            DA0 = (*it).DA0;  
            DA1 = (*it).DA1;  
            DA2 = (*it).DA2;  
            DA3 = (*it).DA3;  
            vx  = (*it).vx;
            vy  = (*it).vy;
            vz  = (*it).vz;
            etas= (*it).etas;

            gamma=1.0/sqrt(max(1.0e-30,1.0-vx*vx-vy*vy-vz*vz));

            for(int k=0; k!=Nry; k++){
                Y= rylo + dry*k;
                for(int i=0; i!=Npt; i++){
                    pt=invslope*pgaula[i]+ptlo;
                    mt=sqrt(pt*pt+mass*mass);

                    mtcy=mt*cosh(Y-etas);
                    mtsy=mt*sinh(Y-etas);
                    for(int j=0; j!=Nphi; j++){
                        ptcp=pt*COSPHI[j];
                        ptsp=pt*SINPHI[j];
                        efkt=(gamma*(mtcy-ptcp*vx-ptsp*vy-mtsy*vz)-mu_b)/t_dec;
                        volum=mtcy*DA0-ptcp*DA1-ptsp*DA2-mtsy*DA3;
                        //pvf=volum/(exp(efkt)-1.0);   //Gev.fm^3
                        pvf=volum/(exp(efkt)+fermi_boson);   //Gev.fm^3
                        dNdYPtdPtdPhi[k][i][j] += g*pvf;
                    }//end for phi
                }//end for pt
            }//end for Y
        }//end for deque

        for(int k=0; k!=Nry; k++){
            for(int i=0; i!=Npt; i++){
                for(int j=0; j!=Nphi; j++){
                    f_dNdYdPtdPhi<<dNdYPtdPtdPhi[k][i][j]<<' ';
                }
                f_dNdYdPtdPhi<<'\n';
            }
        }
        f_dNdYdPtdPhi.close();

        double dNdY[Nry]={0.0};
        char filename_dNdY[128];
        if(par.monval>=0)sprintf(filename_dNdY,"results/event%d/dNdY_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdY,"results/event%d/dNdY_A%d.dat",IEVENT,MC);
        ofstream f_dNdY(filename_dNdY);
        for(int k=0; k!=Nry; k++){
            Y=(k-Nry/2)*dry;
            for(int i=0; i!=Npt; i++){
                double sumphi=0.0;
                pt=invslope*pgaula[i]+ptlo;
                for(int j=0; j!=Nphi/2; j++){
                    sumphi = sumphi + wgauss[j]*(dNdYPtdPtdPhi[k][i][j]+dNdYPtdPtdPhi[k][i][Nphi-1-j]);//Integral over phi
                }
                sumphi = sumphi*0.5*(phihi-philo);
                dNdY[k] += dof*Gevfm3*sumphi*pt*wgaula[i]*invslope; //Integral over pt 
            }
            f_dNdY<<Y<<' '<<dNdY[k]<<endl;
        }
        f_dNdY.close();


        //    CParticle anti_baryon;  
        //    if(par.baryon){
        //        anti_baryon.monval=-par.monval;
        //        char fantib[64];
        //        sprintf(fantib,"results/dNdYPtdPtdPhi_A%d.dat",abs(par.monval));
        //        ofstream f_AdNdYdPtdPhi(fantib);
        //        for(int k=0; k!=Nry; k++){
        //            for(int i=0; i!=Npt; i++)
        //                for(int j=0; j!=Nphi; j++){
        //                    f_AdNdYdPtdPhi<<dNdYPtdPtdPhi[k][i][j]<<' ';
        //                }
        //            f_AdNdYdPtdPhi<<'\n';
        //        }
        //        f_AdNdYdPtdPhi.close();
        //    }/*The anti baryon has the same spactra with baryon 
        //       we don't need to re-calculation */

    }
}


    void CHydro3D::dNdEta(CI & DIR_OR_RESO, CParticle &par){//calc dNdEta for particle i
        //DIR_OR_RESO = 0 for direct produced particles
        //DIR_OR_RESO = 1 for after resonance decay
        CD mass=par.mass;
        int MC=abs(par.monval);
        CD dof=1.0/pow(2.0*PhyConst::pi,3.0);
        char filename[128];

        double rylo, ryhi, etalo, etahi;

        if(DIR_OR_RESO == 1){
            rylo=-5.5;   
            ryhi= 5.5;
            if(par.monval>=0)sprintf(filename,"results/event%d/spec_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename,"results/event%d/spec_A%d.dat",IEVENT,MC);
        }
        else{
            rylo=-5.5;   
            ryhi= 5.5;
            if(par.monval>=0)sprintf(filename,"results/event%d/dNdYPtdPtdPhi_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename,"results/event%d/dNdYPtdPtdPhi_A%d.dat",IEVENT,MC);
        }

        ifstream f_dNdYdPtdPhi(filename);


        CD ptlo =0.0;
        CD philo=0.0;
        CD phihi=2.0*PhyConst::pi;
        CD Gevfm3=pow(1.0/PhyConst::hbarc,3.0);

        CI Nry=NETA;   //rapidity
        CI Npt=NPT;
        CI Nphi=NPHI;

        double PHI[Nphi];
        double WPHI[Nphi];//Gaussian weight for integral
        double COSPHI[Nphi];
        double SINPHI[Nphi];
        double RY[Nry];

        //double invslope=0.25;
        double invslope=1.0/12.0;

        double *pgauss, *wgauss;
        double *pgaula, *wgaula;
        pgauss = gaulep48; wgauss=gaulew48;
        pgaula = gala15x;  wgaula=gala15w;

        double dry, dpt, dphi;
        double mt, pt, phi, Y, gamma;
        double dNdYPtdPtdPhi[Nry][Npt][Nphi]={0.0};
        double dNdEtaPtdPtdPhi[Nry][Npt][Nphi]={0.0};

        double ptce, dYdEta, mtcy, mtsy, ptcp, ptsp, efkt, volum, pvf;

        dry = (ryhi-rylo)/double(Nry-1.0);
        for(int i=0; i<Nry; i++){
            RY[i] = rylo + i*dry;
        }

        for(int i=0;i<24;i++){
            PHI[47-i] = PhyConst::pi*(pgauss[i] + 1.0);
            PHI[i] = PhyConst::pi*(1.0 - pgauss[i]);
        }//PHI from 0 to 2pi
        for(int k=0; k!=Nry; k++){
            for(int i=0; i!=Npt; i++){
                for(int j=0; j!=Nphi; j++){
                    f_dNdYdPtdPhi>>dNdYPtdPtdPhi[k][i][j];

// Qin: Reading differential spectra here
                }
            }
        }
        f_dNdYdPtdPhi.close();

        /********Calc dNdY and dNdPt ************************/
        double dNdYPtdPt[Nry][Npt]={0.0};
        double v2_dot_dN[Nry][Npt]={0.0};

        double dNdY[Nry]={0.0};

        char filename_dNdY[128];

        if(DIR_OR_RESO == 1){ //After resonance decay
            if(par.monval>=0)sprintf(filename_dNdY,"results/event%d/dNdY_Reso_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdY,"results/event%d/dNdY_Reso_A%d.dat",IEVENT,MC);
        }
        else{
            if(par.monval>=0)sprintf(filename_dNdY,"results/event%d/dNdY_Dirc_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdY,"results/event%d/dNdY_Dirc_A%d.dat",IEVENT,MC);
        }
        ofstream f_dNdY(filename_dNdY);
        for(int k=0; k!=Nry; k++){
            for(int i=0; i!=Npt; i++){
                double sumphi=0.0;
                pt=invslope*pgaula[i]+ptlo;
                for(int j=0; j!=Nphi/2; j++){

                    sumphi = sumphi + wgauss[j]*(dNdYPtdPtdPhi[k][i][j]+dNdYPtdPtdPhi[k][i][Nphi-1-j]);//Integral over phi
                    dNdYPtdPt[k][i] += wgauss[j]*(dNdYPtdPtdPhi[k][i][j]+dNdYPtdPtdPhi[k][i][Nphi-1-j]);//Integral over phi
                    v2_dot_dN[k][i] += wgauss[j]*(cos(2.0*PHI[j])*dNdYPtdPtdPhi[k][i][j]+
                            cos(2.0*PHI[Nphi-1-j])*dNdYPtdPtdPhi[k][i][Nphi-1-j]);//Integral over phi
                }
                dNdYPtdPt[k][i] = dNdYPtdPt[k][i]*0.5*(phihi-philo);
                v2_dot_dN[k][i] = v2_dot_dN[k][i]*0.5*(phihi-philo);
                sumphi = sumphi*0.5*(phihi-philo);
                dNdY[k] += dof*Gevfm3*sumphi*pt*wgaula[i]*invslope; //Integral over pt 
            }
            Y=RY[k];
            f_dNdY<<Y<<' '<<dNdY[k]<<endl;
        }
        f_dNdY.close();
        ///////////////////////////////End dNdY()/////////////////////////////////////////

        /////////////////////////////dNdPt()//////////////////////////////////////////////////
        CD CY_L=-1.3;                //low bound for integral over central rapidity
        CD CY_H= 1.3;                //high bound for integral over central rapidity
        CI CNY = NPHI;               //Using gaussian integral with 48 nodes.
        double CY[CNY];              //Yi nodes from [-1.3,1.3]
        double CdNdYPtdPt[NPT][CNY]; //interpolation value at Yi
        double Cv2_dot_dN[NPT][CNY]; //v2*dNdptdy
        double dNdPt[NPT]={0.0};
        double v2dNdPt[NPT]={0.0};   //v2*dNdpt
        char filename_dNdPt[128];
        char filename_v2[128];

        for(int i=0;i<CNY/2;i++){
            CY[i] = CY_L + CY_H*(1.0 - pgauss[i]);      //from -1.3 to 0
            CY[CNY-i-1] = CY_L + CY_H*(pgauss[i] + 1.0);      // from 0 to 1.3
        }// CY from -1.3 to 1.3

        if(DIR_OR_RESO == 1){ //After resonance decay
            if(par.monval>=0)sprintf(filename_dNdPt,"results/event%d/dNdPt_Reso_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdPt,"results/event%d/dNdPt_Reso_A%d.dat",IEVENT,MC);
        }
        else{
            if(par.monval>=0)sprintf(filename_dNdPt,"results/event%d/dNdPt_Dirc_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdPt,"results/event%d/dNdPt_Dirc_A%d.dat",IEVENT,MC);
        }
        ofstream f_dNdPt(filename_dNdPt);

        if(DIR_OR_RESO == 1){ //After resonance decay
            if(par.monval>=0)sprintf(filename_v2,"results/event%d/v2_Reso_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_v2,"results/event%d/v2_Reso_A%d.dat",IEVENT,MC);
        }
        else{
            if(par.monval>=0)sprintf(filename_v2,"results/event%d/v2_Dirc_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_v2,"results/event%d/v2_Dirc_A%d.dat",IEVENT,MC);
        }
        ofstream f_v2(filename_v2);

        for(int i=0; i!=Npt; i++){

            for(int j=0; j!=CNY; j++){
                Y=CY[j];
                int IY=1;
                while(Y>RY[IY] && IY<Nry-1){
                    IY ++;
                }

                assert(IY>=0 && IY<=Nry-1);
                double Y1=RY[IY-1];
                double Y2=RY[IY];
                double FY1=dNdYPtdPt[IY-1][i];
                double FY2=dNdYPtdPt[IY][i];
                double a = (Y-Y1)/(Y2-Y1);
                CdNdYPtdPt[i][j] = (1-a)*FY1 + a*FY2;

                FY1=v2_dot_dN[IY-1][i];
                FY2=v2_dot_dN[IY][i];
                Cv2_dot_dN[i][j] = (1-a)*FY1 + a*FY2;
            }//Interpolatin to calculate CdNdYPtdPt at center rapidity

            for(int j=0;j<CNY/2;j++){
                dNdPt[i] += wgauss[j]*(CdNdYPtdPt[i][j]+CdNdYPtdPt[i][CNY-1-j]);
                v2dNdPt[i] += wgauss[j]*(Cv2_dot_dN[i][j]+Cv2_dot_dN[i][CNY-1-j]);
            }
            dNdPt[i] = dof*Gevfm3*0.5*(CY_H-CY_L)*dNdPt[i];//Do gaussian integral to cal dNdPt
            v2dNdPt[i]=dof*Gevfm3*0.5*(CY_H-CY_L)*v2dNdPt[i];
            double v2=v2dNdPt[i]/dNdPt[i];

            pt=invslope*pgaula[i]+ptlo;
            f_dNdPt<<pt<<' '<<1.0/(2.0*PhyConst::pi)/(CY_H-CY_L)*dNdPt[i]<<endl;
            f_v2<<pt<<' '<<v2<<endl;
            par.PT[i]=pt;
            par.dNdPt[i]=1.0/(2.0*PhyConst::pi)/(CY_H-CY_L)*dNdPt[i];
        }
        f_dNdPt.close();
        /*************End dNdPt********************************/

        /***************dNdPhi***************************/
        double dNdPhi[Nphi]={0.0};
        double dNdPhidY[Nphi][Nry]={0.0};
        char filename_dNdPhi[128];
        if(DIR_OR_RESO == 1){ //After resonance decay
            if(par.monval>=0)sprintf(filename_dNdPhi,"results/event%d/dNdPhi_Reso_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdPhi,"results/event%d/dNdPhi_Reso_A%d.dat",IEVENT,MC);
        }
        else{
            if(par.monval>=0)sprintf(filename_dNdPhi,"results/event%d/dNdPhi_Dirc_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdPhi,"results/event%d/dNdPhi_Dirc_A%d.dat",IEVENT,MC);
        }
        ofstream f_dNdPhi(filename_dNdPhi);
        for(int j=0; j!=Nphi; j++){
            for(int k=0; k!=Nry; k++){
                Y=RY[k];
                for(int i=0; i!=Npt; i++){
                    pt=invslope*pgaula[i]+ptlo;
                    dNdPhidY[j][k] += dof*Gevfm3*dNdYPtdPtdPhi[k][i][j]*pt*wgaula[i]*invslope; //Integral over pt 
                }//pt integral
                dNdPhi[j] += dry*dNdPhidY[j][k];
            }//rapidity integral
            f_dNdPhi<<PHI[j]<<' '<<dNdPhi[j]<<endl;
            par.PHI[j] = PHI[j];
            par.dNdPhi[j] = dNdPhi[j];
        }
        f_dNdPhi.close();


        //Calc dNdEta *****************************/ 
        char filename_dNdEta[128];
        if(DIR_OR_RESO == 1){ //After resonance decay
            if(par.monval>=0)sprintf(filename_dNdEta,"results/event%d/dNdEta_Reso_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdEta,"results/event%d/dNdEta_Reso_A%d.dat",IEVENT,MC);
        }
        else{
            if(par.monval>=0)sprintf(filename_dNdEta,"results/event%d/dNdEta_Dirc_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdEta,"results/event%d/dNdEta_Dirc_A%d.dat",IEVENT,MC);
        }
        ofstream f_dNdEta(filename_dNdEta);


        char filename_dNdEtaPtdPtdPhi[128];
        if(DIR_OR_RESO == 1){ //After resonance decay
            etalo=-5.5;
            etahi= 5.5;
            if(par.monval>=0)sprintf(filename_dNdEtaPtdPtdPhi,"results/event%d/dNdEtaPtdPtdPhi_Reso_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdEtaPtdPtdPhi,"results/event%d/dNdEtaPtdPtdPhi_Reso_A%d.dat",IEVENT,MC);
        }
        else{
            etalo=-5.5;
            etahi= 5.5;
            if(par.monval>=0)sprintf(filename_dNdEtaPtdPtdPhi,"results/event%d/dNdEtaPtdPtdPhi_Dirc_%d.dat",IEVENT,MC);
            else if(par.monval<0)sprintf(filename_dNdEtaPtdPtdPhi,"results/event%d/dNdEtaPtdPtdPhi_Dirc_A%d.dat",IEVENT,MC);
        }
        ofstream f_dNdEtaPtdPtdPhi(filename_dNdEtaPtdPtdPhi);

        double deta=(etahi-etalo)/double(NETA-1.0);
        double dNdEta[NETA]={0.0};
        for(int k=0; k!=NETA; k++){
            double eta=etalo+k*deta;
            for(int i=0; i!=Npt; i++){
                pt=invslope*pgaula[i]+ptlo;
                mt=sqrt(pt*pt+mass*mass);
                double ptce=pt*cosh(eta);
                double dYdEta=ptce/sqrt(mass*mass+ptce*ptce);//   P/E  
                Y=atanh(dYdEta*tanh(eta));//  dY/deta=P/E
                if(Y<rylo){
                    cout<<"Y="<<Y<<endl;	
                    Y=rylo;
                }
                else if(Y>ryhi){
                    cout<<"Y="<<Y<<endl;	
                    Y=ryhi;
                }
                int IY=1;
                while(Y>RY[IY] && IY<Nry-1){
                    IY ++;
                }

                assert(IY>=0 && IY<=Nry-1);
                for(int j=0; j!=Nphi; j++){
                    double Y1= RY[IY-1];
                    double Y2= Y1 + dry;
                    double FY1=dNdYPtdPtdPhi[IY-1][i][j];
                    double FY2=dNdYPtdPtdPhi[IY][i][j];
                    double a=(Y-Y1)/(Y2-Y1);
                    dNdEtaPtdPtdPhi[k][i][j] = ((1.0-a)*FY1+a*FY2)*dYdEta;
                    f_dNdEtaPtdPtdPhi<<dNdEtaPtdPtdPhi[k][i][j]<<' ';
                }
            }
            f_dNdEtaPtdPtdPhi<<endl;
        }

        f_dNdEtaPtdPtdPhi.close();

        for(int k=0; k!=NETA; k++){
            double eta=etalo+k*deta;
            for(int i=0; i!=Npt; i++){
                pt=invslope*pgaula[i]+ptlo;
                double sumphi=0.0;
                for(int j=0; j!=Nphi/2; j++){
                    sumphi = sumphi + wgauss[j]*(dNdEtaPtdPtdPhi[k][i][j]+dNdEtaPtdPtdPhi[k][i][Nphi-1-j]);//Integral over phi
                }
                sumphi = sumphi*0.5*(phihi-philo);
                dNdEta[k] += dof*Gevfm3*sumphi*pt*wgaula[i]*invslope; //Integral over pt 
            }
            f_dNdEta<<eta<<' '<<dNdEta[k]<<endl;
            par.ETA[k]=eta;
            par.dNdEta[k]=dNdEta[k];
        }
        f_dNdEta.close();
        /***********************************************************/
    }



    /**calculate dNdEta for all the charged particles*/
    void CHydro3D::dNdEta_Charged(){
        double Eta[NETA];
        double dNdEta_Charged[NETA]={0.0};
        double dNdEtaPtdPtdPhi_Charged[NETA][NPT][NPHI]={0.0};

        /*dNdEta for charged particles */
        char fname[128];
        sprintf(fname,"results/event%d/dNdEta_Charged.dat",IEVENT);
        ofstream f_dNdEta_Charged(fname);
        /////////////////////////////////

        /*dNdEtaPtdPtdPhi for charged particles */
        sprintf(fname,"results/event%d/dNdEtaPtdPtdPhi_Charged.dat",IEVENT);
        ofstream f_dNdEtaPtdPtdPhi_Charged(fname);
        /////////////////////////////////////////

        for(deque<CParticle>::size_type i=0; i!=particles.size();++i){

            if(particles[i].stable && abs(particles[i].charge)){
                for(int k=0; k!=NETA; k++){ 
                    Eta[k] = particles[i].ETA[k];
                    dNdEta_Charged[k] += particles[i].dNdEta[k];
                }

                char fcname[128];
                sprintf(fcname,"results/event%d/spec_%d.dat",IEVENT,particles[i].monval);
                ifstream fin(fcname);

                for(int k=0; k!=NETA; k++){ 
                    for(int i=0; i!=NPT; i++){
                        for(int j=0; j!=NPHI; j++){
                            double dN;
                            fin>>dN;
                            dNdEtaPtdPtdPhi_Charged[k][i][j] += dN;
                        }
                    }
                }//Calc the dNdEtaPtdPtdPhi for all charged particles

            }
        }

        for (int i=0; i != NETA; i++){
            f_dNdEta_Charged<<Eta[i]<<' '<<dNdEta_Charged[i]<<endl;
        }
        f_dNdEta_Charged.close();

        for(int k=0; k!=NETA; k++){ 
            for(int i=0; i!=NPT; i++){
                for(int j=0; j!=NPHI; j++){
                    f_dNdEtaPtdPtdPhi_Charged<<dNdEtaPtdPtdPhi_Charged[k][i][j]<<' ';
                }
            }
            f_dNdEtaPtdPtdPhi_Charged<<endl;
        }//output the dNdEtaPtdPtdPhi for all charged particles
        f_dNdEtaPtdPtdPhi_Charged.close();
    }





    void CHydro3D::Freezeout4()
    {
        char FileOut[60];
        sprintf(FileOut,"results/dNdEtadPhi%d.dat",IEVENT);
        char FdNdEta[60];
        sprintf(FdNdEta,"results/dNdEta%d.dat",IEVENT);

        char FHyperSF[60];
        sprintf(FHyperSF,"results/HyperSf%d.dat",IEVENT);
        ifstream f_hypsf(FHyperSF);

        char FProgress[60];
        sprintf(FProgress,"results/Frz_progress%d.dat",IEVENT);
        //    ofstream f_progress(FProgress);

        ofstream f_dNdEdPhi(FileOut);	
        ofstream f_dNdEta(FdNdEta);
        CD t_dec=0.130;//freeze out temperature
        //CD m_pi=0.140; //pion mass Gev
        CD m_pi=0.0; //pion mass Gev
        CD dof=2.0/pow(2.0*PhyConst::pi,3.0);

        CD etalo=-10.0;
        CD etahi=10.0;
        CD ptlo =1.0;
        CD pthi =3.0;
        CD philo=0.0;
        CD phihi=2.0*PhyConst::pi;
        CD Gevfm3=pow(1.0/PhyConst::hbarc,3.0);

        CI Neta=41;
        CI Npt=16;
        CI Nphi=41;

        double invslope=0.25;

        double *pgauss, *wgauss;
        double *pgaula, *wgaula;
        pgauss = gaulep16; wgauss=gaulew16;
        pgaula = gala15x;  wgaula=gala15w;

        double deta, dpt, dphi;
        double mt, pt, phi, Y, eta, gamma;
        double DA0, DA1, DA2, DA3, vx, vy, vz, time, x, y, etas;
        double dNdEdPtdPhi[Neta][Npt][Nphi]={0.0};
        double dNdEdPhi[Neta][Nphi]={0.0};
        double dNdEta[Neta]={0.0};
        deta = (etahi-etalo)/(Neta-1.0);
        dphi = (phihi-philo)/(Nphi-1.0);
        double oldtime;

        double ptce, dYdEta, mtcy, mtsy, ptcp, ptsp, efkt, volum, pvf;
        while(!f_hypsf.eof()){
            f_hypsf>>DA0>>DA1>>DA2>>DA3>>vx>>vy>>vz>>etas;

            gamma=1.0/sqrt(max(1.0e-8,1.0-vx*vx-vy*vy-vz*vz));

            for(int k=0; k!=Neta; k++){
                eta=(k-Neta/2)*deta;                //eta from -5 to 5
                for(int i=0; i!=Npt; i++){
                    pt=ptlo+(pthi-ptlo)*pgauss[i];
                    mt=sqrt(pt*pt+m_pi*m_pi);
                    ptce=pt*cosh(eta);
                    dYdEta=ptce/sqrt(m_pi*m_pi+ptce*ptce);//   P/E  
                    Y=atanh(dYdEta*tanh(eta));//  dY/deta=P/E

                    mtcy=mt*cosh(Y-etas);
                    mtsy=mt*sinh(Y-etas);
                    for(int j=0; j!=Nphi; j++){
                        phi=philo + j*dphi;
                        ptcp=pt*cos(phi);
                        ptsp=pt*sin(phi);
                        efkt=gamma*(mtcy-ptcp*vx-ptsp*vy-mtsy*vz)/t_dec;
                        volum=mtcy*DA0-ptcp*DA1-ptsp*DA2-mtsy*DA3;
                        pvf=volum/(exp(efkt)-1.0);   //Gev.fm^3
                        dNdEdPtdPhi[k][i][j] += dYdEta*pvf;
                    }//end for phi
                }//end for pt
            }//end for Eta
        }//end while

        for(int k=0; k!=Neta; k++){
            for(int j=0; j!=Nphi; j++){
                for(int i=0; i!=Npt; i++){
                    pt=ptlo+(pthi-ptlo)*pgauss[i];
                    dNdEdPhi[k][j] += Gevfm3*wgauss[i]*pt*dNdEdPtdPhi[k][i][j];
                }
                f_dNdEdPhi<<" "<<dof*dNdEdPhi[k][j];
                dNdEta[k] += dof*dNdEdPhi[k][j]*dphi;
            }
            f_dNdEdPhi<<"\n";
            eta=(k-Neta/2)*deta;                
            f_dNdEta<<eta<<" "<<dNdEta[k]<<endl;
        }
    }//End Freeze out 4 ( 1<pT<3 )

    void CHydro3D::FrzIsoTime()
    {
        char FileOut[60];
        sprintf(FileOut,"results/dNdEtadPhi%d.dat",IEVENT);
        char FdNdEta[60];
        sprintf(FdNdEta,"results/dNdEta%d.dat",IEVENT);


        ofstream f_dNdEdPhi(FileOut);	
        ofstream f_dNdEta(FdNdEta);	

        CD t_dec=0.130; //freeze out temperature
        CD m_pi=0.140; //pion mass Gev
        //CD m_pi=0.0;
        CD dof=2.0/pow(2.0*PhyConst::pi,3.0);
        double T;

        CD etalo=-10.0;
        CD etahi=10.0;
        CD ptlo =1.0;
        CD pthi =3.0;
        CD philo=0.0;
        CD phihi=2.0*PhyConst::pi;
        CD Gevfm3=pow(1.0/PhyConst::hbarc,3.0);

        CI Neta=41;
        CI Npt=16;
        CI Nphi=41;

        double invslope=0.25;

        double *pgauss, *wgauss;
        double *pgauss1, *wgauss1;
        double *pgaula, *wgaula;
        pgauss = gaulep48; wgauss=gaulew48;
        pgaula = gala15x;  wgaula=gala15w;
        pgauss1 = gaulep16; wgauss1=gaulew16;

        double deta, dpt, dphi;
        double mt, pt, phi, Y, eta, gamma;
        double DA0, DA1, DA2, DA3, vx, vy, vz, time, x, y, etas;
        double dNdEdPtdPhi[Neta][Npt][Nphi]={0.0};
        double dNdEdPhi[Neta][Nphi]={0.0};
        deta = (etahi-etalo)/(Neta-1);
        dphi = (phihi-philo)/(Nphi-1);



        for(int IX=NX0+2; IX<=NX-2; IX++)
            for(int IY=NY0+2; IY<=NY-2; IY++)
                for(int IZ=NZ0+2; IZ<=NZ-2; IZ++){
                    vx = Vx(IX, IY, IZ);
                    vy = Vy(IX, IY, IZ);
                    vz = Vz(IX, IY, IZ);
                    etas=IZ*dz;
                    T = Eos.T( Ed(IX,IY,IZ),Bd(IX,IY,IZ) );

                    gamma=Gamma(IX,IY,IZ);

                    for(int k=0; k!=Neta; k++){
                        eta=(k-Neta/2)*deta;                //eta from -5 to 5
                        for(int i=0; i!=Npt; i++){
                            pt=(pthi-ptlo)*pgauss1[i]+ptlo;
                            mt=sqrt(pt*pt+m_pi*m_pi);
                            double ptce=pt*cosh(eta);
                            double dYdEta=ptce/sqrt(m_pi*m_pi+ptce*ptce);//   P/E  
                            Y=atanh(dYdEta*tanh(eta));//  dY/deta=P/E
                            double mtcy=mt*cosh(Y-etas);
                            double mtsy=mt*sinh(Y-etas);
                            for(int j=0; j!=Nphi; j++){
                                phi=philo + j*dphi;
                                double ptcp=pt*cos(phi);
                                double ptsp=pt*sin(phi);
                                double efkt=gamma*max(mtcy-ptcp*vx-ptsp*vy-mtsy*tau*vz,PhyConst::acu)\
                                            /max(T,PhyConst::acu);
                                DA0 = tau*dx*dy*dz;
                                double volum=mtcy*DA0;
                                //double pvf=volum/max(exp(efkt)-1.0,PhyConst::acu);
                                double pvf=volum*exp(-efkt);
                                dNdEdPtdPhi[k][i][j] += dYdEta*pvf;
                            }//end for phi
                        }//end for pt
                    }//end for Eta
                }


        double dNdEta[Neta]={0.0};
        for(int k=0; k!=Neta; k++){
            for(int j=0; j!=Nphi; j++){
                for(int i=0; i!=Npt; i++){
                    pt=ptlo+(pthi-ptlo)*pgauss1[i];
                    dNdEdPhi[k][j] += Gevfm3*pt*wgauss1[i]*dNdEdPtdPhi[k][i][j];
                }
                dNdEta[k] += dof*dNdEdPhi[k][j]*dphi;
                f_dNdEdPhi<<" "<<dof*dNdEdPhi[k][j];
            }
            f_dNdEdPhi<<"\n";
            eta=(k-Neta/2)*deta;                
            f_dNdEta<<eta<<" "<<dNdEta[k]<<endl;
        }
    }
    ////////////////////////////////////////////////////////////////////////////////
    //
    //Qin: Evolution profiles for jet calculation




///////////////////////////////////////////////////////
/* Yingru:  Data output procedure  */
///////////////////////////////////////////////////////

void CHydro3D::Output(const int & NT)
{
     char HOME[]=".";

        char fbuf[256];
        sprintf(fbuf, "%s/results/event%d/CCNU_hydro.dat", HOME, IEVENT);
        ofstream f_data(fbuf, ios::app);

        char fbuf2[256];
        sprintf(fbuf2, "%s/results/event%d/Anisotropy.dat", HOME, IEVENT);
        ofstream f_ecc(fbuf2, ios::app);

        char fbuf3[256];
        sprintf(fbuf3, "%s/results/event%d/Hydro_steps.dat", HOME, IEVENT);
        ofstream f_time(fbuf3, ios::app);
        f_time <<tau<<"\n";
        f_time.close();
        
        
        cout<<"time="<<tau<<"\n";
        cout<<"After evolution Ed="<<"\n";
        for (int I=0; I<=NX; I+=10)
        {
            for (int J=0; J<=NY; J+=10)
            {
                cout<<setw(12)<<Ed(I,J,0);
            }
            cout<<"\n";
        }

        f_data<<"#tau="<<tau<<"\n"; 
        f_data<<"#           x          y    rapidity          Vx          Vy          Vz    Edensity        Temp         TJT        TJAT"<<"\n";

        for (int i=NX0; i<=NX; i++)
        {
            double eps_const= 0.0;
            double r=1;
            for (int j=NY0; j<=NY; j++)
            {
                for (int k=NZ0; k<=NZ; k++)
                {
                    f_data<<setprecision(5)<<setw(12)<<dx*i<<setw(12)<<dy*j<<setw(12)<<dz*k<<setw(12)<<Vx(i,j,k)<<setw(12)<<Vy(i,j,k)<<setw(12)<<Vz(i,j,k)*tau<<setw(12)<<Ed(i,j,k)<<setw(12)<<Eos.T(Ed(i,j,k),Bd(i,j,k))<<setw(12)<<TJT(i,j,k)<<setw(12)<<TJAT(i,j,k)<<"\n";
                }
            }
        }

        double epsx, epsp;
        Accecentricity(epsx, epsp);
        f_ecc<<"#tau="<<tau<<"\n";
        f_ecc<<setprecision(5)<<setw(12)<<epsx<<setw(12)<<epsp<<"\n";


        f_data.close();
        f_ecc.close();

	
}//End Output()
////////////////////////////////////////////Yingru
    ////////////////////////////////////////////////////////////////////////////////

// Qin: Change initial condition here
    /* Hydrodynamic Initialization *********/
    void CHydro3D::Initialization(char* hydro_config )
    {
        char HOME[]=".";
        char buff[200];       //Store comment
        ifstream fconfig(hydro_config,ios_base::in);
        if(!fconfig.is_open()){
            cout<<"Error openning hydro configure file!\n";
        }
        else{
            fconfig>>setting.Ieos;         fconfig.getline(buff,200);
            fconfig>>setting.b;            fconfig.getline(buff,200);
            fconfig>>setting.t0;           fconfig.getline(buff,200);
            fconfig>>setting.e0;           fconfig.getline(buff,200);
            fconfig>>setting.edec;         fconfig.getline(buff,200);
            fconfig>>setting.Ieos2dec;     fconfig.getline(buff,200);

            fconfig.getline(buff,200);     //Line separator
            fconfig.getline(buff,200);     //Line separator


            fconfig>>setting.Iglaubercgc;  fconfig.getline(buff,200);
            fconfig>>setting.Iein;         fconfig.getline(buff,200);
            fconfig>>setting.Hwn;          fconfig.getline(buff,200);
            fconfig>>setting.A;            fconfig.getline(buff,200);         
            fconfig>>setting.Si0;          fconfig.getline(buff,200);       	
            fconfig>>setting.Ro0;          fconfig.getline(buff,200);           
            fconfig>>setting.Eta;          fconfig.getline(buff,200);               
            fconfig>>setting.Ra;           fconfig.getline(buff,200);                     
            fconfig>>setting.Eta_flat;     fconfig.getline(buff,200);                     
            fconfig>>setting.Eta_gw;       fconfig.getline(buff,200);                     

            fconfig.getline(buff,200);     //Line separator
            fconfig.getline(buff,200);     //Line separator

            fconfig>>setting.ViscousC;     fconfig.getline(buff,200);               
            fconfig>>setting.VisBeta;      fconfig.getline(buff,200);                
            fconfig>>setting.VisBulk;      fconfig.getline(buff,200);          
            fconfig>>setting.IRelaxBulk;   fconfig.getline(buff,200);       	
            fconfig>>setting.BulkTau;      fconfig.getline(buff,200);         

            fconfig.getline(buff,200);     //Line separator
            fconfig.getline(buff,200);     //Line separator

            fconfig>>setting.NXD;          fconfig.getline(buff,200);           
            fconfig>>setting.NYD;          fconfig.getline(buff,200);       
            fconfig>>setting.NZD;          fconfig.getline(buff,200);       
            fconfig>>setting.NTD;          fconfig.getline(buff,200);          
            fconfig>>setting.dt;           fconfig.getline(buff,200);
            fconfig>>setting.dx;           fconfig.getline(buff,200);
            fconfig>>setting.dy;           fconfig.getline(buff,200);
            fconfig>>setting.deta;         fconfig.getline(buff,200);
        }
        fconfig.close();

        //if(IEVENT==5){setting.e0=40.0; setting.t0=0.511795;}
        //else if(IEVENT==6){setting.e0=35.0; setting.t0=0.567756;}
        //else if(IEVENT==7){setting.e0=30.0; setting.t0=0.640234;}
        //else if(IEVENT==8){setting.dt=0.02; setting.dx=0.2; setting.dy=0.2; setting.deta=0.2;}
        //else if(IEVENT==9){setting.NTD=1; setting.NXD=1; setting.NYD=1; setting.NZD=1;}
        //else if(IEVENT==10){setting.NTD=5; setting.NXD=2; setting.NYD=2; setting.NZD=2;}
        //else if(IEVENT==11){setting.NTD=20; setting.NXD=4; setting.NYD=4; setting.NZD=2;}

        /* b got from YiFei */
//        if(IEVENT==0){setting.b=3.156;}
//        else if(IEVENT==1){setting.b=5.66;}
//        else if(IEVENT==2){setting.b=7.33;}
//        else if(IEVENT==3){setting.b=8.68;}
//        else if(IEVENT==4){setting.b=4.5;}
//        else if(IEVENT==5){setting.b=6.3;}
//        else if(IEVENT==6){setting.b=7.5;}
//        else if(IEVENT==7){setting.b=2.4;}

        Eos.SetEos(setting.Ieos);

        CNucleus A1(setting.A, setting.Ro0, setting.Ra, setting.Eta);      //projectile nucleu
        CNucleus A2(setting.A, setting.Ro0, setting.Ra, setting.Eta);      //target nucl
        CAAcollision AA(A1, A2, setting.Si0, setting.b, setting.Hwn, 1.0-setting.Hwn, setting.Eta_flat, setting.Eta_gw);                 //collision
        CAAcollision AACentral(A1, A2, setting.Si0, 0.0, setting.Hwn, 1.0-setting.Hwn, setting.Eta_flat, setting.Eta_gw);                //central collisions for kFactor calculation

        double kFactor=setting.e0/AACentral.NCollision(0.0, 0.0, 0.0);
        t0=setting.t0;
        dt=setting.dt;
        dx=setting.dx;
        dy=setting.dy;
        dz=setting.deta;
        tau=t0;
        ED_MAX=0.0;

        //PPInitial3();


// Qin: Change initial condition here (z is eta_s)
// Qin: Nx=152, Ny=152, Neta=42 see main.cpp for details
// dx = 0.1fm, dy = 0.1fm, deta = 0.3 see Vishydro1.inp for details
// Qin: Nx=77, Ny=77, Neta=42 see main.cpp for details
// dx = 0.2fm, dy = 0.2fm, deta = 0.3 see Vishydro1.inp for details
// Qin: Nx=51, Ny=51, Neta=31 see main.cpp for details
// dx = 0.3fm, dy = 0.3fm, deta = 0.4 see Vishydro1.inp for details


/*
             ofstream fglb("ED_Glauber.dat");
               for(int i=NX0; i<=NX; i++)
               for(int j=NY0; j<=NY; j++)
               for(int k=NZ0; k<=NZ; k++){
               double x= i*dx;
               double y= j*dy;
               double z= k*dz;
               Vx(i,j,k)= 0.0;
               Vy(i,j,k)= 0.0;
               Vz(i,j,k)= 0.0;


        //if(abs(z)<=5.0){
        Ed(i,j,k) = kFactor*AA.NCollision(x,y,z);  
        if(ED_MAX<Ed(i,j,k))ED_MAX=Ed(i,j,k);
        //Bd(i,j,k)= 0.08*Ed(i,j,k);
        Bd(i,j,k)= 0.0;
        Pr(i,j,k)=Eos.P(Ed(i,j,k),Bd(i,j,k));
        //}
//        else{
//        	Ed(i,j,k)=0.0;
//        	Bd(i,j,k)=0.0;
//        	Pr(i,j,k)=0.0;
//        }

//	write energy density to data file

        if(Ed(i,j,k)>PhyConst::acu)
        fglb<<i<<' '<<j<<' '<<k<<' '<<Ed(i,j,k)<<endl;



        double U0=Gamma(i,j,k);
        double U1=U0*Vx(i,j,k);
        double U2=U0*Vy(i,j,k);
        double U3=U0*Vz(i,j,k);

        double epp=Ed(i,j,k)+Pr(i,j,k);
        TT00(i,j,k)=t0*(epp*U0*U0-Pr(i,j,k));
        TT01(i,j,k)=t0*(epp*U0*U1);
        TT02(i,j,k)=t0*(epp*U0*U2);
        TT03(i,j,k)=t0*(epp*U0*U3);
        TJT(i,j,k)=t0*Bd(i,j,k)*U0;
        TJAT(i,j,k)=t0*Bd(i,j,k)*U0;
        }
        fglb.close();
*/

//	read my hydro input (Qin)

	bool iread;
    	double E0, B0, T00, T01, T02, T03, JT0, JA0, EPV, M; // Root finding parameters
	CD acu = PhyConst::acu;

        ifstream finput("dt_hydro_input.dat",ios_base::in);
        if(!finput.is_open()){
            cout<<"Error openning initial condition file!\n";
            exit (EXIT_FAILURE);
        }

        for(int I=NX0; I<=NX; I++)
        for(int J=NY0; J<=NY; J++)
        for(int K=NZ0; K<=NZ; K++)
	{
        	double x = I*dx;
        	double y = J*dy;
        	double z = K*dz;

		//cout << "Here1: " << I << setw(16) << J << setw(16) << K << endl;
		//cout << "Here2: " << dx << setw(16) << dy << setw(16) << dz << endl;
		//cout << "Here3: " << x << setw(16) << y << setw(16) << z << endl;

//		read I, J, K, TT00(I,J,K), TT01(I,J,K), TT02(I,J,K), TT03(I,J,K)
//		iread = finput >> I >> J >> K >> TT00(I, J, K) >> TT01(I,J,K) >> TT02(I,J,K) >> TT03(I,J,K) >> TJT(I,J,K);
//		read I, J, K, TT00(I,J,K), TT01(I,J,K), TT02(I,J,K), TT03(I,J,K)
		iread = finput >> I >> J >> K >> T00 >> T01 >> T02 >> T03 >> JT0;

		if (!iread) break;

                TT00(I,J,K) = T00 * tau;
                TT01(I,J,K) = T01 * tau;
                TT02(I,J,K) = T02 * tau;
                TT03(I,J,K) = T03 * tau;
                TJT(I,J,K) = JT0 * tau;

		//if (I == 0 && J == 0)
		if (I == 0 && J == 0 && K == 0)
		{
		cout << "Here4: " << I << setw(16) << J << setw(16) << K << endl;
		cout << "Here5: " << T00 << setw(16) << T01 << setw(16) << T02 << setw(16) << T03 << setw(16) << JT0 << endl;
		}

//		find the local energy density
		T03 = T03 * tau;
                Rootfinding(E0, B0, T00, T01, T02, T03, JT0);

                Ed(I,J,K)=max(0.0, E0);
                if(ED_MAX<Ed(I,J,K))ED_MAX=Ed(I,J,K);

                Bd(I,J,K)=max(0.0, B0);
                Pr(I,J,K)=Eos.P( Ed(I,J,K),Bd(I,J,K) );
                EPV=max( acu, T00+Pr(I,J,K) );

//		get the local flow velocity
                Vx(I,J,K)=T01/EPV;
                Vy(I,J,K)=T02/EPV;
                Vz(I,J,K)=T03/EPV/tau;

                double NORM=(Vx(I,J,K)*Vx(I,J,K)+Vy(I,J,K)*Vy(I,J,K)+tau*tau*
                        Vz(I,J,K)*Vz(I,J,K));
                if(NORM>1.0)
		{
                    NORM = sqrt(NORM);
                    Vx(I,J,K) = Vx(I,J,K)/NORM;
                    Vy(I,J,K) = Vy(I,J,K)/NORM;
                    Vz(I,J,K) = Vz(I,J,K)/NORM;
                }

		//if (I == 0 && J == 0)
		if (I == 0 && J == 0 && K == 0)
		{
		cout << "Here6: " << Ed(I,J,K) << setw(16) << Bd(I,J,K) << setw(16) << Pr(I,J,K) << setw(16) << Vx(I,J,K) << setw(16) << Vy(I,J,K) << setw(16) << Vz(I,J,K) << endl;
		}

	}
	finput.close();

        cout<<"ED_MAX="<<ED_MAX<<endl;

    }//Initialization of energy density , velocity and energy momentum tensor
    ////////////////////////////////////////////////////////////////////////////////

    CHydro3D::CHydro3D(int Xl, int Xh, int Yl, int Yh, int Zl, int Zh)\
        :Ed(Xl,Xh, Yl,Yh, Zl,Zh), \
        Bd(Xl,Xh, Yl,Yh, Zl,Zh), \
        Pr(Xl,Xh, Yl,Yh, Zl,Zh), \
        Vx(Xl,Xh, Yl,Yh, Zl,Zh), \
        Vy(Xl,Xh, Yl,Yh, Zl,Zh), \
        Vz(Xl,Xh, Yl,Yh, Zl,Zh), \
        TT00(Xl,Xh, Yl,Yh, Zl,Zh), \
        TT01(Xl,Xh, Yl,Yh, Zl,Zh), \
        TT02(Xl,Xh, Yl,Yh, Zl,Zh), \
        TT03(Xl,Xh, Yl,Yh, Zl,Zh), \
        TJT(Xl,Xh, Yl,Yh, Zl,Zh),\
        TJAT(Xl,Xh, Yl,Yh, Zl,Zh)
        {
            NX0 = Xl;
            NX  = Xh;
            NY0 = Yl;
            NY  = Yh;
            NZ0 = Zl;
            NZ  = Zh;
        }//Construct an 3D hydro class

    /*Calc Ed and Bd through T^\mu\nu */
    void CHydro3D::Rootfinding(double &E0, double &B0, CD &T00IJ, CD &T01IJ, CD &T02IJ, CD &T03IJ, CD &JTIJ)
    {
        CI IMax=500000;
        CD acu=PhyConst::acu;
        int i=0;
        double Pr, E1, B1;
        E1=T00IJ;
        B1=JTIJ;
        double K= T01IJ*T01IJ +T02IJ*T02IJ +T03IJ*T03IJ;
        do{
            E0=E1;
            B0=B1;
            Pr=Eos.P(E0,B0);
            E1=T00IJ-K/(T00IJ+Pr);
            B1=JTIJ*sqrt(max(0.0,(E0+Pr)/(T00IJ+Pr)));
            i++;
            //assert(i<=IMax);
        }while((i<IMax)&&(abs(E1-E0)>acu));
        //}while((i<IMax)&&(abs(E1-E0)>acu)||(abs(B1-B0)>acu));
        if(i>=IMax){
            cout<<"T00IJ="<<T00IJ<<" T01IJ="<<T01IJ<<" T02IJ="<<T02IJ<<" T03IJ="<<T03IJ<<endl;
            cout<<"E1="<<E1<<", E0="<<E0<<endl;
            cout<<"B1="<<B1<<", B0="<<B0<<endl;
            //assert(i<IMax);
        }

}//Solve Ed, Bd through iteration with an initial value of T00 and JT.

void CHydro3D::StepUpdate(){
    const double acu=PhyConst::acu;
    M3D ScT0(NX0,NX, NY0,NY, NZ0,NZ);     //source terms for T^{\tau\tau}
    M3D ScT1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScT2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScT3(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScJT(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScJA(NX0,NX, NY0,NY, NZ0,NZ);
    ED_MAX = 0.0;

    double E0, B0, T00, T01, T02, T03, JT0,JA0, EPV, M; // Root finding parameters
    for(int I=NX0+2; I<=NX-2; I++)
        for(int J=NY0+2; J<=NY-2; J++)
            for(int K=NZ0+2; K<=NZ-2; K++){
                T00=TT00(I,J,K)/tau;
                T01=TT01(I,J,K)/tau;
                T02=TT02(I,J,K)/tau;
                T03=TT03(I,J,K);      //T03=(e+P)g^2*v_{etas}*tau
                JT0= TJT(I,J,K)/tau;
                JA0= TJAT(I,J,K)/tau;

                T00=max(0.0,T00);
                T01=(abs(T01)>acu)?T01:0.0;
                T02=(abs(T02)>acu)?T02:0.0;
                T03=(abs(T03)>acu)?T03:0.0;
                JT0=(abs(JT0)>acu)?JT0:0.0;
                JA0=(abs(JA0)>acu)?JA0:0.0;
                M=sqrt(T01*T01+T02*T02+T03*T03);
                T00=(M>T00) ? (1.000001*M) : T00;

                Rootfinding(E0, B0, T00, T01, T02, T03, JT0);
                Ed(I,J,K)=max(0.0, E0);//here acu changed to 0.0 4/27
                //Ed(I,J,K)=max(acu, E0);  //here 0.0 changed to acu 1/27
                if(ED_MAX<Ed(I,J,K))ED_MAX=Ed(I,J,K);

                Bd(I,J,K)=max(0.0, B0);
                Pr(I,J,K)=Eos.P( Ed(I,J,K),Bd(I,J,K) );
                EPV=max( acu, T00+Pr(I,J,K) );
                Vx(I,J,K)=T01/EPV;
                Vy(I,J,K)=T02/EPV;
                Vz(I,J,K)=T03/EPV/tau;

                double NORM=(Vx(I,J,K)*Vx(I,J,K)+Vy(I,J,K)*Vy(I,J,K)+tau*tau*
                        Vz(I,J,K)*Vz(I,J,K));
                if(NORM>1.0){//Add on 02/22/2011
                    NORM = sqrt(NORM);
                    Vx(I,J,K) = Vx(I,J,K)/NORM;
                    Vy(I,J,K) = Vy(I,J,K)/NORM;
                    Vz(I,J,K) = Vz(I,J,K)/NORM;
                }

                TT00(I,J,K)=T00*tau;
                TT01(I,J,K)=T01*tau;
                TT02(I,J,K)=T02*tau;
                TT03(I,J,K)=T03;
                TJT (I,J,K)=JT0*tau;
                //TJAT(I,J,K)=JA0*tau;
            }


    SetBoundary(Ed);
    SetBoundary(Bd);
    SetBoundary(Pr);
    SetBoundary(Vx);
    SetBoundary(Vy);
    SetBoundary(Vz);
    VzSmooth();

    //Moved from blow CalcSources
    //ATsplit(TT00, TT01, TT02, TT03, TJT, TJAT, Ed, Bd, Vx, Vy, Vz, dt, dx, dy, dz, 0.125);
    QuickTsplit(TT00, TT01, TT02, TT03, TJT, Ed, Bd, Vx, Vy, Vz, dt, dx, dy, dz, 0.125);
    //Tsplit(TT00, TT01, TT02, TT03, TJT, Ed, Bd, Vx, Vy, Vz, dt, dx, dy, dz, 0.125);

    CalcSources(ScT0,ScT1,ScT2,ScT3,ScJT,ScJA);
    //AnomalySources(ScT0,ScT1,ScT2,ScT3,ScJT,ScJA);


    for(int I=NX0+2; I<=NX-2; I++)
        for(int J=NY0+2; J<=NY-2; J++)
            for(int K=NZ0+2; K<=NZ-2; K++){
                TT00(I,J,K)=TT00(I,J,K)-dt*ScT0(I,J,K);
                TT01(I,J,K)=TT01(I,J,K)-dt*ScT1(I,J,K);
                TT02(I,J,K)=TT02(I,J,K)-dt*ScT2(I,J,K);
                TT03(I,J,K)=TT03(I,J,K)-dt*ScT3(I,J,K);
                TJT(I,J,K) =TJT(I,J,K) -dt*ScJT(I,J,K);
                //TJAT(I,J,K)=TJAT(I,J,K)-dt*ScJA(I,J,K);
            }

    SetBoundary(TT00);
    SetBoundary(TT01);
    SetBoundary(TT02);
    SetBoundary(TT03);
    SetBoundary(TJT);
    //SetBoundary(TJAT);
}


void CHydro3D::StepUpdate2(M3D & H00, M3D & H01, M3D & H02, M3D & H03, M3D & HJT,
        M3D & F00, M3D & F01, M3D & F02, M3D & F03, M3D & FJT,
        CD & DT){
    const double acu=PhyConst::acu;
    //const double acu=PhyConst::acu;
    M3D ScT0(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScT1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScT2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScT3(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScJT(NX0,NX, NY0,NY, NZ0,NZ);
    M3D ScJA(NX0,NX, NY0,NY, NZ0,NZ);
    ED_MAX = 0.0;

    CalcSources2(ScT0,ScT1,ScT2,ScT3,ScJT,ScJA,F00,F03,DT);

    QuickTsplit(H00, H01, H02, H03, HJT, Ed, Bd, Vx, Vy, Vz, DT, dx, dy, dz, 0.125);
    //Tsplit(H00, H01, H02, H03, HJT, Ed, Bd, Vx, Vy, Vz, DT, dx, dy, dz, 0.125);

    for(int I=NX0+2; I<=NX-2; I++)
        for(int J=NY0+2; J<=NY-2; J++)
            for(int K=NZ0+2; K<=NZ-2; K++){
                H00(I,J,K)=H00(I,J,K)-DT*ScT0(I,J,K);
                H01(I,J,K)=H01(I,J,K)-DT*ScT1(I,J,K);
                H02(I,J,K)=H02(I,J,K)-DT*ScT2(I,J,K);
                H03(I,J,K)=H03(I,J,K)-DT*ScT3(I,J,K);
                HJT(I,J,K)=HJT(I,J,K)-DT*ScJT(I,J,K);
            }

    double time=tau+DT;
    double E0, B0, T00, T01, T02, T03, JT0,JA0, EPV, M; // Root finding parameters
    for(int I=NX0+2; I<=NX-2; I++)
        for(int J=NY0+2; J<=NY-2; J++)
            for(int K=NZ0+2; K<=NZ-2; K++){
                T00=H00(I,J,K)/time;
                T01=H01(I,J,K)/time;
                T02=H02(I,J,K)/time;
                T03=H03(I,J,K);      //T03=(e+P)g^2*time*v_{eta}
                JT0=HJT(I,J,K)/time;

                T00=max(0.0,T00);
                T01=(abs(T01)>acu)?T01:0.0;
                T02=(abs(T02)>acu)?T02:0.0;
                T03=(abs(T03)>acu)?T03:0.0;
                JT0=(abs(JT0)>acu)?JT0:0.0;
                M=sqrt(T01*T01+T02*T02+T03*T03);
                T00=(M>T00) ? (1.000001*M) : T00;
                //T00=(M>T00) ? (1.0+acu)*M : T00;

                Rootfinding(E0, B0, T00, T01, T02, T03, JT0);
                Ed(I,J,K)=max(0.0, E0);//here acu changed to 0.0 4/27

                if(ED_MAX<Ed(I,J,K))ED_MAX=Ed(I,J,K);

                Bd(I,J,K)=max(0.0, B0);
                Pr(I,J,K)=Eos.P( Ed(I,J,K),Bd(I,J,K) );
                Pr(I,J,K)=max(0.0,Pr(I,J,K));

                EPV=max( acu, T00+Pr(I,J,K) );
                Vx(I,J,K)=T01/EPV;
                Vy(I,J,K)=T02/EPV;
                Vz(I,J,K)=T03/EPV/time;


                double U0=1.0/sqrt(max(1.0e-30,1.0-Vx(I,J,K)*Vx(I,J,K)-Vy(I,J,K)*Vy(I,J,K)-time*time*Vz(I,J,K)*Vz(I,J,K)));
                double U1=U0*Vx(I,J,K);
                double U2=U0*Vy(I,J,K);
                double U3=U0*Vz(I,J,K);

                double epp=Ed(I,J,K)+Pr(I,J,K);


                H00(I,J,K)=time*(epp*U0*U0-Pr(I,J,K));
                H01(I,J,K)=time*(epp*U0*U1);
                H02(I,J,K)=time*(epp*U0*U2);
                H03(I,J,K)=time*(epp*U0*U3);
                HJT (I,J,K)=time*Bd(I,J,K);

                //	H00(I,J,K)=T00*time;
                //	H01(I,J,K)=T01*time;
                //	H02(I,J,K)=T02*time;
                //	H03(I,J,K)=T03;
                //	HJT (I,J,K)=JT0*time;
            }


    SetBoundary(Ed);
    SetBoundary(Bd);
    SetBoundary(Pr);
    SetBoundary(Vx);
    SetBoundary(Vy);
    SetBoundary(Vz);
    //VzSmooth();

    SetBoundary(H00);
    SetBoundary(H01);
    SetBoundary(H02);
    SetBoundary(H03);
    SetBoundary(HJT);
}

double minmod(CD & x1, CD &x2, CD &x3);

void CHydro3D::CalcSources2(M3D &ScT0, M3D &ScT1, M3D &ScT2, M3D &ScT3, M3D &ScJT, M3D &ScJA, M3D &F00, M3D &F03, CD &DT)
{/* Calculate the source contribution by using second order Runge Kuta method */
    enum mode{minm, central};
    double PX, PY, PZ, PVX, PVY, PVZ;
    CD s0=1.1;
    const int m=minm;  //source cal method 

    for(int I=NX0+2; I<=NX-2; I++)
        for(int J=NY0+2; J<=NY-2; J++)
            for(int K=NZ0+2; K<=NZ-2; K++){
                if(m==minm){
                    PX=minmod(s0*(Pr(I+1,J,K)-Pr(I,J,K))/dx, (Pr(I+1,J,K)-Pr(I-1,J,K))/(2.0*dx), s0*(Pr(I,J,K)-Pr(I-1,J,K))/dx );
                    PY=minmod(s0*(Pr(I,J+1,K)-Pr(I,J,K))/dy, (Pr(I,J+1,K)-Pr(I,J-1,K))/(2.0*dy), s0*(Pr(I,J,K)-Pr(I,J-1,K))/dy );
                    PZ=minmod(s0*(Pr(I,J,K+1)-Pr(I,J,K))/dz, (Pr(I,J,K+1)-Pr(I,J,K-1))/(2.0*dz), s0*(Pr(I,J,K)-Pr(I,J,K-1))/dz );

                    PVX=minmod(s0*(Pr(I+1,J,K)*Vx(I+1,J,K)-Pr(I,J,K)*Vx(I,J,K))/dx,\
                            (Pr(I+1,J,K)*Vx(I+1,J,K)-Pr(I-1,J,K)*Vx(I-1,J,K))/(2.0*dx),\
                            s0*(Pr(I,J,K)*Vx(I,J,K)-Pr(I-1,J,K)*Vx(I-1,J,K))/dx );

                    PVY=minmod(s0*(Pr(I,J+1,K)*Vy(I,J+1,K)-Pr(I,J,K)*Vy(I,J,K))/dy,\
                            (Pr(I,J+1,K)*Vy(I,J+1,K)-Pr(I,J-1,K)*Vy(I,J-1,K))/(2.0*dy),\
                            s0*(Pr(I,J,K)*Vy(I,J,K)-Pr(I,J-1,K)*Vy(I,J-1,K))/dy );

                    PVZ=minmod(s0*(Pr(I,J,K+1)*Vz(I,J,K+1)-Pr(I,J,K)*Vz(I,J,K))/dz,\
                            (Pr(I,J,K+1)*Vz(I,J,K+1)-Pr(I,J,K-1)*Vz(I,J,K-1))/(2.0*dz),\
                            s0*(Pr(I,J,K)*Vz(I,J,K)-Pr(I,J,K-1)*Vz(I,J,K-1))/dz );
                }
                else if(m==central){
                    PX=(Pr(I+1,J,K)-Pr(I-1,J,K))/(2.0*dx);
                    PY=(Pr(I,J+1,K)-Pr(I,J-1,K))/(2.0*dy);
                    PZ=(Pr(I,J,K+1)-Pr(I,J,K-1))/(2.0*dz);
                    PVX=(Pr(I+1,J,K)*Vx(I+1,J,K)-Pr(I-1,J,K)*Vx(I-1,J,K))/(2.0*dx);
                    PVY=(Pr(I,J+1,K)*Vy(I,J+1,K)-Pr(I,J-1,K)*Vy(I,J-1,K))/(2.0*dy);
                    PVZ=(Pr(I,J,K+1)*Vz(I,J,K+1)-Pr(I,J,K-1)*Vz(I,J,K-1))/(2.0*dz);
                }

                double time=tau+DT-0.5*dt;
                ScT0(I,J,K)=Pr(I,J,K)+pow(Vz(I,J,K),2.0)*time*(F00(I,J,K)+Pr(I,J,K)*time)+time*(PVX+PVY+PVZ);
                ScT1(I,J,K)=time*PX;
                ScT2(I,J,K)=time*PY;
                ScT3(I,J,K)=1.0/time*PZ +2.0/time*F03(I,J,K);
                ScJT(I,J,K)=0.0;
                //ScJA(I,J,K)=0.0;
            }
}
////////Calc the source terms through minmod method/////////////////////////////////
void CHydro3D::CalcSources(M3D &ScT0, M3D &ScT1, M3D &ScT2, M3D &ScT3, M3D &ScJT, M3D &ScJA)
{
    enum mode{minm, central};
    double PX, PY, PZ, PVX, PVY, PVZ;
    CD s0=1.1;
    const int m=minm;  //source cal method 

    for(int I=NX0+2; I<=NX-2; I++)
        for(int J=NY0+2; J<=NY-2; J++)
            for(int K=NZ0+2; K<=NZ-2; K++){
                if(m==minm){
                    PX=minmod(s0*(Pr(I+1,J,K)-Pr(I,J,K))/dx, (Pr(I+1,J,K)-Pr(I-1,J,K))/(2.0*dx), s0*(Pr(I,J,K)-Pr(I-1,J,K))/dx );
                    PY=minmod(s0*(Pr(I,J+1,K)-Pr(I,J,K))/dy, (Pr(I,J+1,K)-Pr(I,J-1,K))/(2.0*dy), s0*(Pr(I,J,K)-Pr(I,J-1,K))/dy );
                    PZ=minmod(s0*(Pr(I,J,K+1)-Pr(I,J,K))/dz, (Pr(I,J,K+1)-Pr(I,J,K-1))/(2.0*dz), s0*(Pr(I,J,K)-Pr(I,J,K-1))/dz );

                    PVX=minmod(s0*(Pr(I+1,J,K)*Vx(I+1,J,K)-Pr(I,J,K)*Vx(I,J,K))/dx,\
                            (Pr(I+1,J,K)*Vx(I+1,J,K)-Pr(I-1,J,K)*Vx(I-1,J,K))/(2.0*dx),\
                            s0*(Pr(I,J,K)*Vx(I,J,K)-Pr(I-1,J,K)*Vx(I-1,J,K))/dx );

                    PVY=minmod(s0*(Pr(I,J+1,K)*Vy(I,J+1,K)-Pr(I,J,K)*Vy(I,J,K))/dy,\
                            (Pr(I,J+1,K)*Vy(I,J+1,K)-Pr(I,J-1,K)*Vy(I,J-1,K))/(2.0*dy),\
                            s0*(Pr(I,J,K)*Vy(I,J,K)-Pr(I,J-1,K)*Vy(I,J-1,K))/dy );

                    PVZ=minmod(s0*(Pr(I,J,K+1)*Vz(I,J,K+1)-Pr(I,J,K)*Vz(I,J,K))/dz,\
                            (Pr(I,J,K+1)*Vz(I,J,K+1)-Pr(I,J,K-1)*Vz(I,J,K-1))/(2.0*dz),\
                            s0*(Pr(I,J,K)*Vz(I,J,K)-Pr(I,J,K-1)*Vz(I,J,K-1))/dz );
                }
                else if(m==central){
                    PX=(Pr(I+1,J,K)-Pr(I-1,J,K))/(2.0*dx);
                    PY=(Pr(I,J+1,K)-Pr(I,J-1,K))/(2.0*dy);
                    PZ=(Pr(I,J,K+1)-Pr(I,J,K-1))/(2.0*dz);
                    PVX=(Pr(I+1,J,K)*Vx(I+1,J,K)-Pr(I-1,J,K)*Vx(I-1,J,K))/(2.0*dx);
                    PVY=(Pr(I,J+1,K)*Vy(I,J+1,K)-Pr(I,J-1,K)*Vy(I,J-1,K))/(2.0*dy);
                    PVZ=(Pr(I,J,K+1)*Vz(I,J,K+1)-Pr(I,J,K-1)*Vz(I,J,K-1))/(2.0*dz);
                }

                ScT0(I,J,K)=Pr(I,J,K)+pow(Vz(I,J,K),2.0)*tau*(TT00(I,J,K)+Pr(I,J,K)*tau)+tau*(PVX+PVY+PVZ);
                ScT1(I,J,K)=tau*PX;
                ScT2(I,J,K)=tau*PY;
                ScT3(I,J,K)=1.0/tau*PZ +2.0/tau*TT03(I,J,K);
                ScJT(I,J,K)=0.0;
                //ScJA(I,J,K)=0.0;
            }
}

void CHydro3D::AnomalySources(M3D &ScT0, M3D &ScT1, M3D &ScT2, M3D &ScT3, M3D &ScJT, M3D &ScJA)
{
    for(int I=NX0+2; I<=NX-2; I++)
        for(int J=NY0+2; J<=NY-2; J++)
            for(int K=NZ0+2; K<=NZ-2; K++){
                ScT0(I,J,K) -= 0.0;
                ScT1(I,J,K) -= 0.0;
                ScT2(I,J,K) -= 0.0;
                ScT3(I,J,K) -= 0.0;
                ScJT(I,J,K) -= 0.0;
                ScJA(I,J,K) -= 0.0;
            }
}

inline double minmod(CD & x1, CD &x2, CD &x3){
    double result;
    if(x1>0.0 && x2>0.0 && x3>0.0)result=min(min(x1,x2),x3);
    else if(x1<0.0 && x2<0.0 && x3<0.0)result=max(max(x1,x2),x3);
    else result=0.0;

    return result;
}/*arXiv:1004.1408v1 hep-ph Appendix  */

void CHydro3D::SetBoundary(M3D &M0)
{//For a large enough box, we can set the out boundary is vacum
    for(int I=NX0; I<=NX; I++)
        for(int J=NY0; J<=NY; J++)
            for(int K=0; K<=1; K++){
                assert(NZ0+K<NZ);       //for 2D hydro where NZ0=NZ=0, NZ0+1 will be out of range
                //M0(I,J,NZ0+K)=0.0;
                //M0(I,J, NZ-K)=0.0;
                //M0(I,J,NZ0+K)=2.0*M0(I,J,NZ0+2)-M0(I,J,NZ0+4-K);
                //M0(I,J, NZ-K)=2.0*M0(I,J,NZ -2)-M0(I,J,NZ -4+K);
                M0(I,J,NZ0+K)=M0(I,J,NZ0+2);
                M0(I,J, NZ-K)=M0(I,J,NZ-2);
            }

    for(int J=NY0; J<=NY; J++)
        for(int K=NZ0; K<=NZ; K++)
            for(int I=0; I<=1; I++){
                assert(NX0+I<NX);
                //M0(NX0+I,J,K)=0.0;
                //M0(NX -I,J,K)=0.0;
                //                    M0(NX0+I,J,K)=2.0*M0(NX0+2,J,K)-M0(NX0+4-I,J,K);
                //                    M0(NX -I,J,K)=2.0*M0(NX -2,J,K)-M0(NX -4+I,J,K);
                M0(NX0+I,J,K)=M0(NX0+2,J,K);
                M0(NX -I,J,K)=M0(NX -2,J,K);
            }

    for(int I=NX0; I<=NX; I++)
        for(int K=NZ0; K<=NZ; K++)
            for(int J=0; J<=1; J++){
                assert(NY0+J<NY);
                //M0(I,NY0+J,K)=0.0;
                //M0(I,NY -J,K)=0.0;
                //                    M0(I,NY0+J,K)=2.0*M0(I,NY0+2,K)-M0(I,NY0+4-J,K);
                //                    M0(I,NY -J,K)=2.0*M0(I,NY -2,K)-M0(I,NY -4+J,K);
                M0(I,NY0+J,K)=M0(I,NY0+2,K);
                M0(I,NY -J,K)=M0(I,NY -2,K);
    }
}//extrapolation

inline double CHydro3D::Gamma(CI &i, CI &j, CI &k)
{//calc gamma=1.0/sqrt(1-vx^2-vy^2-vz^2);
    return 1.0/sqrt( max(1.0e-30, 1.0-Vx(i,j,k)*Vx(i,j,k)-Vy(i,j,k)*Vy(i,j,k)-tau*tau*Vz(i,j,k)*Vz(i,j,k) ));
}

void CHydro3D::Accecentricity(double & epsx, double & epsp){
    double EpsX1=0.0;
    double EpsX2=0.0;
    double EpsP1=0.0;
    double EpsP2=0.0;
    double U0, U1, U2, U3;
    double xx, yy, x2, y2, ep;
    double TXX, TYY;

    const int k=0;
    for(int i=0; i<=NX; i++)
        for(int j=0; j<=NY; j++){
            xx=i*dx;
            yy=j*dy;
            x2=xx*xx;
            y2=yy*yy;
            ep=Ed(i,j,k)+Pr(i,j,k);
            U0=Gamma(i,j,k);
            U1=U0*Vx(i,j,k);
            U2=U0*Vy(i,j,k);
            U3=U0*Vz(i,j,k);
            TXX=ep*U1*U1+Pr(i,j,k);
            TYY=ep*U2*U2+Pr(i,j,k);

            EpsX1=EpsX1 + (y2-x2)*Ed(i,j,k)*U0;
            EpsX2=EpsX2 + (y2+x2)*Ed(i,j,k)*U0;

            EpsP1=EpsP1 + TXX-TYY;
            EpsP2=EpsP2 + TXX+TYY;
        }
    epsx=EpsX1/EpsX2;
    epsp=EpsP1/EpsP2;
}//Calc the Accecentricity epsx and epsp

/////////////////////////////////////////////////////////////////
void CHydro3D::CalcHyperSf(M3D &Ed0,  M3D &Vx0, M3D &Vy0, M3D &Vz0) 
{
    double tau0=tau-dt;//previous time
    CD Edec = setting.edec;

    char FHyperSF[60];
    sprintf(FHyperSF,"results/HyperSf%d.dat",IEVENT);
    ofstream f_frz(FHyperSF,ios::app);
    double DA0,DA1,DA2,DA3,v1,v2,v3,time,x,y,z;
    double sign;
    double ratio;


    for(int i=NX0+2; i<=NX-2; i++)
        for(int j=NY0+2; j<=NY-2; j++)
            for(int k=NZ0+2; k<=NZ-2; k++){
                sign = ((Ed0(i,j,k)-Edec)>=0) ? 1.0 : -1.0;
                if( (Ed0(i,j,k)-Edec)*(Ed(i,j,k)-Edec) <0.0 ){
                    ratio=abs(Edec-Ed0(i,j,k))/abs(Ed(i,j,k)-Ed0(i,j,k));
                    time=tau0+ratio*dt;       
                    DA0=sign*time*dx*dy*dz;
                    DA1=0.0;
                    DA2=0.0;
                    DA3=0.0;
                    v1=(1.0-ratio)*Vx0(i,j,k)+ratio*Vx(i,j,k);    
                    v2=(1.0-ratio)*Vy0(i,j,k)+ratio*Vy(i,j,k);    
                    v3=(1.0-ratio)*Vz0(i,j,k)+ratio*Vz(i,j,k);
                    z=k*dz;

                    //f_frz<<DA0 <<' '<<DA1 <<' '<<DA2 <<' '<<DA3;
                    //f_frz<<' '<<v1 <<' '<<v2 <<' '<<v3*time;
                    //f_frz<<' '<<z <<"\n";

                    SF A;
                    A.DA0 = DA0;
                    A.DA1 = DA1;
                    A.DA2 = DA2;
                    A.DA3 = DA3;
                    A.vx  = v1;
                    A.vy  = v2;
                    A.vz  = v3*time;
                    A.etas= z;
                    hypsf.push_back(A);

                }

                if( (Ed0(i+1,j,k)-Edec)*(Ed0(i,j,k)-Edec) <0.0 ){
                    ratio=abs(Edec-Ed0(i,j,k))/abs(Ed0(i+1,j,k)-Ed0(i,j,k));
                    DA0=0.0;
                    DA1=-sign*tau0*dt*dy*dz;
                    DA2=0.0;
                    DA3=0.0;

                    v1=(1.0-ratio)*Vx0(i,j,k)+ratio*Vx0(i+1,j,k);     
                    v2=(1.0-ratio)*Vy0(i,j,k)+ratio*Vy0(i+1,j,k);     
                    v3=(1.0-ratio)*Vz0(i,j,k)+ratio*Vz0(i+1,j,k);     
                    z=k*dz;
                    time=tau0;

                    //						f_frz<<DA0 <<' '<<DA1 <<' '<<DA2 <<' '<<DA3;
                    //						f_frz<<' '<<v1 <<' '<<v2 <<' '<<v3*time;
                    //						f_frz<<' '<<z <<"\n";
                    SF A;
                    A.DA0 = DA0;
                    A.DA1 = DA1;
                    A.DA2 = DA2;
                    A.DA3 = DA3;
                    A.vx  = v1;
                    A.vy  = v2;
                    A.vz  = v3*time;
                    A.etas= z;
                    hypsf.push_back(A);

                }

                if( (Ed0(i,j+1,k)-Edec)*(Ed0(i,j,k)-Edec) <0.0 ){
                    ratio=abs(Edec-Ed0(i,j,k))/abs(Ed0(i,j+1,k)-Ed0(i,j,k));
                    DA0=0.0;
                    DA1=0.0;
                    DA2=-sign*tau0*dt*dx*dz;
                    DA3=0.0;

                    v1=(1.0-ratio)*Vx0(i,j,k)+ratio*Vx0(i,j+1,k);     //    v1=Vx0(i,j,k);     		//
                    v2=(1.0-ratio)*Vy0(i,j,k)+ratio*Vy0(i,j+1,k);     //    v2=Vy0(i,j,k);     		//
                    v3=(1.0-ratio)*Vz0(i,j,k)+ratio*Vz0(i,j+1,k);     //    v3=Vz0(i,j,k);     		//
                    z=k*dz;
                    time=tau0;

                    //						f_frz<<DA0 <<' '<<DA1 <<' '<<DA2 <<' '<<DA3;
                    //						f_frz<<' '<<v1 <<' '<<v2 <<' '<<v3*time;
                    //						f_frz<<' '<<z <<"\n";

                    SF A;
                    A.DA0 = DA0;
                    A.DA1 = DA1;
                    A.DA2 = DA2;
                    A.DA3 = DA3;
                    A.vx  = v1;
                    A.vy  = v2;
                    A.vz  = v3*time;
                    A.etas= z;
                    hypsf.push_back(A);

                }

                if( (Ed0(i,j,k+1)-Edec)*(Ed0(i,j,k)-Edec) <0.0 ){
                    ratio=abs(Edec-Ed0(i,j,k))/abs(Ed0(i,j,k+1)-Ed0(i,j,k));
                    DA0=0.0;
                    DA1=0.0;
                    DA2=0.0;
                    DA3=-sign*dt*dx*dy;  
                    v1=(1.0-ratio)*Vx0(i,j,k)+ratio*Vx0(i,j,k+1);  
                    v2=(1.0-ratio)*Vy0(i,j,k)+ratio*Vy0(i,j,k+1);  
                    v3=(1.0-ratio)*Vz0(i,j,k)+ratio*Vz0(i,j,k+1);  
                    time=tau0;
                    z=(k+ratio)*dz;         

                    //						f_frz<<DA0 <<' '<<DA1 <<' '<<DA2 <<' '<<DA3;
                    //						f_frz<<' '<<v1 <<' '<<v2 <<' '<<v3*time;
                    //						f_frz<<' '<<z <<"\n";
                    SF A;
                    A.DA0 = DA0;
                    A.DA1 = DA1;
                    A.DA2 = DA2;
                    A.DA3 = DA3;
                    A.vx  = v1;
                    A.vy  = v2;
                    A.vz  = v3*time;
                    A.etas= z;
                    hypsf.push_back(A);

                }
            }

}//Calc the hydper surface 
//////////////////////////////////////////////////////////////////

/////////////VzSmooth////////////////////////////////////////////
void CHydro3D::VzSmooth()
{
    const double acu=PhyConst::acu;
    for(int i=Vz.pag_l; i<=Vz.pag_h; i++)
        for(int j=Vz.row_l; j<=Vz.row_h; j++){
            bool Vzfound=true;
            for(int k=Vz.col_l+4; k<=0; k++){
                if(Vzfound && abs(Vz(i,j,k))>acu){
                    Vzfound=false;
                    for(int kn=1; kn<=4; kn++){
                        Vz(i,j,k-kn)=max(-1.0/tau,2.0*Vz(i,j,k)-Vz(i,j,k+kn));
                        Vz(i,j,k-kn)= (Vz(i,j,k-kn)<0.0) ? Vz(i,j,k-kn) : 0.0;
                    }
                }
            }

            Vzfound=true;
            for(int k=Vz.col_h-4; k>=0; k--){
                if(Vzfound && abs(Vz(i,j,k))>acu){
                    Vzfound=false;
                    for(int kn=1; kn<=4; kn++){
                        Vz(i,j,k+kn)=min(1.0/tau,2.0*Vz(i,j,k)-Vz(i,j,k-kn));
                        Vz(i,j,k+kn)= (Vz(i,j,k+kn)>0.0) ? Vz(i,j,k+kn) : 0.0;
                    }
                }
            }
        }
}
///////////////////////////VzSmooth end///////////////////////////


void CHydro3D::CalcHyperSf2(const int &NT, M3D &Ed0,  M3D &Vx0, M3D &Vy0, M3D &Vz0){
    char FHyperSF[60];
    sprintf(FHyperSF,"results/HyperSf%d.dat",IEVENT);

    char FHyperSF_tx[60];
    sprintf(FHyperSF_tx,"results/Tau_x%d.dat",IEVENT);
    char FHyperSF_ty[60];
    sprintf(FHyperSF_ty,"results/Tau_y%d.dat",IEVENT);
    char FHyperSF_th[60];
    sprintf(FHyperSF_th,"results/Tau_h%d.dat",IEVENT);
    ofstream f_frz_taux(FHyperSF_tx,ios::app);
    ofstream f_frz_tauy(FHyperSF_ty,ios::app);
    ofstream f_frz_tauh(FHyperSF_th,ios::app);

    ofstream f_frz(FHyperSF,ios::app);
    double x,y,z, DA0, DA1, DA2, DA3;
    double v1, v2, v3;
    DA0=0.0;
    DA1=0.0;
    DA2=0.0;
    DA3=0.0;
    CD acu=PhyConst::acu;
    CD Edec = setting.edec;
    bool intersection;
    double EDCorners[16];
    double VLCorners[16][3];   //vx, vy, vz for the 16 corners
    double DTD=setting.NTD*dt;
    double DXD=setting.NXD*dx;
    double DYD=setting.NYD*dy;
    double DZD=setting.NZD*dz;


//    cout<<"tau="<<tau<<endl;  //Modified by Yingru: add as a comment
//    cout<<"NT="<<NT<<endl;
    cout<<endl;
    if(ED_MAX<Edec)DTD= (NT%setting.NTD)*dt;//for the last several time steps

    double tau0=tau-DTD;
    Cube cube1(Edec,DTD,DXD,DYD,DZD);    // 3+1 d energy density cube

    if( ((ED_MAX<Edec)||(NT%setting.NTD==0)) ){
        //Do hypersurface calc here
        int DI=setting.NXD;
        int DJ=setting.NYD;
        int DK=setting.NZD;
        for(int i=NX0+2; i<=NX-2-DI; i=i+DI)
            for(int j=NY0+2; j<=NY-2-DJ; j=j+DJ)
                for(int k=NZ0+2; k<=NZ-2-DK; k=k+DK){
                    if(     (Edec-Ed0(i,j,k))*(Edec-Ed(i+DI,j+DJ,k+DK))>0 && \
                            (Edec-Ed0(i+DI,j,k))*(Edec-Ed(i,j+DJ,k+DK))>0 && \
                            (Edec-Ed0(i,j+DJ,k))*(Edec-Ed(i+DI,j,k+DK))>0 && \
                            (Edec-Ed0(i,j,k+DK))*(Edec-Ed(i+DI,j+DJ,k))>0 && \
                            (Edec-Ed0(i+DI,j+DJ,k+DK))*(Edec-Ed(i,j,k))>0 && \
                            (Edec-Ed0(i,j+DJ,k+DK))*(Edec-Ed(i+DI,j,k))>0 && \
                            (Edec-Ed0(i+DI,j,k+DK))*(Edec-Ed(i,j+DJ,k))>0 && \
                            (Edec-Ed0(i+DI,j+DJ,k))*(Edec-Ed(i,j,k+DK))>0   )
                        intersection=false;
                    else intersection=true;

                    if(intersection==true){
                        //	cout<<"corner value"<<"\n";
                        for(int i1=0; i1<2; i1++)
                            for(int i2=0; i2<2; i2++)
                                for(int i3=0; i3<2; i3++){
                                    int id=i1*4+i2*2+i3;
                                    int I = i+i1*DI;
                                    int J = j+i2*DJ;
                                    int K = k+i3*DK;

                                    EDCorners[id]=Ed0(I,J,K);
                                    VLCorners[id][0]=Vx0(I,J,K);
                                    VLCorners[id][1]=Vy0(I,J,K);
                                    //VLCorners[id][2]=Vz0(I,J,K);  //Changed on 03/03/2011
                                    VLCorners[id][2]=Vz0(I,J,K)*tau0;

                                    EDCorners[8+id]=Ed(I,J,K);
                                    VLCorners[8+id][0]=Vx(I,J,K);
                                    VLCorners[8+id][1]=Vy(I,J,K);
                                    //VLCorners[8+id][2]=Vz(I,J,K);
                                    VLCorners[8+id][2]=Vz(I,J,K)*tau;
                                }
                        cube1.SetCornersValue(EDCorners);           //set 16 corners' energy density
                        cube1.SetOrigin(tau0, i*dx, j*dy, k*dz);    //set pos of the cube 
                        cube1.CalcInts();                           //Calc intersection points
                        cube1.SetCornersVelocity(VLCorners);
                        cube1.CalcVelocity();	 //calc the velocity through intersect
                        cube1.CalcSF2(DA0, DA1, DA2, DA3);   //calc hyper surface 

                        v1=cube1.V.x;
                        v2=cube1.V.y;
                        v3=cube1.V.z;
                        double tau_int = tau0 + DTD*cube1.MassCenter.t;
                        x = i*dx + DXD*cube1.MassCenter.x;   //Cube is 1*1*1*1
                        y = j*dy + DYD*cube1.MassCenter.y;
                        z = k*dz + DZD*cube1.MassCenter.z;
                        f_frz<<DA0 <<' '<<-DA1 <<' '<<-DA2 <<' '<<-DA3;
                        //f_frz<<' '<<v1 <<' '<<v2 <<' '<<v3*tau_int;
                        f_frz<<' '<<v1 <<' '<<v2 <<' '<<v3;
                        f_frz<<' '<<z <<"\n";

                        SF A;
                        A.DA0 = DA0;
                        A.DA1 = -DA1;
                        A.DA2 = -DA2;
                        A.DA3 = -DA3;
                        A.vx  = v1;
                        A.vy  = v2;
                        //A.vz  = v3*tau_int;
                        A.vz  = v3;
                        A.etas= z;
                        hypsf.push_back(A);


                        if(abs(z)<2*dz && abs(y)<2*dy ){
                            f_frz_taux<<x<<' '<<tau_int<<'\n';
                        }
                        if(abs(z)<2*dz && abs(x)<2*dx ){
                            f_frz_tauy<<y<<' '<<tau_int<<'\n';
                        }
                        if(abs(x)<2*dx && abs(y)<2*dy ){
                            f_frz_tauh<<z<<' '<<tau_int<<'\n';
                        }
                    }
                }
    }
}
//end CalcHyperSf2 /////////////////////////////////////////////


struct MS{ double E, Px, Py, Pz, t, x, y, z, etas; };

void CHydro3D::PPInitial()
{   char file_p4xy[100];
    //sprintf(file_p4xy,"HijingIni/moment%d.dat",IEVENT);
    sprintf(file_p4xy,"ampt_Ini/P%d.dat",IEVENT);
    double Pmu[4];
    double g[4]={1.0, -1.0, -1.0, -tau*tau};
    double sigr;
    double sigz;
    int NTOTAL;
    //double sigz=sigr;
    ED_MAX = 0.0;

    ifstream fin(file_p4xy);
    deque<MS> JETS;
    MS A;

    fin>>NTOTAL;
    while(!fin.eof()){
        //    double x,y;
        fin>>A.E>>A.Px>>A.Py>>A.Pz;
        fin>>A.t>>A.x>>A.y>>A.z;
        //    generate_xy(x,y,A.E,A.Px,A.Py,A.Pz);
        //    A.x=x;
        //    A.y=y;
        JETS.push_back(A);
    }
    JETS.pop_back();
    fin.close();

    char file_EDXY_ETA0[50];
    sprintf(file_EDXY_ETA0,"results/ED_XY%d.dat",IEVENT);
    ofstream fout(file_EDXY_ETA0);
    fout.precision(3);

    for(int i=NX0; i<=NX; i++)
        for(int j=NY0; j<=NY; j++){
            double ed_Ave=0.0;
            for(int k=NZ0; k<=NZ; k++){

                double E=0.0; 
                double Px=0.0;
                double Py=0.0;
                double Pz=0.0;
                double x= i*dx;
                double y= j*dy;
                double etas= k*dz;

                for(deque<MS>::iterator it=JETS.begin(); it!=JETS.end(); ++it){
                    double X0  =(*it).x;
                    double Y0  =(*it).y;
                    //double ETANN=max((*it).E+(*it).Pz,1.0e-15);
                    //double ETADD=max((*it).E-(*it).Pz,1.0e-15);
                    //double ETA0= 0.5*log(ETANN/ETADD);

                    double ETANN=max((*it).t+(*it).z,1.0e-15);
                    double ETADD=max((*it).t-(*it).z,1.0e-15);
                    double ETA0= 0.5*log(ETANN/ETADD);

                    double x1=(x-X0);
                    double y1=(y-Y0);
                    double etas1=etas-ETA0;
                    double z1=tau*(sinh(etas)-sinh(ETA0));
                    sigr = dx;
                    sigz = 1.0;
                    double w1 = 1.0/(2.0*PhyConst::pi*sigr*sigr);
                    double w2 = 1.0/(2.0*PhyConst::pi*sigz*sigz);  
                    //double weight=w1*sqrt(w2)*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr)-z1*z1/(2.0*sigz*sigz));
                    //double weight=w1*sqrt(w2)/tau*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr)-etas1*etas1/(2.0*sigz*sigz));
                    double acu=PhyConst::acu;
                    double weight=(x1>=0 && y1>=0 && etas1>=0 && abs(x1)<dx && abs(y1)<dy && abs(etas1)<dz)?1.0/(tau*dx*dy*dz):0.0;
                    //double weight=(abs(etas1)<=sigz) ? w1*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr))/(2.0*sigz*tau) : 0.0;
                    double mt=sqrt(max( (*it).E*(*it).E - (*it).Pz*(*it).Pz, 0.0));
                    //E += mt*weight;

                    E += (*it).E*weight;
                    Px+= (*it).Px*weight;
                    Py+= (*it).Py*weight;
                    Pz+= (*it).Pz*weight;
                }//End Summation of all T^{\mu\nu}_i
                double MT, ETA0;
                double T[4][4];
                double E2PZ = max(E*E-Pz*Pz,Px*Px+Py*Py);
                double ETANN=max(E+Pz,1.0e-15);
                double ETADD=max(E-Pz,1.0e-15);
                MT = sqrt(E2PZ);
                ETA0= 0.5*log(ETANN/ETADD);
                Pmu[0] = MT;
                Pmu[1] = Px;
                Pmu[2] = Py;
                Pmu[3] = 0.0;
                for(int m=0; m<4; m++)
                    for(int n=0; n<4; n++){
                        if(Pmu[0]<PhyConst::acu) T[m][n]=0.0;
                        else T[m][n] = Pmu[m]*Pmu[n]/Pmu[0];
                    }

                double ed;
                double Umu[4];

                SolveUmu2(T, g, Umu, ed);
                //cout<<ed<<endl;
                /*Test for Umu={1,0,0,0} case*/
                Umu[0]=1.0;
                for(int m=1; m<=3; m++)Umu[m]=0.0;
                /*0 Initial velocity test*/

                if(ED_MAX<ed)ED_MAX=ed;
                ed_Ave += ed*dz;
                Ed(i,j,k) = max(0.0, ed);

                Bd(i,j,k) = 0.0;
                Pr(i,j,k) = Eos.P(Ed(i,j,k),Bd(i,j,k));
                Umu[0] = sqrt(1.0+Umu[1]*Umu[1]+Umu[2]*Umu[2]);
                Vx(i,j,k)=Umu[1]/Umu[0];
                Vy(i,j,k)=Umu[2]/Umu[0];
                Vz(i,j,k)=Umu[3]/Umu[0];

                //fout<<x<<' '<<y<<' '<<etas<<' '<<Ed(i,j,k)<<' '<<Vx(i,j,k)<<' '<<Vy(i,j,k)<<endl;

                double epp=Ed(i,j,k)+Pr(i,j,k);
                TT00(i,j,k)=t0*(epp*Umu[0]*Umu[0]-Pr(i,j,k));
                TT01(i,j,k)=t0*(epp*Umu[0]*Umu[1]);
                TT02(i,j,k)=t0*(epp*Umu[0]*Umu[2]);
                TT03(i,j,k)=t0*(epp*Umu[0]*Umu[3]);
                TJT(i,j,k)=t0*Bd(i,j,k)*Umu[0];
            }
            //fout<<ed_Ave<<" ";
        }
    fout.close();
}



void CHydro3D::PPInitial2()//For ampt initial condition with t,x,y,z
{   char file_p4xy[100];
    //sprintf(file_p4xy,"../ampt-v1.21-v2.21/ana/ampt_b0-5_Ini/P%d.dat",IEVENT);
    sprintf(file_p4xy,"../ampt-v1.21-v2.21/ana/ampt_b20-30_Ini/P%d.dat",IEVENT);

    double Pmu[4];
    double g[4]={1.0, -1.0, -1.0, -tau*tau};
    /*sfactor=1.5, sigr=1.8dx, sigz=sigr is a group of good parameters*/
    double sfactor, sigr, sigz;

    /* For test */
    //int IE=0;
    //sprintf(file_p4xy,"../ampt-v1.21-v2.21/ana/ampt_b20-30_Ini/P%d.dat",IE);

    //if(IEVENT==0){
    // sfactor = 1.5;
    // sigr = 2.0*dx;
    // sigz = sigr;
    //}
    //else if(IEVENT==1){
    // sfactor = 1.5;
    // sigr = 2.0*dx;
    // sigz = 0.4;
    //}
    //else if(IEVENT==2){
    // sfactor = 1.5;
    // sigr = 2.0*dx;
    // sigz = 0.2;
    //}
    //else if(IEVENT==3){
    // sfactor = 1.5;
    // sigr = 2.0*dx;
    // sigz = 0.1;
    //}
    //else if(IEVENT==4){
    // sfactor = 1.5;
    // sigr = 2.0*dx;
    // sigz = 0.05;
    //}
    //else if(IEVENT==5){
    // sfactor = 2.5;
    // sigr = 2.0*dx;
    // sigz = 0.1;
    //}
    //else if(IEVENT==6){
    // sfactor = 3.0;
    // sigr = 2.0*dx;
    // sigz = 0.1;
    //}
    //else if(IEVENT==7){
    // sfactor = 3.5;
    // sigr = 2.0*dx;
    // sigz = 0.1;
    //}
    //else if(IEVENT==8){
    // sfactor = 2.0;
    // sigr = 3.0*dx;
    // sigz = 0.2;
    //}
    //else if(IEVENT==9)
    {
     sfactor = 1.5;
     sigr = 1.8*dx;
     sigz = sigr;
    }
    ////////////////////////////////////////////////




    double w1 = 1.0/(2.0*PhyConst::pi*sigr*sigr);
    double w2 = 1.0/(tau*sqrt(2.0*PhyConst::pi*sigz*sigz));  

    int NTOTAL;
    ED_MAX = 0.0;
    double acu=PhyConst::acu;

    ifstream fin(file_p4xy);
    deque<MS> JETS(20000);
    MS A;

    fin>>NTOTAL;
    while(!fin.eof()){
        fin>>A.E>>A.Px>>A.Py>>A.Pz;
        fin>>A.t>>A.x>>A.y>>A.z;
        double ETANN=max(A.t+A.z,1.0e-15);
        double ETADD=max(A.t-A.z,1.0e-15);
        A.etas = 0.5*log(ETANN/ETADD);//spatial rapidity from ampt


        //double pt=sqrt(A.Px*A.Px+A.Py*A.Py);
        //double phi=2.0*PhyConst::pi*ran2(&NSEED);
        //A.Px = pt*cos(phi);
        //A.Py = pt*sin(phi);
        //double r=sqrt(A.x*A.x+A.y*A.y);
        //A.x = r*cos(phi);
        //A.y = r*sin(phi);
        JETS.push_back(A);
    }
    JETS.pop_back();
    fin.close();


    char file_EDXY_ETA0[50];
    sprintf(file_EDXY_ETA0,"results/ED_XY%d.dat",IEVENT);
    ofstream fout(file_EDXY_ETA0);
    fout.precision(3);

    for(int i=NX0; i<=NX; i++)
        for(int j=NY0; j<=NY; j++){
            double ed_Ave=0.0;
            for(int k=NZ0; k<=NZ; k++){

                double E=0.0; 
                double Px=0.0;
                double Py=0.0;
                double Pz=0.0;
                double x= i*dx;
                double y= j*dy;
                double etas= k*dz;

                for(deque<MS>::iterator it=JETS.begin(); it!=JETS.end(); ++it){

                    double x1=(x-(*it).x);
                    double y1=(y-(*it).y);
                    double etas1=etas-(*it).etas;
                    //double z1=tau*(sinh(etas)-sinh(ETA0));
                    //double pt=sqrt((*it).Px*(*it).Px + (*it).Py*(*it).Py);
                    //sigr = PhyConst::hbarc/max(0.4,pt);
                    //sigz = dz;



                    double weight = ((x1*x1 + y1*y1 +tau*tau*etas1*etas1)<4.0) ? 1.0 : 0.0;
                    weight *= w1*w2*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr)-etas1*etas1/(2.0*sigz*sigz));

                    //double weight=w1*sqrt(w2)*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr)-z1*z1/(2.0*sigz*sigz));
                    //double weight=w1*sqrt(w2)/tau*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr)-etas1*etas1/(2.0*sigz*sigz));
                    //double weight=(abs(x1)<=0.5*dx && abs(y1)<=0.5*dy && abs(etas1)<=0.5*dz)?1.0/(tau*dx*dy*dz):0.0;
                    //double weight=(abs(x1)<=dx && abs(y1)<=dy && abs(etas1)<=dz)?1.0/(8.0*tau*dx*dy*dz):0.0;
                    //double weight=(abs(etas1)<=sigz) ? w1*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr))/(2.0*sigz*tau) : 0.0;
                    E += sfactor*(*it).E*weight;
                    Px+= sfactor*(*it).Px*weight;
                    Py+= sfactor*(*it).Py*weight;
                    Pz+= sfactor*(*it).Pz*weight;
                }//End Summation of all T^{\mu\nu}_i
                double MT, RY0;
                double T[4][4];
                double E2PZ = max(E*E-Pz*Pz,Px*Px+Py*Py);
                double ETANN=max(E+Pz,1.0e-15);
                double ETADD=max(E-Pz,1.0e-15);
                MT = sqrt(E2PZ);
                RY0= 0.5*log(ETANN/ETADD);


                Pmu[0] = MT*cosh(RY0-etas);
                Pmu[1] = Px;
                Pmu[2] = Py;
                Pmu[3] = MT*sinh(RY0-etas)/tau;
                //if(Pmu[0]>PhyConst::acu)
                //fout<<Pmu[0]<<' '<<Pmu[1]<<' '<<Pmu[2]<<' '<<Pmu[3]<<endl;
                for(int m=0; m<4; m++)
                    for(int n=0; n<4; n++){
                        if(Pmu[0]<PhyConst::acu) T[m][n]=0.0;
                        else T[m][n] = Pmu[m]*Pmu[n]/Pmu[0];
                    }

                double ed;
                double Umu[4];

                SolveUmu3(T, g, Umu, ed);

                //Umu[1]=0.0;     /*No initial velocity to test the slope of the pt spectra*/
                //Umu[2]=0.0;
                //Umu[3]=0.0;


                //for(int m=1; m<=2; m++)Umu[m]=0;
                //cout<<ed<<endl;
                if(ED_MAX<ed)ED_MAX=ed;
                ed_Ave += ed*dz;
                Ed(i,j,k) = max(0.0, ed);

                Bd(i,j,k) = 0.0;
                Pr(i,j,k) = Eos.P(Ed(i,j,k),Bd(i,j,k));
                Umu[0] = sqrt(1.0-g[1]*Umu[1]*Umu[1]-g[2]*Umu[2]*Umu[2]-g[3]*Umu[3]*Umu[3]);
                Vx(i,j,k)=Umu[1]/Umu[0];
                Vy(i,j,k)=Umu[2]/Umu[0];
                Vz(i,j,k)=Umu[3]/Umu[0];  //here Vz=Vetas/tau

                //if(Ed(i,j,k)>PhyConst::acu)
                //	fout<<i<<' '<<j<<' '<<k<<' '<<Ed(i,j,k)<<endl;

                double epp=Ed(i,j,k)+Pr(i,j,k);
                TT00(i,j,k)=t0*(epp*Umu[0]*Umu[0]-Pr(i,j,k));
                TT01(i,j,k)=t0*(epp*Umu[0]*Umu[1]);
                TT02(i,j,k)=t0*(epp*Umu[0]*Umu[2]);
                TT03(i,j,k)=t0*(epp*Umu[0]*Umu[3]);
                TJT(i,j,k)=t0*Bd(i,j,k)*Umu[0];
            }
        }
    fout.close();
}


struct NWN{ double x,y,z,ncol; };
//For MCG initial condition with x,y,z and ncol
void CHydro3D::PPInitial4()
{
    char file_p4xy[100];
    sprintf(file_p4xy,"../MCG/MCG_INI_b0-10/nwn_%d.dat",IEVENT);
    ifstream fin(file_p4xy);
    deque<NWN> NWN0(1000);
     double sfactor, sigr, sigz, gw;
    {
     sfactor = 120;
     sigr = 0.4;
     sigz = 0.4;
     gw = 2.95;
    }
 
    double w1 = 1.0/(2.0*PhyConst::pi*sigr*sigr);
    double w2 = 1.0/(2.0*PhyConst::pi*sigz*sigz);
    ED_MAX = 0.0;
    double acu=PhyConst::acu;

    NWN A;
    while(!fin.eof()){
       char nucleus;
       fin>>nucleus>>A.x>>A.y>>A.z>>A.ncol;
       NWN0.push_back(A);
    }
    NWN0.pop_back();

     for(int i=NX0; i<=NX; i++)
        for(int j=NY0; j<=NY; j++)
            for(int k=NZ0; k<=NZ; k++){

                double x= i*dx;
                double y= j*dy;
                double etas= k*dz;

                for(deque<NWN>::iterator it=NWN0.begin(); it!=NWN0.end(); ++it){
                    double x1=(x-(*it).x);
                    double y1=(y-(*it).y);
                    double weight = (abs(etas)<gw) ? 1.0 : exp(-(abs(etas)-gw)*(abs(etas)-gw)/(2.0*sigz*sigz));
                    weight *= w1*w2*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr));
                    Ed(i,j,k) += sfactor*weight;
                }
                if(ED_MAX<Ed(i,j,k))ED_MAX=Ed(i,j,k);

                Pr(i,j,k)=Eos.P( Ed(i,j,k),Bd(i,j,k) );

                Vx(i,j,k)=0.0;    //Test if a smaller v2 is caused by initial velocity
                Vy(i,j,k)=0.0;
                Vz(i,j,k)=0.0;
                double Umu[4];
                Umu[0]=Gamma(i,j,k);
                Umu[1]=Umu[0]*Vx(i,j,k);
                Umu[2]=Umu[0]*Vy(i,j,k);
                Umu[3]=Umu[0]*Vz(i,j,k);

                double epp=Ed(i,j,k)+Pr(i,j,k);
                TT00(i,j,k)=t0*(epp*Umu[0]*Umu[0]-Pr(i,j,k));
                TT01(i,j,k)=t0*(epp*Umu[0]*Umu[1]);
                TT02(i,j,k)=t0*(epp*Umu[0]*Umu[2]);
                TT03(i,j,k)=t0*(epp*Umu[0]*Umu[3]);
                TJT(i,j,k)=t0*Bd(i,j,k)*Umu[0];

                }
 
}
////////////////////////////////////////////////////////////////////////////
void CHydro3D::PPInitial3()//For ampt initial condition with t,x,y,z
{   char file_p4xy[100];
    //sprintf(file_p4xy,"../ampt-v1.21-v2.21/ana/ampt_b0-5_Ini/P%d.dat",IEVENT);
    sprintf(file_p4xy,"../ampt-v1.21-v2.21/ana/ampt_b0-80_Ini/P%d.dat",IEVENT);

    double Pmu[4];
    double g[4]={1.0, -1.0, -1.0, -tau*tau};
    /*sfactor=1.5, sigr=1.8dx, sigz=sigr is a group of good parameters*/
    double sfactor, sigr, sigz;
    {
     sfactor = 1.6;
     sigr = 2.0*dx;
     sigz = sigr;
    }
    ////////////////////////////////////////////////

    double w1 = 1.0/(2.0*PhyConst::pi*sigr*sigr);
    double w2 = 1.0/(tau*sqrt(2.0*PhyConst::pi*sigz*sigz));  

    int NTOTAL;
    ED_MAX = 0.0;
    double acu=PhyConst::acu;

    ifstream fin(file_p4xy);
    deque<MS> JETS(20000);
    MS A;

    fin>>NTOTAL;
    while(!fin.eof()){
        fin>>A.E>>A.Px>>A.Py>>A.Pz;
        fin>>A.t>>A.x>>A.y>>A.z;
        double ETANN=max(A.t+A.z,1.0e-15);
        double ETADD=max(A.t-A.z,1.0e-15);
        A.etas = 0.5*log(ETANN/ETADD);//spatial rapidity from ampt
        double pt=sqrt(A.Px*A.Px+A.Py*A.Py);
        A.x = A.x + 1.0*A.Px/pt;
        A.y = A.y + 1.0*A.Py/pt;
        //if(abs(A.etas)<3.0)A.etas *= 1.3;
        double sign= (A.etas>0.0) ? 1.0 : -1.0;
        //if(abs(A.etas)<3.0)A.etas *= 1.0+0.5*cos(A.etas*PhyConst::pi/6.0);
        if(abs(A.etas)<4.0)A.etas *= 1.0+0.75*(1.0-abs(A.etas)/4.0);
        //else if(abs(A.etas)>3.0)A.etas = sign*(3.0 + (abs(A.etas)-3.0)*0.2);
        //else if(abs(A.etas)>5.0)A.etas *= 1.0+0.3*cos((4.5-abs(A.etas))/3.0*PhyConst::pi);


        //double pt=sqrt(A.Px*A.Px+A.Py*A.Py);
        //double phi=2.0*PhyConst::pi*ran2(&NSEED);
        //A.Px = pt*cos(phi);
        //A.Py = pt*sin(phi);
        //double r=sqrt(A.x*A.x+A.y*A.y);
        //A.x = r*cos(phi);
        //A.y = r*sin(phi);
        JETS.push_back(A);
    }
    JETS.pop_back();
    fin.close();

    double Psi[7]={0};
    double Ecc[7]={0};
    double Ar2cos[7]={0.0};   //<r^2 * cos(n\phi)>
    double Ar2sin[7]={0.0};   //<r^2 * sin(n\phi)>

    for(int i=NX0; i<=NX; i++)
        for(int j=NY0; j<=NY; j++){
            double ed_Ave=0.0;
            for(int k=NZ0; k<=NZ; k++){

                double E=0.0; 
                double Px=0.0;
                double Py=0.0;
                double Pz=0.0;
                double x= i*dx;
                double y= j*dy;
                double etas= k*dz;

                for(deque<MS>::iterator it=JETS.begin(); it!=JETS.end(); ++it){

                    double x1=(x-(*it).x);
                    double y1=(y-(*it).y);
                    double etas1=etas-(*it).etas;
                    double z1=tau*sinh(etas)-(*it).z;
                    //double z1=tau*(sinh(etas)-sinh(ETA0));
                    //double pt=sqrt((*it).Px*(*it).Px + (*it).Py*(*it).Py);
                    //sigr = PhyConst::hbarc/max(0.4,pt);
                    //sigz = dz;



                    double weight = ((x1*x1 + y1*y1 +tau*tau*etas1*etas1)<4.0) ? 1.0 : 0.0;
                    //double weight=1.0;
                    weight *= w1*w2*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr)-etas1*etas1/(2.0*sigz*sigz));

                    //double weight=w1*sqrt(w2)*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr)-z1*z1/(2.0*sigz*sigz));
                    //double weight=w1*sqrt(w2)/tau*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr)-etas1*etas1/(2.0*sigz*sigz));
                    //double weight=(abs(x1)<=0.5*dx && abs(y1)<=0.5*dy && abs(etas1)<=0.5*dz)?1.0/(tau*dx*dy*dz):0.0;
                    //double weight=(abs(x1)<=dx && abs(y1)<=dy && abs(etas1)<=dz)?1.0/(8.0*tau*dx*dy*dz):0.0;
                    //double weight=(abs(etas1)<=sigz) ? w1*exp(-(x1*x1+y1*y1)/(2.0*sigr*sigr))/(2.0*sigz*tau) : 0.0;
                    E += sfactor*(*it).E*weight;
                    Px+= sfactor*(*it).Px*weight;
                    Py+= sfactor*(*it).Py*weight;
                    Pz+= sfactor*(*it).Pz*weight;
                }//End Summation of all T^{\mu\nu}_i
                double MT, RY0;
                double T[4][4];
                double E2PZ = max(E*E-Pz*Pz,Px*Px+Py*Py);
                double ETANN=max(E+Pz,1.0e-15);
                double ETADD=max(E-Pz,1.0e-15);
                MT = sqrt(E2PZ);
                RY0= 0.5*log(ETANN/ETADD);


                Pmu[0] = MT*cosh(RY0-etas);
                Pmu[1] = Px;
                Pmu[2] = Py;
                Pmu[3] = MT*sinh(RY0-etas)/tau;
                //if(Pmu[0]>PhyConst::acu)
                //fout<<Pmu[0]<<' '<<Pmu[1]<<' '<<Pmu[2]<<' '<<Pmu[3]<<endl;
                for(int m=0; m<4; m++)
                    for(int n=0; n<4; n++){
                        if(Pmu[0]<PhyConst::acu) T[m][n]=0.0;
                        else T[m][n] = Pmu[m]*Pmu[n]/Pmu[0];
                    }

                double E0=0.0; 
                double B0=0.0;
                double T00=T[0][0];
                double T01=T[0][1];
                double T02=T[0][2];
                double T03=T[0][3]*tau*tau;
                double JT0=0.0;
                Rootfinding(E0, B0, T00, T01, T02, T03, JT0);
                Ed(i,j,k)=max(0.0, E0);//here acu changed to 0.0 4/27
                if(ED_MAX<Ed(i,j,k))ED_MAX=Ed(i,j,k);

                Bd(i,j,k)=max(0.0, B0);
                Pr(i,j,k)=Eos.P( Ed(i,j,k),Bd(i,j,k) );
                double EPV=max( acu, T00+Pr(i,j,k) );
                Vx(i,j,k)=T01/EPV;
                Vy(i,j,k)=T02/EPV;
                Vz(i,j,k)=T03/EPV/tau;

                //Vx(i,j,k)=0.0;    //Test if a smaller v2 is caused by initial velocity
                //Vy(i,j,k)=0.0;
                //Vz(i,j,k)=0.0;
                double Umu[4];
                Umu[0]=Gamma(i,j,k);
                Umu[1]=Umu[0]*Vx(i,j,k);
                Umu[2]=Umu[0]*Vy(i,j,k);
                Umu[3]=Umu[0]*Vz(i,j,k);

/*
                double ed;
                double Umu[4];
                SolveUmu3(T, g, Umu, ed);
                //for(int m=1; m<=2; m++)Umu[m]=0;
                //cout<<ed<<endl;
                if(ED_MAX<ed)ED_MAX=ed;
                ed_Ave += ed*dz;
                Ed(i,j,k) = max(0.0, ed);

                Bd(i,j,k) = 0.0;
                Pr(i,j,k) = Eos.P(Ed(i,j,k),Bd(i,j,k));
                Umu[0] = sqrt(1.0-g[1]*Umu[1]*Umu[1]-g[2]*Umu[2]*Umu[2]-g[3]*Umu[3]*Umu[3]);
                Vx(i,j,k)=Umu[1]/Umu[0];
                Vy(i,j,k)=Umu[2]/Umu[0];
                Vz(i,j,k)=Umu[3]/Umu[0];  //here Vz=Vetas/tau

                //if(Ed(i,j,k)>PhyConst::acu)
                //	fout<<i<<' '<<j<<' '<<k<<' '<<Ed(i,j,k)<<endl;
*/
                //if(Ed(i,j,k)>PhyConst::acu)
                double epp=Ed(i,j,k)+Pr(i,j,k);
                TT00(i,j,k)=t0*(epp*Umu[0]*Umu[0]-Pr(i,j,k));
                TT01(i,j,k)=t0*(epp*Umu[0]*Umu[1]);
                TT02(i,j,k)=t0*(epp*Umu[0]*Umu[2]);
                TT03(i,j,k)=t0*(epp*Umu[0]*Umu[3]);
                TJT(i,j,k)=t0*Bd(i,j,k)*Umu[0];

                double r2 = x*x + y*y;
                double phiS=acos(x/max(1.0E-15,sqrt(r2)));
                if(y<0) phiS = 2.0*pi-phiS;
                for(int n=1; n<=6; n++){
                    Ar2sin[n] += r2*sin(n*phiS)*Ed(i,j,k);
                    Ar2cos[n] += r2*cos(n*phiS)*Ed(i,j,k);
                }
            }
        }

    for(int n=1; n<=6; n++){
        Psi[n] = (atan(2.0*Ar2sin[n]/Ar2cos[n]))/double(n);  
    }

    double Ar2sin_ecc[7]={0.0};
    double Ar2cos_ecc[7]={0.0};
    double Ar2_ecc=0.0;

    for(int i=NX0; i<=NX; i++){
        double x = i*dx;
        for(int j=NY0; j<=NY; j++){
            double y = j*dy;
            for(int k=NZ0; k<=NZ; k++){
                double r2 = x*x + y*y;
                double phiS=acos(x/max(1.0E-15,sqrt(r2)));
                if(y<0) phiS = 2.0*pi-phiS;
                for(int n=1; n<=6; n++){
                    Ar2sin_ecc[n] += r2*sin(n*(phiS-Psi[n]))*Ed(i,j,k);
                    Ar2cos_ecc[n] += r2*cos(n*(phiS-Psi[n]))*Ed(i,j,k);
                }
                Ar2_ecc += r2*Ed(i,j,k);
    }
    }
    }//End for  i,j,k


    char file_Psi[50];
    sprintf(file_Psi,"results/event%d/Psi.dat",IEVENT);
    ofstream fPsi(file_Psi);
    fPsi.precision(6);

    for(int n=1; n<=6; n++){
       Ecc[n] = sqrt( Ar2sin_ecc[n]*Ar2sin_ecc[n] + Ar2cos_ecc[n]*Ar2cos_ecc[n] )/Ar2_ecc;
       fPsi<<Psi[n]<<' '<<Ecc[n]<<endl;
    }
    fPsi.close();


    char file_EDXY_ETA0[50];
    sprintf(file_EDXY_ETA0,"results/event%d/ED_XY.dat",IEVENT);
    ofstream fout(file_EDXY_ETA0);
    fout.precision(3);

    for(int i=NX0; i<=NX; i++)
        for(int j=NY0; j<=NY; j++){
                	fout<<i<<' '<<j<<' '<<0<<' '<<Ed(i,j,0)<<endl;
    }
    fout.close();


    char file_EDXH_ETA0[50];
    sprintf(file_EDXH_ETA0,"results/event%d/ED_XH.dat",IEVENT);
    ofstream fout0(file_EDXH_ETA0);
    fout0.precision(3);

    for(int i=NX0; i<=NX; i++)
        for(int j=NZ0; j<=NZ; j++){
                	fout0<<i<<' '<<j<<' '<<0<<' '<<Ed(i,0,j)<<endl;
    }
    fout0.close();

    char file_VT_ETA0[50];
    sprintf(file_VT_ETA0,"results/event%d/VT_XY.dat",IEVENT);
    ofstream fout1(file_VT_ETA0);
    fout1.precision(3);

    for(int i=NX0; i<=NX; i++)
        for(int j=NY0; j<=NY; j++){
                	fout1<<i<<' '<<j<<' '<<Vx(i,j,0)<<' '<<Vy(i,j,0)<<endl;
    }
    fout1.close();

    char file_VXH_ETA0[50];
    sprintf(file_VXH_ETA0,"results/event%d/VXH_XH.dat",IEVENT);
    ofstream fout2(file_VXH_ETA0);
    fout2.precision(3);

    for(int i=NX0; i<=NX; i++)
        for(int j=NZ0; j<=NZ; j++){
                	fout2<<i<<' '<<j<<' '<<Vx(i,j,0)<<' '<<Vz(i,0,j)<<endl;
    }
    fout2.close();



}
