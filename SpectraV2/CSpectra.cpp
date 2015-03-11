#include"CSpectra.h"
using namespace std;

CSpectra::CSpectra()
{
    NEVENT = 1;
    frz_temp = 0.137; //Default frz out temp value
}

CSpectra::CSpectra(CI & NEVENT_IN)
{
    NEVENT = NEVENT_IN; 
    frz_temp = 0.137; //Default frz out temp value
}

void CSpectra::Set_FrzT(CD &TFRZ)
{
    frz_temp = TFRZ;
}/* Set Frz out temp */


void CSpectra::ReadHyperSF(CI & IEVENT)
{
    char fhypsf[256];
    sprintf(fhypsf,"../results/HyperSf%d.dat",IEVENT);
    ifstream fin(fhypsf);
    SF s;
    while(!fin.eof()){
        fin>>s.DA0>>s.DA1>>s.DA2>>s.DA3>>s.vx>>s.vy>>s.vz>>s.etas;
        hypsf.push_back(s);
    }
    hypsf.pop_back();
    fin.close();
}/*Read Frz out Hyper surface data */



void CSpectra::Calc_dNdYPtdPtdPhi(CI &IEVENT, CParticle &par)
{
    /*Calc dNdYptdptdphi for particle par Before Reso Decay*/
    if(par.monval > 0){
        //Only cal positive particles, the anti-particles has the same spectra

        CD mass=par.mass;
        CD fermi_boson = (abs(par.baryon)) ? +1.0 : -1.0; // Fermi or boson
        //CD mu_b=(abs(par.baryon)) ? 0.130 : 0.0;              // To fix the pi/p ratio
        CD mu_b=0.0;

        int MC=abs(par.monval);
        CD t_dec=frz_temp;//freeze out temperature
        //CD t_dec=Eos.T(setting.edec,0.0);
        cout<<"t_dec="<<t_dec<<endl;

        CD g = par.gspin;
        CD dof = 1/pow(2.0*PhyConst::pi,3.0);
        char filename[128];
        if(par.monval>=0)sprintf(filename,"../results/event%d/dNdYPtdPtdPhi_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename,"../results/event%d/dNdYPtdPtdPhi_A%d.dat",IEVENT,MC);

        ofstream f_dNdYdPtdPhi(filename);

        CD rylo=-8.0;
        CD ryhi= 8.0;
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
        if(par.monval>=0)sprintf(filename_dNdY,"../results/event%d/dNdY_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdY,"../results/event%d/dNdY_A%d.dat",IEVENT,MC);
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
        //        sprintf(fantib,"../results/dNdYPtdPtdPhi_A%d.dat",abs(par.monval));
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


enum DIR_OR_RESO{DIRECT, RESO}; //{0, 1}
void CSpectra::dNdEta(CI &IEVENT, CI &DIR_OR_RESO, CParticle &par)
{
    /*Calc dNdEta, dNdPt, v2 for particle par 
      Before Reso Decay
      or After Reso Decay*/
    //DIR_OR_RESO = 0 for direct produced particles
    //DIR_OR_RESO = 1 for after resonance decay
    CD mass=par.mass;
    int MC=abs(par.monval);
    CD dof=1.0/pow(2.0*PhyConst::pi,3.0);
    char filename[128];

    double rylo, ryhi, etalo, etahi;

    if(DIR_OR_RESO == 1){
        rylo=-8.0;   
        ryhi= 8.0;
        if(par.monval>=0)sprintf(filename,"../results/event%d/spec_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename,"../results/event%d/spec_A%d.dat",IEVENT,MC);
    }
    else{
        rylo=-8.0;   //Before resance decay y from -8 to 8
        ryhi= 8.0;
        if(par.monval>=0)sprintf(filename,"../results/event%d/dNdYPtdPtdPhi_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename,"../results/event%d/dNdYPtdPtdPhi_A%d.dat",IEVENT,MC);
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
        if(par.monval>=0)sprintf(filename_dNdY,"../results/event%d/dNdY_Reso_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdY,"../results/event%d/dNdY_Reso_A%d.dat",IEVENT,MC);
    }
    else{
        if(par.monval>=0)sprintf(filename_dNdY,"../results/event%d/dNdY_Dirc_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdY,"../results/event%d/dNdY_Dirc_A%d.dat",IEVENT,MC);
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
        if(par.monval>=0)sprintf(filename_dNdPt,"../results/event%d/dNdPt_Reso_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdPt,"../results/event%d/dNdPt_Reso_A%d.dat",IEVENT,MC);
    }
    else{
        if(par.monval>=0)sprintf(filename_dNdPt,"../results/event%d/dNdPt_Dirc_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdPt,"../results/event%d/dNdPt_Dirc_A%d.dat",IEVENT,MC);
    }
    ofstream f_dNdPt(filename_dNdPt);

    if(DIR_OR_RESO == 1){ //After resonance decay
        if(par.monval>=0)sprintf(filename_v2,"../results/event%d/v2_Reso_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_v2,"../results/event%d/v2_Reso_A%d.dat",IEVENT,MC);
    }
    else{
        if(par.monval>=0)sprintf(filename_v2,"../results/event%d/v2_Dirc_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_v2,"../results/event%d/v2_Dirc_A%d.dat",IEVENT,MC);
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
        if(par.monval>=0)sprintf(filename_dNdPhi,"../results/event%d/dNdPhi_Reso_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdPhi,"../results/event%d/dNdPhi_Reso_A%d.dat",IEVENT,MC);
    }
    else{
        if(par.monval>=0)sprintf(filename_dNdPhi,"../results/event%d/dNdPhi_Dirc_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdPhi,"../results/event%d/dNdPhi_Dirc_A%d.dat",IEVENT,MC);
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
        if(par.monval>=0)sprintf(filename_dNdEta,"../results/event%d/dNdEta_Reso_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdEta,"../results/event%d/dNdEta_Reso_A%d.dat",IEVENT,MC);
    }
    else{
        if(par.monval>=0)sprintf(filename_dNdEta,"../results/event%d/dNdEta_Dirc_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdEta,"../results/event%d/dNdEta_Dirc_A%d.dat",IEVENT,MC);
    }
    ofstream f_dNdEta(filename_dNdEta);


    char filename_dNdEtaPtdPtdPhi[128];
    if(DIR_OR_RESO == 1){ //After resonance decay
        etalo=-8.0;
        etahi= 8.0;
        if(par.monval>=0)sprintf(filename_dNdEtaPtdPtdPhi,"../results/event%d/dNdEtaPtdPtdPhi_Reso_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdEtaPtdPtdPhi,"../results/event%d/dNdEtaPtdPtdPhi_Reso_A%d.dat",IEVENT,MC);
    }
    else{
        etalo=-8.0;
        etahi= 8.0;
        if(par.monval>=0)sprintf(filename_dNdEtaPtdPtdPhi,"../results/event%d/dNdEtaPtdPtdPhi_Dirc_%d.dat",IEVENT,MC);
        else if(par.monval<0)sprintf(filename_dNdEtaPtdPtdPhi,"../results/event%d/dNdEtaPtdPtdPhi_Dirc_A%d.dat",IEVENT,MC);
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

//////////////////////////////////////////////////////////////////////
void CSpectra::ReadParticles(char * particle_data_table)
{
    double buff;
    char oneline[256];

    ifstream fin(particle_data_table);
    CParticle p;

    while(!fin.eof()){
        p.stable = 0;
        fin>>p.monval>>p.name>>p.mass>>p.width         \
            >>p.gspin>>p.baryon>>p.strange>>p.charm     \
            >>p.bottom>>p.gisospin>>p.charge>>p.decays;

        if(particle_data_table == "../Reso3D/resoweak.dat"){
            fin>>buff>>buff>>buff; 
        }//Not used in s95p-v1/pdg05.dat
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
//////////////////////////////////////////////////////////////////////

void CSpectra::CalcSpectra_BeforeReso(CI & IEVENT)
{

    for(deque<CParticle>::size_type i=0; i!=particles.size();++i){
        Calc_dNdYPtdPtdPhi(IEVENT, particles[i]);
        dNdEta(IEVENT, DIRECT, particles[i]);
    }

}


void CSpectra::CalcSpectra_AfterReso(CI & IEVENT)
{
    for(deque<CParticle>::size_type i=0; i!=particles.size();++i){
        dNdEta(IEVENT, RESO, particles[i]);
    }
}


/**calculate dNdEta for all the charged particles*/
void CSpectra::dNdEta_Charged(CI &IEVENT)
{
    double Eta[NETA];
    double dNdEta_Charged[NETA]={0.0};
    double dNdEtaPtdPtdPhi_Charged[NETA][NPT][NPHI]={0.0};

    /*dNdEta for charged particles */
    char fname[128];
    sprintf(fname,"../results/event%d/dNdEta_Charged.dat",IEVENT);
    ofstream f_dNdEta_Charged(fname);
    /////////////////////////////////

    /*dNdEtaPtdPtdPhi for charged particles */
    sprintf(fname,"../results/event%d/dNdEtaPtdPtdPhi_Charged.dat",IEVENT);
    ofstream f_dNdEtaPtdPtdPhi_Charged(fname);
    /////////////////////////////////////////

    for(deque<CParticle>::size_type i=0; i!=particles.size();++i){

        if(particles[i].stable && abs(particles[i].charge)){
            for(int k=0; k!=NETA; k++){ 
                Eta[k] = particles[i].ETA[k];
                dNdEta_Charged[k] += particles[i].dNdEta[k];
            }

            char fcname[128];
            sprintf(fcname,"../results/event%d/spec_%d.dat",IEVENT,particles[i].monval);
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

//#define __TEST_SPECTRA__
#ifdef __TEST_SPECTRA__
int main(int argc, char** argv)
{
    CSpectra spec(1);
    spec.ReadParticles("../Reso3D/resoweak.dat");
    for(deque<CParticle>::size_type i=0; i!=spec.particles.size();++i){
        int m=spec.particles[i].monval;
        cout<<spec.particles[i].monval<<' '<<endl;
    }

    spec.Set_FrzT(0.137);
    spec.ReadHyperSF(6);
    spec.CalcSpectra_BeforeReso(6);
    //spec.CalcSpectra_AfterReso(6);
    //spec.dNdEta_Charged(6);
}
#endif


