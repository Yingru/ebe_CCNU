#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<cassert>
#include<deque>
#include"Int.h"

//#define GET_PSI_FROM_HYDRO


typedef const int CI;
typedef const double CD;
CI NEVENT=100;       //Number of events
CI NH = 6;                              //For harmonic flow vn, n from 1 through NH
//double eps_n[NH+1];                     //Initial geometric irregularities \eps_n , n=1,2,...NH
//double psi_n[NH+1];                     //Event plan angle for each harmonic flow   n=1,2,...NH

CI NPT=15;                            //From hydro output
CI NPHI=48;                           //From hydro output
CI NETA=41;                           //From hydro output
CD eta_lo=-5.0;
CD eta_hi=5.0;
CD pi=3.141592653;
CD deta = (eta_hi-eta_lo)/(NETA-1.0);
CD dphi = (2.0*pi)/40.0;

CD trig_lo=3.0;
CD trig_hi=4.0;
CD asso_lo=2.0;
CD asso_hi=3.0;

CD slope=12.0;

double ecc2=0.0;
double pt[NPT];
double phi[NPHI];
double eta[NETA];
double VN[NH+1][NPT]={0.0};
using namespace std;

/***********************************************************************/
struct CJet{ double E, Px, Py, Pz, t, x, y, z; };



void Calc_Event_Plan(const int &IEVENT, double psi_n[NH+1], double eps_n[NH+1]){
    /* Calc eps_n and psi_n in paper arXiv:1011.5249  */
    /* eps_n is the initial geometric eccentricity    */
    /* psi_n is event plan angle of the harmonic flow */
    char file_ini[64];
    sprintf(file_ini,"../ampt-v1.21-v2.21/ana/ampt_b0-10_Ini/P%d.dat",IEVENT);
    ifstream fin(file_ini);
    int NTOTAL;

    deque<CJet> jets;
    fin>>NTOTAL;
    CJet A;
    double Ar2cos[NH+1]={0.0};   //<r^2 * cos(n\phi)>
    double Ar2sin[NH+1]={0.0};   //<r^2 * sin(n\phi)>
    double Ar2   =0.0;   //<r^2 >

    double Ar2cos_ecc[NH+1]={0.0};   //<r^2 * cos(n\phi-n\Psi_n)>
    double Ar2sin_ecc[NH+1]={0.0};   //<r^2 * sin(n\phi-n\Psi_n)>

    double Apc[NH+1]={0.0};   //<pt*cos(n phi)>
    double Aps[NH+1]={0.0};   //<pt*cos(n phi)>

    for(int i=0; i!=NTOTAL; i++){
        fin>>A.E>>A.Px>>A.Py>>A.Pz>>A.t>>A.x>>A.y>>A.z;
        jets.push_back(A);
        double r2 = A.x*A.x + A.y*A.y;
        double phiS=acos(A.x/max(1.0E-15,sqrt(r2)));
        if(A.y<0) phiS = 2.0*pi-phiS;

        double pt=A.Px*A.Px + A.Py*A.Py;
        double phi=acos(A.Px/max(1.0E-15,sqrt(pt)));
        if(A.Py<0) phi=2.0*pi-phi;

        double tmz = max(1.0E-15, A.t-A.z);
        double tpz = max(1.0E-15, A.t+A.z);
        double etas=0.5*log(tpz/tmz);

        //if( abs(etas)<1.0 ){
        for(int n=1; n<=NH; n++){
            double cnp=cos(n*phi);
            double snp=sin(n*phi);
            double pta=pow(pt,0.25);
            //double pta=pow(pt,1.00);
            Ar2cos[n] += r2*cos(n*phiS);
            Ar2sin[n] += r2*sin(n*phiS);
            Apc[n] += pta*cnp;
            Aps[n] += pta*snp;
        }
        Ar2  += r2;
        //}
    }

#ifdef GET_PSI_FROM_HYDRO /*Calced from smoothed hydro initial condition*/
    char Psi_new[64];
    sprintf(Psi_new,"results/event%d/Psi_new.dat",IEVENT);
    ifstream fpsi(Psi_new);
    for(int n=1; n<=NH; n++){
        fpsi>>psi_n[n]>>eps_n[n];
    }
    fpsi.close();

#else          /*Calc from AMPT initial condition */
    for(int n=1; n<=NH; n++){
        psi_n[n] = (atan(2.0*Ar2sin[n]/Ar2cos[n]))/double(n);  
        //fpsi>>psi_n[n];
        //psi_n[n] = 0.0;
        //eps_n[n] = sqrt(Ar2cos[n]*Ar2cos[n]+Ar2sin[n]*Ar2sin[n])/Ar2;
    }

    /* Calc ecc_n **************************************************/
    for(deque<CJet>::size_type i=0; i!=NTOTAL; i++){
        A = jets[i];
        double r2 = A.x*A.x + A.y*A.y;
        double phiS=acos(A.x/max(1.0E-15,sqrt(r2)));
        if(A.y<0) phiS = 2.0*pi-phiS;
        for(int n=1; n<=NH; n++){
            Ar2cos_ecc[n] += r2*cos(n*(phiS-psi_n[n]));
            Ar2sin_ecc[n] += r2*sin(n*(phiS-psi_n[n]));
        }
    }


    for(int n=1; n<=NH; n++){
        eps_n[n] = sqrt(Ar2cos_ecc[n]*Ar2cos_ecc[n]+Ar2sin_ecc[n]*Ar2sin_ecc[n])/Ar2;
    }
#endif   /* End ifdef GET_PSI_FROM_HYDRO */
}
/***********************************************************************/

double dNdp3(CD &eta1, CD &pt1, CD &phi1, double M[NETA][NPT][NPHI]);

void ReadSpectra(CI &IEVENT, double dNdEtaPtdPtdPhi[NETA][NPT][NPHI]){
    char fname[128];
    //sprintf(fname,"results/ampt1/event%d/dNdEtaPtdPtdPhi_Charged.dat",IEVENT);
    //sprintf(fname,"results/event%d/dNdEtaPtdPtdPhi_Charged.dat",IEVENT);
    //sprintf(fname,"results/event%d/dNdYPtdPtdPhi_211.dat",IEVENT);
    sprintf(fname,"results/event%d/spec_211.dat",IEVENT);
    ifstream fin(fname);
    double dum;
    for(int i=0; i<NETA; i++)
        for(int j=0; j<NPT; j++)
            for(int k=0; k<NPHI; k++){
                fin>>dum;
                dNdEtaPtdPtdPhi[i][j][k] = dum;
            }
}// Read spectra from the hydro output
//////////////////////////////////////////////////////////////////////

/*interpolation dNdYptdptdphi in phi direction and integral from ptCut_l to ptCut_h to
  generate triger and associate dNdEtadPhi with constant deta and dphi .
  Note the input pt and phi are gaussian nodes while the output eta and phi are constant 
  separations*/
void dNdEtadPhi(double M[NETA][NPT][NPHI], double N[41][41], CD & ptCut_l, CD & ptCut_h){

    const int npt=24;
    double *p_pt = gaulep48;                           //Gaussian integrate nodes
    double *w_pt = gaulew48;                           //Gaussian integrate weight
    double  dn[npt]={0.0};

    for(int j=0; j<41; j++)
        for(int k=0; k<41; k++){
            double intp_eta=eta_lo + deta *j;    // -8 to 8 equal-distance
            double intp_phi=k*dphi;              //  0 to 2pi equal-distance
            double sumpt=0.0;

            for(int i=0; i<npt; i++){
                double intp_pt=ptCut_l + (ptCut_h-ptCut_l)*(1.0-p_pt[i]);
                dn[i] = dNdp3(intp_eta, intp_pt, intp_phi, M);
                sumpt += w_pt[i]*intp_pt*dn[i];
            }
            sumpt = sumpt*(ptCut_h-ptCut_l);
            N[j][k] = sumpt;
        }
}//////////////////////////////////////////////////////////////////////

void Init_coordinates(){
    //Calc pt[], phi[], eta[] 
    double *p_pt = gala15x;
    double *p_phi = gaulep48;
    for(int i=0; i<NPT; i++)pt[i]=p_pt[i]/slope;
    for(int j=0; j<NPHI/2;j++){
        phi[j] = pi*(1.0-p_phi[j]);
        phi[NPHI-1-j] = pi*(1.0+p_phi[j]);
    }
    for(int k=0; k<NETA; k++)eta[k]=eta_lo + k*deta;
}
//////////////////////////////////////////////////////////////////////

void Calc_Event_Plan2(const int &IEVENT, double dNdEtaPtdPtdPhi[NETA][NPT][NPHI], double psi_n[NH+1], double eps_n[NH+1]){
    double AptsinA[NH+1]={0.0};
    double AptcosA[NH+1]={0.0};

    for(int i=0; i<NETA; i++)
        for(int j=0; j<NPT; j++)
            for(int k=0; k<NPHI; k++){
                for(int n=1; n<=NH; n++){
                    AptsinA[n] += dNdEtaPtdPtdPhi[i][j][k]*pt[j]*sin(n*phi[k]);
                    AptcosA[n] += dNdEtaPtdPtdPhi[i][j][k]*pt[j]*cos(n*phi[k]);
                }
            }//End for

    for(int n=1; n<=NH; n++){
        psi_n[n] = atan(AptsinA[n]/max(1.0E-15,AptcosA[n]))/double(n);
    }//calc the event plane from emitted particles
}

//////////////////////////////////////////////////////////////////////
//
//
double lin_int(CD &x1, CD &x2, CD &y1, CD &y2, CD &x){
    /* linear interpolation  */
    assert(x1 != x2);
    double r=(x-x1)/(x2-x1);
    return y1*(1.0-r) + y2*r;
}
//////////////////////////////////////////////////////////////////////
double dNdp3(CD &eta1, CD &pt1, CD &phi1, double M[NETA][NPT][NPHI]){

    assert( phi1>=0.0 && phi1<=2.0*pi);
    int i=1;   // eta position
    int j=1;   // pt  position 
    int k=1;   // phi position

    while( eta1>eta[i] && i<NETA-1) i++;
    while( pt1>pt[j] && j<NPT-1) j++;
    while( phi1>phi[k] && k<NPHI-1) k++;

    double f1, f2, fy1, fy2;

    fy1 = lin_int(phi[k-1], phi[k], 
            M[i-1][j-1][k-1], M[i-1][j-1][k], phi1); 
    fy2 = lin_int(phi[k-1], phi[k], 
            M[i][j-1][k-1], M[i][j-1][k], phi1); 

    f1 = lin_int(eta[i-1], eta[i], fy1, fy2, eta1);         //interpolate for eta and phi at pt[k-1]

    fy1 = lin_int(phi[k-1], phi[k], 
            M[i-1][j][k-1], M[i-1][j][k], phi1); 
    fy2 = lin_int(phi[k-1], phi[k], 
            M[i][j][k-1], M[i][j][k], phi1); 

    f2 = lin_int(eta[i-1], eta[i], fy1, fy2, eta1);         //interpolate for eta and phi at pt[k]

    if(pt1 > 1.0){                                          //for pt>1.0, use log to linerize spectra
        f1 = (f1>1.0E-300)?f1:1.0E-300;
        f2 = (f2>1.0E-300)?f2:1.0E-300;
        f1 = log(f1); 
        f2 = log(f2);
    }

    double val = lin_int(pt[j-1], pt[j], f1, f2, pt1);     //interpolate for pt

    if(pt1 > 1.0) val = exp(val);

    val = (val>1.0E-300)?val:1.0E-300;

    return val;
}
//////////////////////////////////////////////////////////////////////

/*Calc Vn(pt) for -0.35<eta<0.35  */

void Calc_Vn1(CI &IEVENT, double dNdEtaPtdPtdPhi[NETA][NPT][NPHI], double Vn_Pt[NH+1][NPT])
{ /*Calc vn for |eta|<0.35   */
    double eps_n[NH+1]={0.0};                     //Initial geometric irregularities \eps_n , n=1,2,...NH
    double psi_n[NH+1]={0.0};                     //Event plan angle for each harmonic flow   n=1,2,...NH
    double *wgauss = gaulew48;
    //Calc_Event_Plan2(IEVENT, dNdEtaPtdPtdPhi, psi_n, eps_n);
    Calc_Event_Plan(IEVENT,  psi_n, eps_n);

    ecc2 += eps_n[2]/double(NEVENT);

    cout<<"Event "<<IEVENT<<":"<<endl;
    for(int i=1; i<=NH; i++){
        cout<<"eps["<<i<<"]="<<eps_n[i]<<"    ";
        cout<<"psi["<<i<<"]="<<psi_n[i]<<endl;
        cout<<endl;
    }

    for(int j=0; j!=NPT; j++){
        double dN=0.0;
        double Vn_dot_dN[NH+1]={0.0};

        for(int i=0; i!=NETA; i++){
            double dh =0.70/(NETA-1.0);
            double eta=-0.35+dh*i;             //integral over -0.35<eta<0.35
            //double dh =2.60/(NETA-1.0);
            //double eta=-1.3+dh*i;             //integral over -0.35<eta<0.35

            for(int k=0; k!=NPHI/2; k++){
                dN += wgauss[k]*(dNdp3(eta,pt[j],phi[k], dNdEtaPtdPtdPhi) + dNdp3(eta,pt[j],phi[NPHI-1-k], dNdEtaPtdPtdPhi));
                for(int n=1; n<=NH; n++){
                    Vn_dot_dN[n] += wgauss[k]*( cos( n*(phi[k]-psi_n[n]) )*dNdp3(eta,pt[j],phi[k], dNdEtaPtdPtdPhi) +
                            cos( n*(phi[NPHI-1-k]-psi_n[n]) ) * dNdp3(eta,pt[j],phi[NPHI-1-k], dNdEtaPtdPtdPhi));
                }
                // dN += wgauss[k]*(dNdEtaPtdPtdPhi[i][j][k]+dNdEtaPtdPtdPhi[i][j][NPHI-1-k]);
                // for(int n=1; n<=NH; n++){
                //     Vn_dot_dN[n] += wgauss[k]*( cos( n*(phi[k]-psi_n[n]) )*dNdEtaPtdPtdPhi[i][j][k]+
                //             cos( n*(phi[NPHI-1-k]-psi_n[n]) ) * dNdEtaPtdPtdPhi[i][j][NPHI-1-k]);
                // }
            }//End phi
        }//End Eta

        for(int n=1; n<=NH; n++){
            Vn_Pt[n][j] = Vn_dot_dN[n]/dN;
            //VN[n][j] += Vn_Pt[n][j]/double(NEVENT);
        }//Calc Vn
    }//End pt

    for(int n=1; n<=NH; n++){
        char fvn[128];
        sprintf(fvn,"C12/v%d_%d.dat",n,IEVENT);
        ofstream fout(fvn);
        for(int j=0; j!=NPT; j++){
            fout<<pt[j]<<' '<<Vn_Pt[n][j]<<endl;
        }
        fout.close();
    }

}
//////////////////////////////////////////////////////////////////////
void Calc_Vn(CI &IEVENT, double dNdEtaPtdPtdPhi[NETA][NPT][NPHI], double dN[NETA][NPT], double Vn_Eta_Pt[NH+1][NETA][NPT])
{
    double eps_n[NH+1]={0.0};                     //Initial geometric irregularities \eps_n , n=1,2,...NH
    double psi_n[NH+1]={0.0};                     //Event plan angle for each harmonic flow   n=1,2,...NH
    double *wgauss = gaulew48;
    Calc_Event_Plan(IEVENT, psi_n, eps_n);
    //Calc_Event_Plan2(IEVENT,dNdEtaPtdPtdPhi, psi_n, eps_n);

    //cout<<"Event 0:"<<endl;
    //for(int i=1; i<=NH; i++){
    //    cout<<"eps[n]="<<eps_n[i]<<endl;
    //    cout<<"psi[n]="<<psi_n[i]<<endl;
    //    cout<<endl;
    //}

    for(int i=0; i!=NETA; i++)
        for(int j=0; j!=NPT; j++){
            double dN_tmp=0.0;
            double Vn_dot_dN[NH+1]={0.0};

            for(int k=0; k!=NPHI/2; k++){
                dN_tmp += wgauss[k]*(dNdEtaPtdPtdPhi[i][j][k]+dNdEtaPtdPtdPhi[i][j][NPHI-1-k]);
                for(int n=1; n<=NH; n++){
                    Vn_dot_dN[n] += wgauss[k]*( cos( n*(phi[k]-psi_n[n]) )*dNdEtaPtdPtdPhi[i][j][k]+
                            cos( n*(phi[NPHI-1-k]-psi_n[n]) ) * dNdEtaPtdPtdPhi[i][j][NPHI-1-k]);
                }
            }//End phi

            for(int n=1; n<=NH; n++){
                //Vn_Eta_Pt[n][i][j] = Vn_dot_dN[n]/dN_tmp;
                Vn_Eta_Pt[n][i][j] = Vn_dot_dN[n];
            }//Calc Vn
            dN[i][j] = dN_tmp; 
        }//End pt

    /*    for(int n=1; n<=NH; n++){
          char fvn[128];
          sprintf(fvn,"C12/v%d_%d.dat",n,IEVENT);
          ofstream fout(fvn);
          for(int j=0; j!=NPT; j++){
          fout<<pt[j]<<' '<<Vn_Eta_Pt[n][20][j]<<endl;
          }
          fout.close();
          }//Not eta range integral
     */
}
//////////////////////////////////////////////////////////////////////
void Calc_Vn_ptcut(double dN[NETA][NPT], double Vn_Eta_Pt[NH+1][NETA][NPT], double Vn_trig[NH+1][NETA],CD & ptCut_l,CD & ptCut_h)
{
    const int npt=24;
    double *p_pt = gaulep48;                           //Gaussian integrate nodes
    double *w_pt = gaulew48;                           //Gaussian integrate weight
    for(int n=1; n<=NH; n++){
        double  vn[npt]={0.0};
        double  dn[npt]={0.0};

        for(int j=0; j<NETA; j++){

            double sumpt1=0.0;  //sum over vn*dN/dPt
            double sumpt2=0.0;  //sum over dN/dPt

            for(int i=0; i<npt; i++){
                double pt1=ptCut_l + (ptCut_h-ptCut_l)*(1.0-p_pt[i]);
                int k=1;
                while( pt1>pt[k] && k<NPT-1) k++;

                vn[i] = lin_int(pt[k-1], pt[k], Vn_Eta_Pt[n][j][k-1], Vn_Eta_Pt[n][j][k], pt1);
                sumpt1 += w_pt[i]*pt1*vn[i];

                dn[i] = lin_int(pt[k-1], pt[k], dN[j][k-1], dN[j][k], pt1);
                sumpt2 += w_pt[i]*pt1*dn[i];

            }//End pt
            sumpt1 = sumpt1*(ptCut_h-ptCut_l);
            sumpt2 = sumpt2*(ptCut_h-ptCut_l);
            Vn_trig[n][j] = sumpt1/sumpt2;
        }//End eta
    }//End vn
}
//////////////////////////////////////////////////////////////////////
void Generate_Trig_Assos(int IEVENT, double Ntrig[41][41], double Nasso[41][41],\
        double Vn_trig[NH+1][NETA], double Vn_asso[NH+1][NETA])
{
    double dNdEtaPtdPtdPhi[NETA][NPT][NPHI];               //spectra

    double Vn_Eta_Pt[NH+1][NETA][NPT];                     //Vn(eta, pt)

    double dN[NETA][NPT];

    double Vn_Pt[NH+1][NPT];                     //Vn(eta, pt)

    Init_coordinates();                                    //Initialize eta, pt, phi from hydro output

    ReadSpectra(IEVENT,dNdEtaPtdPtdPhi);                   //

    dNdEtadPhi(dNdEtaPtdPtdPhi, Ntrig, trig_lo, trig_hi);          //Calc Ntrig
    dNdEtadPhi(dNdEtaPtdPtdPhi, Nasso, asso_lo, asso_hi);          //Calc Nasso

    Calc_Vn(IEVENT, dNdEtaPtdPtdPhi, dN, Vn_Eta_Pt);                   //Calc Vn(eta, pt)
    Calc_Vn1(IEVENT, dNdEtaPtdPtdPhi, Vn_Pt);                   //Calc Vn(eta, pt)

    Calc_Vn_ptcut(dN, Vn_Eta_Pt, Vn_trig, trig_lo, trig_hi);           //Calc Vn_trig(eta)
    Calc_Vn_ptcut(dN, Vn_Eta_Pt, Vn_asso, asso_lo, asso_hi);           //Calc Vn_asso(eta)

    for(int n=1; n<=NH; n++){
        for(int j=0; j!=NPT; j++){
            //VN[n][j] += Vn_Pt[n][j]*Vn_Pt[n][j]/double(NEVENT);
            VN[n][j] += Vn_Pt[n][j]/double(NEVENT);
        }
    }

}

//////////////////////////////////////////////////////////////////////

CI NI=41;          //for eta
CI NJ=41;          //for phi

int main(int argc, char **argv){



    double C12[NI][NJ]={0.0};
    double AN1N2[NI][NJ]={0.0};
    double AN1AN2[NI][NJ]={0.0};
    double AN1[NI][NJ]={0.0};  //<Ntrig>
    double AN2[NI][NJ]={0.0};  //<Nasso>

    double f[NI][NJ]={0.0};    
    /* f(deta, dphi) = B(1+sum_{n=1}^{6}
     *<vn_trig(eta1)vn_asso(eta2) + vn_trig(eta2)vn_asso(eta1)>*cos(n*dphi))
     where deta=(eta1-eta2);
     < > is event average;
     */

    double Ntrig[NI][NJ];
    double Nasso[NI][NJ];

    double V2trig[NH+1][NI];
    double V2asso[NH+1][NI];

    for(int n=0; n<NEVENT; n++){
        Generate_Trig_Assos(n, Ntrig, Nasso, V2trig, V2asso);

        for(int dEta=0; dEta<NI; dEta++)
            for(int dPhi=0; dPhi<NJ-1; dPhi++){
                f[dEta][dPhi] += 1.0;
                int NRE = NI-dEta; //Record how many times does dEta repeat in one event   
                for(int Eta1=0; Eta1<NRE; Eta1++){
                    int Eta2 = Eta1 + dEta;
                    for(int i=2; i<=NH; i++){ //sum over harmonic flow
                        f[dEta][dPhi] += (V2trig[i][Eta1]*V2asso[i][Eta2] + V2trig[i][Eta2]*V2asso[i][Eta1]) \
                                         *cos(i*dPhi*dphi)/double(NRE);
                    }//End NH
                }//End Eta1
                /*        for(int i=2; i<=NH; i++){ //sum over harmonic flow
                          f[dEta][dPhi] += 2.0*V2trig[i][NI/2]*V2asso[i][NI/2]*cos(i*dPhi*dphi);
                          }
                 */
                f[dEta][dPhi] /= double(NEVENT);
            }//End dPhi




        cout<<"V2trig(eta=0)="<<V2trig[2][20]<<" ";
        cout<<"V2asso(eta=0)="<<V2asso[2][20]<<endl;

        cout<<"V3trig(eta=0)="<<V2trig[3][20]<<" ";
        cout<<"V3asso(eta=0)="<<V2asso[3][20]<<endl;



        for(int i=0; i!=NI; i++)
            for(int j=0; j!=NJ; j++){
                AN1[i][j] += Ntrig[i][j]/double(NEVENT);
                AN2[i][j] += Nasso[i][j]/double(NEVENT);
            }

        for(int dEta=0; dEta<NI; dEta++)
            for(int dPhi=0; dPhi<NJ-1; dPhi++){
                int NRE = NI-dEta; //Record how many times does dEta repeat in one event   
                int NRP = NJ-1;      //
                for(int Eta1=0; Eta1<NRE; Eta1++)
                    for(int Phi1=0; Phi1<NRP; Phi1++){
                        int Eta2 = Eta1 + dEta;
                        int Phi2 = (Phi1+dPhi)%(NJ-1);
                        AN1N2[dEta][dPhi] += 0.25*Ntrig[Eta1][Phi1]*Nasso[Eta2][Phi2]/double(NRE*NRP*NEVENT);
                        AN1N2[dEta][dPhi] += 0.25*Ntrig[Eta1][Phi2]*Nasso[Eta2][Phi1]/double(NRE*NRP*NEVENT);
                        AN1N2[dEta][dPhi] += 0.25*Nasso[Eta1][Phi2]*Ntrig[Eta2][Phi1]/double(NRE*NRP*NEVENT);
                        AN1N2[dEta][dPhi] += 0.25*Nasso[Eta1][Phi1]*Ntrig[Eta2][Phi2]/double(NRE*NRP*NEVENT);
                    }
            }
        cout<<"n="<<n<<endl;

    } //End for NEVENT, Calc <N1N2>


    for(int dEta=0; dEta<NI; dEta++)
        for(int dPhi=0; dPhi<NJ-1; dPhi++){

            int NRE = NI-dEta; //Record how many times does dEta repeat in one event   
            int NRP = NJ-1;      //
            for(int Eta1=0; Eta1<NRE; Eta1++)
                for(int Phi1=0; Phi1<NRP; Phi1++){
                    int Eta2 = Eta1 + dEta;
                    int Phi2 = (Phi1+dPhi)%(NJ-1);
                    AN1AN2[dEta][dPhi] += 0.25*AN1[Eta1][Phi1]*AN2[Eta2][Phi2]/double(NRE*NRP);
                    AN1AN2[dEta][dPhi] += 0.25*AN1[Eta1][Phi2]*AN2[Eta2][Phi1]/double(NRE*NRP);
                    AN1AN2[dEta][dPhi] += 0.25*AN2[Eta1][Phi2]*AN1[Eta2][Phi1]/double(NRE*NRP);
                    AN1AN2[dEta][dPhi] += 0.25*AN2[Eta1][Phi1]*AN1[Eta2][Phi2]/double(NRE*NRP);
                }
        } // <N1><N2>



    double fmin=100.0;
    int fmin_Eta=100;
    int fmin_Phi=100;

    double C12_min=100000.0;
    int C12m_Eta=100;
    int C12m_Phi=100;


    for(int dEta=0; dEta<NI; dEta++)
        for(int dPhi=0; dPhi<NJ-1; dPhi++){
            C12[dEta][dPhi] = AN1N2[dEta][dPhi]/AN1AN2[dEta][dPhi];
        }

    //for(int dEta=0; dEta<NI; dEta++)
    int dEta=0;
    for(int dPhi=0; dPhi<NJ-1; dPhi++){
        if(f[dEta][dPhi]<fmin){
            fmin=f[dEta][dPhi];
            fmin_Eta=dEta;
            fmin_Phi=dPhi;
        }

        if(C12[dEta][dPhi]<C12_min){
            C12_min=C12[dEta][dPhi];
            C12m_Eta=dEta;
            C12m_Phi=dPhi;
        }
    }//Calc using ZYAM 
    //double B=C12_min/fmin;
    double B=(C12_min-1.0)/f[C12m_Eta][C12m_Phi];

    cout<<"fmin("<<fmin_Eta<<","<<fmin_Phi<<")="<<fmin<<endl;
    cout<<"C12_min("<<C12m_Eta<<","<<C12m_Phi<<")="<<C12_min<<endl;



    ofstream fout("C12/C12_bin0.dat");
    fout.precision(10);
    for(int dEta=0; dEta<NI; dEta++)
        for(int dPhi=0; dPhi<NJ-1; dPhi++){
            //assert(AN1AN2[dEta][dPhi]!=0.0);
            //C12[dEta][dPhi] = AN1N2[dEta][dPhi]/AN1AN2[dEta][dPhi] - 1.0;
            //C12[dEta][dPhi] =AN1N2[dEta][dPhi]-AN1AN2[dEta][dPhi];
            double DeltaEta = dEta*deta;
            double DeltaPhi = dPhi*dphi;
            double ANA=AN1AN2[dEta][dPhi];
            double c12 = (ANA>0.0)?(AN1N2[dEta][dPhi]/ANA-1.0):-1.0;
            fout<<C12[dEta][dPhi]<<' '<<c12<<' '<<f[dEta][dPhi]<<endl;
        }
    fout.close();



    //cout<<"dNdp3(0,0,0,dNdEtaPtdPtdPhi)="<<dNdp3(0.0,0.0077,0,dNdEtaPtdPtdPhi)<<endl;
    cout<<"Ntrig[0][0]="<<Ntrig[10][20]<<endl;
    cout<<"Ntrig[3][0]="<<Ntrig[20][20]<<endl;
    cout<<"Ntrig[5][0]="<<Ntrig[30][20]<<endl;
    cout<<"Nasso[0][0]="<<Nasso[10][0]<<endl;
    cout<<"Nasso[3][0]="<<Nasso[20][0]<<endl;
    cout<<"Nasso[5][0]="<<Nasso[30][0]<<endl;




    for(int n=1; n<=NH; n++){
        char fvn[128];
        sprintf(fvn,"C12/v%d.dat",n);
        ofstream fout(fvn);
        for(int j=0; j!=NPT; j++){
            //fout<<pt[j]<<' '<<sqrt(VN[n][j])<<endl;
            fout<<pt[j]<<' '<<VN[n][j]<<endl;
        }
        fout.close();
    }

    cout<<"ecc2="<<ecc2<<endl;

    return 0;



}
