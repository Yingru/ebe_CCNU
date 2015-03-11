#include<iostream>
#include<cstdio>
#include<cassert>
#include<fstream>
const int NEVENT=20;
const int NETA=41;
const int NPHI=41;
const double pi=3.141592653;
const double etahi=6.0;
const double etalo=-6.0;
const double deta=(etahi-etalo)/(NETA-1.0);
const double dphi=2.0*pi/(NPHI-1.0);


using namespace std;

int main(int argc, char** argv)
{
    assert(argc==2);
    int IE=2*atoi(argv[1]);

    double C12[NETA][NPHI]={0.0};
    double AN1N2[NETA][NPHI]={0.0};     // <N1N2>[\Delat ETA][\Delta PHI]
    double AN1AN2[NETA][NPHI]={0.0};     // <N1><N2>[\Delat ETA][\Delta PHI]
    double AN[NETA][NPHI]={0.0};        // <N>

    double N[NEVENT][NETA][NPHI]={0.0}; //  record N for each event and each bin
    double dN[NEVENT][NETA][NPHI];      //  dN = N-<N>  is the fluctuation in each bin
    double AdN1dN2[NETA][NPHI]={0.0};   //  <dN1dN2> 

    for(int i=0; i<NETA; i++)
        for(int j=0; j<NPHI; j++){
            AN[i][j]=0.0;
            AdN1dN2[i][j]=0.0;
        }

    for(int n=1; n<=NEVENT; n++){
//    for(int n=IE; n<=IE+1; n++){
        char fdn[50];
        sprintf(fdn, "results/dNdEtadPhi%d.dat",n);
        int IM=IE+n-1;
        //sprintf(fdn, "results/dNdEtadPhi%d.dat",IM);
        ifstream fin(fdn);
        for(int i=0; i<NETA; i++)
            for(int j=0; j<NPHI-1; j++){
                fin>>N[n-1][i][j];          //Read dNdEtadPhi in one evevent
                AN[i][j] += N[n-1][i][j]/double(NEVENT);
            }
        for(int dEta=0; dEta<NETA; dEta++)
            for(int dPhi=0; dPhi<NPHI-1; dPhi++){
                int NRE = NETA-dEta; //Record how many times does dEta repeat in one event   
                int NRP = NPHI-1;      //
                for(int Eta1=0; Eta1<NRE; Eta1++)
                    for(int Phi1=0; Phi1<NRP; Phi1++){
                        int Eta2 = Eta1 + dEta;
                        int Phi2 = (Phi1+dPhi)%(NPHI-1);
                        AN1N2[dEta][dPhi] += 0.5*N[n-1][Eta1][Phi1]*N[n-1][Eta2][Phi2]/double(NRE*NRP*NEVENT);
                        AN1N2[dEta][dPhi] += 0.5*N[n-1][Eta1][Phi2]*N[n-1][Eta2][Phi1]/double(NRE*NRP*NEVENT);
                    }
            }
    } //End for NEVENT, Calc <N1N2>

    for(int n=1; n<=NEVENT; n++){
        for(int i=0; i<NETA; i++)
            for(int j=0; j<NPHI-1; j++){
                dN[n-1][i][j] = N[n-1][i][j] - AN[i][j];
            }
        for(int dEta=0; dEta<NETA; dEta++)
            for(int dPhi=0; dPhi<NPHI-1; dPhi++){
                int NRE = NETA-dEta; //Record how many times does dEta repeat in one event   
                int NRP = NPHI-1;
                for(int Eta1=0; Eta1<NRE; Eta1++)
                    for(int Phi1=0; Phi1<NRP; Phi1++){
                        int Eta2 = Eta1 + dEta;
                        int Phi2 = (Phi1 + dPhi)%(NPHI-1);
                        AdN1dN2[dEta][dPhi] += 0.5*dN[n-1][Eta1][Phi1]*dN[n-1][Eta2][Phi2]/double(NRE*NRP*NEVENT);
                        AdN1dN2[dEta][dPhi] += 0.5*dN[n-1][Eta1][Phi2]*dN[n-1][Eta2][Phi1]/double(NRE*NRP*NEVENT);
                    }
            }
    }//End for NEVENT, Calculate dN and <dN1dN2>


    char fc12[50];
    int IR=atoi(argv[1]);
    //sprintf(fc12,"C12/C12_bin%d.dat",IR);
    //ofstream fout(fc12);
    ofstream fout("C12/C12_bin0.dat");
    fout.precision(10);
    for(int dEta=0; dEta<NETA; dEta++)
        for(int dPhi=0; dPhi<NPHI-1; dPhi++){
            int NRE = NETA-dEta; //Record how many times does dEta repeat in one event   
            int NRP = NPHI-1;      //Record how many times does dPhi repeat in one event   
            for(int Eta1=0; Eta1<NRE; Eta1++)
                for(int Phi1=0; Phi1<NRP; Phi1++){
                    int Eta2 = Eta1 + dEta;
                    int Phi2 = (Phi1 + dPhi)%(NPHI-1);
                    AN1AN2[dEta][dPhi] += 0.5*AN[Eta1][Phi1]*AN[Eta2][Phi2]/double(NRE*NRP);
                    AN1AN2[dEta][dPhi] += 0.5*AN[Eta1][Phi2]*AN[Eta2][Phi1]/double(NRE*NRP);
                }
        }//Calce <N1><N2>

    for(int dEta=0; dEta<NETA; dEta++)
        for(int dPhi=0; dPhi<NPHI-1; dPhi++){
            //assert(AN1AN2[dEta][dPhi]!=0.0);
            //C12[dEta][dPhi] = AN1N2[dEta][dPhi]/AN1AN2[dEta][dPhi] - 1.0;
            //C12[dEta][dPhi] =AN1N2[dEta][dPhi]-AN1AN2[dEta][dPhi];
            double DeltaEta = dEta*deta;
            double DeltaPhi = dPhi*dphi;
            double ANA=AN1AN2[dEta][dPhi];
            double c12 = (ANA>0.0)?(AN1N2[dEta][dPhi]/ANA-1.0):-1.0;
            fout<<AN1N2[dEta][dPhi]<<' '<<ANA<<' '<<c12<<endl;
        }
    fout.close();
}



