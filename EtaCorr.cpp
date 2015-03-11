#include<iostream>
#include<cstdio>
#include<cassert>
#include<fstream>
const int NEVENT=100;
const int NETA=41;
const double etahi=6.0;
const double etalo=-6.0;
const double deta=(etahi-etalo)/(NETA-1.0);

using namespace std;

int main(int argc, char** argv)
{
    double C12[NETA]={0.0};
    double AN1N2[NETA]={0.0};     // <N1N2>[\Delat ETA][\Delta PHI]
    double AN1AN2[NETA]={0.0};     // <N1><N2>[\Delat ETA][\Delta PHI]
    double AN[NETA]={0.0};        // <N>

    double N[NETA]={0.0}; //  record N for each event and each bin

    for(int i=0; i<NETA; i++){
            AN[i]=0.0;
        }

    for(int n=1; n<=NEVENT; n++){
        char fdn[50];
        sprintf(fdn, "results/dNdEta%d.dat",n);
        ifstream fin(fdn);
        double eta;
        for(int i=0; i<NETA; i++){
                fin>>eta>>N[i];          //Read dNdEtadPhi in one evevent
                AN[i] += N[i]/double(NEVENT);
            }
        for(int dEta=0; dEta<NETA; dEta++){
                int NRE = NETA-dEta; //Record how many times does dEta repeat in one event   
                for(int Eta1=0; Eta1<NRE; Eta1++){
                        int Eta2 = Eta1 + dEta;
                        AN1N2[dEta] += N[Eta1]*N[Eta2]/double(NRE*NEVENT);
                    }
            }
    } //End for NEVENT, Calc <N1N2>


    ofstream fout("C12_eta.dat");
    fout.precision(10);
    for(int dEta=0; dEta<NETA; dEta++){
            int NRE = NETA-dEta; //Record how many times does dEta repeat in one event   
            for(int Eta1=0; Eta1<NRE; Eta1++){
                    int Eta2 = Eta1 + dEta;
                    AN1AN2[dEta] += AN[Eta1]*AN[Eta2]/double(NRE);
                }
        }//Calce <N1><N2>

    for(int dEta=0; dEta<NETA; dEta++){
            //assert(AN1AN2[dEta][dPhi]!=0.0);
            //C12[dEta][dPhi] = AN1N2[dEta][dPhi]/AN1AN2[dEta][dPhi] - 1.0;
            //C12[dEta][dPhi] =AN1N2[dEta][dPhi]-AN1AN2[dEta][dPhi];
            double DeltaEta = dEta*deta;
            double ANA=AN1AN2[dEta];
            double c12 = (ANA>0.0)?(AN1N2[dEta]/ANA-1.0):-1.0;
            //fout<<AN1N2[dEta]<<' '<<ANA<<' '<<c12<<endl;
            fout<<DeltaEta<<' '<<c12<<endl;
        }
    fout.close();
}



