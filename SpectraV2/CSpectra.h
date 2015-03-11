#ifndef __CSPECTRA__
#define __CSPECTRA__

#include<iostream>
#include<fstream>
#include<string>
#include<deque>
#include<cstdlib>
#include<cmath>
#include<cassert>
#include"../Int.h"
#include"../PhyConst.h"

using namespace std;
typedef const int CI;
typedef const double CD;

CI  NPT=15;
CI  NPHI=48;
CI  NETA=41;
struct SF{
    double DA0,DA1,DA2,DA3,vx,vy,vz,etas;
};//The freeze out hypersurf structure

struct CParticle{
    int monval;
    string  name;
	double	mass;
	double	width;
	int	gspin;      /* spin degeneracy */
	int	baryon;     /* baryon number   */
	int	strange;
	int	charm;
	int	bottom;
	int	gisospin;  /* isospin degeneracy */
	int	charge;    
	int	decays;    /* amount of decays listed for this resonance */
	int	stable;     /* defines whether this particle is considered as stable */

    double PT[NPT];
    double PHI[NPHI];
    double ETA[NETA];
    double dNdPt[NPT];
    double dNdPhi[NPHI];
    double dNdEta[NETA];
};

class CSpectra
{
    private:
        int NEVENT;              /*Event Number */
    	deque<SF> hypsf;         /*Frz out hyper surface*/
        void Calc_dNdYPtdPtdPhi(CI &IEVENT, CParticle &par);/*Before Reso Decay*/
        void dNdEta(CI& IEVENT, CI &DIR_OR_RESO, CParticle &par);

    public:
    	deque<CParticle> particles; /*Particles data table*/
        double frz_temp; /*  Frz out temperature */
        CSpectra();      
        CSpectra(CI & NEVENT);               
        void ReadParticles(char * particle_data_table);  /*Read particle data table */

        void Set_FrzT(CD &TFRZ);  /* Set Frz temperature */
        void ReadHyperSF(CI & IEVENT);           /*Read Frz out Hyper surface data */
        void CalcSpectra_BeforeReso(CI & IEVENT);/*Calc dN/dYPtdPtdPhi before resonance decay */
        void CalcResonanceDecay(CI & IEVENT);    /*Calc resonance decay---not implemented in this class yet */
        void CalcSpectra_AfterReso(CI & IEVENT); /*dNdEta, dNdPt, v2 for identified particles after Resonance decay*/
        void dNdEta_Charged(CI & IEVENT);        /*dNdEta, dNdPt, v2 for charged particles after Resonance decay*/
};

#endif
