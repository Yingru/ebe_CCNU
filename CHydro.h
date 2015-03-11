#ifndef __CHYDROH
#define __CHYDROH

#include<iomanip>
#include<fstream>
#include<string>
#include<deque>
#include<cmath>
#include<cassert>
#include<cstdarg>
#include<cstdlib>
#include<ctime>
#include<sys/types.h>
#include<sys/stat.h>
#include<sys/wait.h>

#include"Algorithm.h"
#include"CEos.h"
#include"CMatrix.h"
#include"CAAcollision.h"
#include"Int.h"
#include"CCube.h"
#include"SolveU_mu.h"
typedef ofstream cof;
using namespace std;

const int NPT=15;
const int NPHI=48;
const int NETA=41;
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



class CHydro3D{
    public:
	HydroConfig setting;             //Input parameters, user can change the input parameters in Vishydro.inp
	CEos Eos;                        //Equation Of State, Eos.P(), Eos.T(), Eos.S()
	int IEVENT;                   //For event by event calculation
	int NX0, NY0, NZ0;            //Lower boundary of numerical grid
	int NX, NY, NZ;               //Uper boundary of numerical grid 
	double tau;                   //Initial time plus evolution time
	double t0, dt, dx, dy, dz;     //Initial time and Grid size in t,x,y,z direction
	double ED_MAX;                 //if ED_MAX<e_dec, the hydro stops

	M3D Ed ;      //Energy density Ed(i,j,k)
	M3D Bd ;      //Baryon density Bd(i,j,k)
	M3D Pr ;      //Pressure density Pr()
	M3D Vx ;      //Velocity Vx
	M3D Vy ;      //Vy
	M3D Vz ;      //Veta

	M3D TT00;      //tau*T^{\tau\tau}(i,j,k)
	M3D TT01;      //tau*T^{\tau x  }
	M3D TT02;      //tau*T^{\tau y  }
	M3D TT03;      //tau*T^{\tau\eta}
	M3D TJT;       //tau*\J^\tau= \tau*nb*\gamma
	M3D TJAT;       //tau*\J^\tau= \tau*nb*\gamma

	deque<CParticle> particles;
	deque<SF> hypsf;

	CHydro3D();
	CHydro3D(int Xl, int Xh, int Yl, int Yh, int Zl, int Zh);
	void Initialization(char* hydro_config);
	void PPInitial();
	void PPInitial2();
	void PPInitial3();
	void PPInitial4();
	void Evolution();
	void Freezeout();
	void Freezeout_xyz_half();
	void Freezeout2();
	void Freezeout3();
	void Freezeout4();

	void FrzIsoTime();//Iso-time freeze out for pp collisions
	void ReadParticles();//Read particles' information
	void dNdYdptdphi(CParticle par);//calc dN/dYPtdPtdPhi for particle i
	void dNdEta(CI & DIR_OR_RESO, CParticle & par);//calc dNdEta for particle i; DIR_OR_RESO=0 for before resonance decay, DIR_OR_RESO=1 for after resonance decay
	void dNdEta_Charged();
	//void dNdPt(CParticle par, CD &DY);//calc dNdPt in |Y|<DY rapidity window
	//void dNdPhi(CParticle par, CD &DY);//calc dNdPhi in |Y|<DY rapidity window
	//void Vn(CI &n, CParticle par, CD &DY);//calc Vn(n=2,3,4,...) in |Y|<DY rapidity window for particle i

	void Output(CI & NT);
	void Rootfinding(double &E0, double &B0, CD &T00IJ, CD &T01IJ, CD &T02IJ, CD &T03IJ, CD &JTIJ);
	void CalcSources(M3D &ScT0, M3D &ScT1, M3D &ScT2, M3D &ScT3, M3D &ScJT, M3D &ScJA);
	void AnomalySources(M3D &ScT0, M3D &ScT1, M3D &ScT2, M3D &ScT3, M3D &ScJT, M3D &ScJA);

	void CalcHyperSf(M3D &Ed0,  M3D &Vx0, M3D &Vy0, M3D &Vz0);

	void CalcHyperSf2(CI &NT, M3D &Ed0,  M3D &Vx0, M3D &Vy0, M3D &Vz0);
	void StepUpdate();
	void StepUpdate2(M3D & H00, M3D & H01, M3D & H02, M3D & H03, M3D & HJT,
                     M3D & F00, M3D & F01, M3D & F02, M3D & F03, M3D & FJT,
                     CD & DT); //using the half time step ed and v
	void CalcSources2(M3D &ScT0, M3D &ScT1, M3D &ScT2, M3D &ScT3, M3D &ScJT, M3D &ScJA,
                      M3D &H00, M3D &H03, CD &DT);
	void SetBoundary(M3D &M0);//set the boundary value of Ed, Vx, Vy ...
	void VzSmooth();//
	void Accecentricity(double &epsx, double &epsp);
	inline double Gamma(CI &i, CI &j, CI &k);//gamma=1.0/sqrt(1-vx^2-vy^2-tau^2*Veta^2)

	inline double STotal();//Calc stotal at time tau to check entropy conservation

};

#endif
